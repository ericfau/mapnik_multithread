/*****************************************************************************
 *
 * Copyright (C) 2011 Andrew Harvey <andrew.harvey4@gmail.com>
 * 
 * Additional work:
 * Copyright (C) 2013 Lucas Madar <lucas.madar@gmail.com>
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 *****************************************************************************/
 
 
/* This program takes a list of OSM meta-tiles, and renders them using Mapnik,
   then chops up the meta-tiles into regular 256x256 size tiles. */
 
 
// IO libraries
#include <iostream>
#include <fstream>
#include <string>
#include <set>
 
#include <math.h>
 
#include <assert.h>
 
// queue for keeping store of which tiles to render
#include <queue>
 
// filesystem libraries for making new directories
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
 
// use Boost to parse program arguments/options
#include <boost/program_options.hpp>
namespace po = boost::program_options;
 
// use pthreads to multithread the rendering engine
#include <pthread.h>
 
// Mapnik includes
#include <mapnik/map.hpp>
#include <mapnik/datasource_cache.hpp>
#include <mapnik/load_map.hpp>
#include <mapnik/font_engine_freetype.hpp>
#include <mapnik/agg_renderer.hpp>
#include <mapnik/image_util.hpp>
#include <mapnik/config_error.hpp>
 
// ImageMagick for splitting meta tiles
#include <Magick++.h>
 
/* Default width and height of each individual tile in pixels.
 
   This value is defined by the specification at 
   http://wiki.openstreetmap.org/wiki/Slippy_map_tilenames */
#define TILE_SIZE 256
 
/* limits of the Spherical Mercator projection that OSM uses on the default web
   map */
#define MAX_X 20037508
#define MAX_Y 20037508
 
/* A Meta-OSM Slippy map tile structure */
struct Metatile {
    Metatile( const unsigned int z, const unsigned int x, const unsigned int y ) : z(z), x(x), y(y) {}
 
    const unsigned int z, x, y;
};
 
struct Metatile_ptr_comp {
    bool operator()( const Metatile* lhs, const Metatile* rhs ) const {
        int diff;
 
        if( lhs->z != rhs->z ) { diff = lhs->z - rhs->z; }
        else if( lhs->x != rhs->x ) { diff = lhs->x - rhs->x; }
        else { diff = lhs->y - rhs->y; };
 
        return diff < 0;
    }
};
 
/* An arbitary geometry bounding box with no coordinate system */
struct Bbox {
    double left, bottom, right, top;
};
 
// function declarations

void removeDupTiles(std::set<Metatile*,Metatile_ptr_comp> &tileSet);

Bbox tileToMercBounds(int z, int x, int y);


void *renderThread(void *argument);
Metatile *getNextTile();
 
// empty global tile queue
std::deque<Metatile*> tileQueue;
pthread_mutex_t tileQueueMutex = PTHREAD_MUTEX_INITIALIZER;
 
// The buffer in pixels to render beyond the extent of the meta tile.
int buffer = 128;
 
// mapnik map.xml file
std::string mapfile;
 
/* metaTiles^2 will equal the number of tiles to render at a time as one large
   meta tile, which will then be split up and save separately */
unsigned int metaTiles = 8;
 
bool verbose = false;
bool png256 = false;
 
std::string tileFileName;
std::string outputDirName;
 
int main ( int argc , char** argv )
{
    try {
 
        /*     configuration     */
 
        // number of rendering threads
        int numThreads = 1;
 
        // read program options
        po::options_description desc("Allowed options");
        desc.add_options()
            ("help,h", "produce help message")
            ("png256,2", "use 256 color png palette")
            ("verbose,v", "show verbose output")
            ("buffer,b", po::value<int>(&buffer)->default_value(128), "set rendering buffer")
            ("metaTiles,t", po::value<unsigned int>(&metaTiles)->default_value(8), "set number of metaTiles")
            ("mapfile,m", po::value(&mapfile), "path to mapnik style file (osm.xml)")
            ("threads,j", po::value<int>(&numThreads)->default_value(4), "set number of rendering threads")
            ("tileFile,f", po::value(&tileFileName), "file with list of tiles to render")
            ("outputdir,o", po::value(&outputDirName)->default_value("tiles"), "output directory")
        ;
 
        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);
 
        if (vm.count("help")) {
            std::cout << desc << "\n";
            return EXIT_SUCCESS;
        }
 
        if (vm.count("png256")) {
            png256 = true;
        }
 
        if (!vm.count("mapfile")) {
            std::cout << desc << "\n";
            return EXIT_FAILURE;
        }
 
        png256 = vm.count("png256") > 0;
        verbose = vm.count("verbose") > 0;
 
        // read all tiles into tile set
        std::set<Metatile*,Metatile_ptr_comp> tileSet;
        std::string line;
        std::ifstream tileFile(tileFileName.c_str());
        if (tileFile.is_open()) {
            while (tileFile.good()) {
                getline(tileFile, line);
 
                // skip empty lines
                if (line.size() > 0) {
                    // parse the line for an z/x/y tile reference
                    int slash1 = line.find('/');
                    int z = atoi(line.substr(0, slash1).c_str());
                    int slash2 = line.find('/', slash1 + 1);
                    int x = atoi(line.substr(slash1 + 1, slash2 - slash1 - 1).c_str());
                    int slash3 = line.find('/', slash2 + 1);
                    int y = atoi(line.substr(slash2 + 1, slash3 - slash2 - 1).c_str());
 
                    // check z/x/y tile looks like a valid tile
                    if (!((z >= 0) && (z <= 30))) {
                        fprintf(stderr, "Unexpected tile reference.Does not satisfy 0 <= z <= 30, found: %d", z);
                        continue;
                    }
 
                    if (!((x >= 0) && (x < pow(2, z)))) {
                        fprintf(stderr, "Unexpected tile reference. Does not satisfy 0 <= x < 2^%d, found: %d", z, x);
                        continue;
                    }
 
                    if (!((y >= 0) && (y < pow(2, z)))) {
                        fprintf(stderr, "Unexpected tile reference. Does not satisfy 0 <= y < 2^%d, found: %d", z, y);
                        continue;
                    }
 
                    Metatile *tile = new Metatile(z, x, y);
 
                    std::pair<std::set<Metatile*,Metatile_ptr_comp>::iterator,bool> ret = tileSet.insert(tile);
 
                    if(ret.second==false) {  					
                        fprintf(stderr, "Duplicate tile in tileSet: %d/%d/%d\n", z, x, y); 
                    }
                }
            }
        }else{
            std::cout << "Unable to open file.\n";
            return EXIT_FAILURE;
        }
 
        std::cout << "Loaded " << tileSet.size() << " tiles.\n";
        removeDupTiles(tileSet);
 
        std::copy(tileSet.begin(), tileSet.end(), std::back_inserter(tileQueue));
 
        // register mapnik datasources
        mapnik::datasource_cache::instance()->register_datasources("/usr/lib/mapnik/input/");
 
        // register all truetype fonts recursively below the Debian font directory
        mapnik::freetype_engine::register_fonts("/usr/share/fonts/truetype/", true);
 
        // span rendering threads
        pthread_t renderingThreads[numThreads];
        int threadID[numThreads];
 
        for (int i = 0; i < numThreads; i++) {
            threadID[i] = i;
            if(verbose)
                printf("Spawning rendering thread %d\n", i);
            int rc = pthread_create(&renderingThreads[i], NULL, renderThread, (void *) &threadID[i]);
            assert(rc == 0);
        }
 
        /* wait for all threads to complete */
        for (int i = 0; i < numThreads; i++) {
            int rc = pthread_join(renderingThreads[i], NULL);
            assert(0 == rc);
        }
 
        std::cout << "done";
        exit(EXIT_SUCCESS);
    }
    catch ( const mapnik::config_error & ex )
    {
        std::cerr << "### Configuration error: " << ex.what() << std::endl;
        return EXIT_FAILURE;
    }
    catch ( const std::exception & ex )
    {
        std::cerr << "### std::exception: " << ex.what() << std::endl;
        return EXIT_FAILURE;
    }
    catch ( ... )
    {
        std::cerr << "### Unknown exception." << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
 
/* removes tiles from the set that would be duplicated by other tiles' metatiles */
void removeDupTiles(std::set<Metatile*,Metatile_ptr_comp> &tileSet) {
    std::set<Metatile*,Metatile_ptr_comp>::iterator it;;
 
    unsigned int lastZoomLevel = 0xFFFFFFFF;
    unsigned int tilesRemoved = 0;
    unsigned int tilesRemovedTotal = 0;
 
    for(it = tileSet.begin(); it != tileSet.end(); it++) {
        Metatile *tile = *it;
        if(lastZoomLevel != tile->z) {
            if(tilesRemoved > 0) {
                if(verbose) printf("Removed %u meta-duplicated tiles at zoom %u\n", tilesRemoved, lastZoomLevel);
            }
            tilesRemoved = 0;
            lastZoomLevel = tile->z;
        }
 
        for(unsigned int rx = tile->x; rx < tile->x + metaTiles; rx++) {
            for(unsigned int ry = tile->y; ry < tile->y + metaTiles; ry++) {
                if (rx==tile->x && ry==tile->y) {
                    continue; // don't remove self
                }
 
                Metatile delTile = Metatile(tile->z, rx, ry);
                unsigned int removed = tileSet.erase( &delTile );
                tilesRemoved += removed;
                tilesRemovedTotal += removed;
            }
        }
    }
 
    if(tilesRemoved > 0) {
        if(verbose) printf("Removed %u meta-duplicated tiles at zoom %u\n", tilesRemoved, lastZoomLevel);
    }
 
    if(tilesRemovedTotal > 0) {
        printf("Removed %u meta-duplicated tiles total, %lu remain\n", tilesRemovedTotal, tileSet.size());
    }
}
 
unsigned int lastTileZoomLevel = 0xffffffff;
 
/* returns the next tile from the tileQueue for a rendering thread */
Metatile *getNextTile() {
    Metatile *tile;
 
    pthread_mutex_lock(&tileQueueMutex);
        if (tileQueue.empty()) {
            tile = NULL;
        } else {
            tile = tileQueue.front();
            tileQueue.pop_front();
        }
    pthread_mutex_unlock(&tileQueueMutex);
 
    if(tile!=NULL && tile->z != lastTileZoomLevel) {
        lastTileZoomLevel = tile->z;
        printf("Working on zoom %u...\n", lastTileZoomLevel);
    }
 
    return tile;
}
 
void *renderThread(void *argument) {
    int tid = *((int *) argument);
 
    int metatile_count=0;
    int tile_count=0;
 
    // create a new meta tile as a mapnik map to render on
    mapnik::Map m(TILE_SIZE * metaTiles, TILE_SIZE * metaTiles);
 
    // load our map stylesheet
    mapnik::load_map(m, mapfile);
 
    Metatile* tile;
    while ((tile = getNextTile())!=NULL) {
        unsigned int x = tile->x;
        unsigned int y = tile->y;
        unsigned int z = tile->z;
 
        delete tile;
        tile = NULL;
 
        if(verbose) printf("%d: Rendering %d/%d/%d * %d\n", tid, z, x, y, metaTiles);
 
        Bbox bound_top_left = tileToMercBounds(z, x, y);
 
        int brx = x + metaTiles - 1;
        int bry = y + metaTiles - 1;
 
        // ensure the metatile doesn't extend beyond normal tile limits
        assert(brx < pow(2, z));
        assert(bry < pow(2, z));
 
        Bbox bound_bottom_right = tileToMercBounds(z, brx, bry);
        m.set_buffer_size(buffer);
        m.zoom_to_box(mapnik::box2d<double>(bound_top_left.left, bound_bottom_right.bottom, bound_bottom_right.right, bound_top_left.top));
 
        mapnik::image_32 buf(m.width(),m.height());
        mapnik::agg_renderer<mapnik::image_32> ren(m,buf);
        ren.apply();
 
        // make parent tile directories
        if (mkdir("tiles", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) && (errno != EEXIST)) {
            fprintf(stderr, "Can't make tile directory: tiles (%s)\n", strerror(errno));
            return NULL;
        }
 
        char tilesDir[FILENAME_MAX];
 
        // make z directory
        snprintf(tilesDir, FILENAME_MAX, "tiles/%d", z);
        if (mkdir(tilesDir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) && (errno != EEXIST)) {
            fprintf(stderr, "Can't make tile directory: tiles (%s)\n", strerror(errno));
            return NULL;
        }
 
        // make x directories (enough for all the split tiles too)
        for (unsigned int i = 0; i < metaTiles; i++) {
            snprintf(tilesDir, FILENAME_MAX, "tiles/%d/%d", z, x + i);
            if (mkdir(tilesDir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) && (errno != EEXIST)) {
                fprintf(stderr, "Can't make tile directory: tiles (%s)\n", strerror(errno));
                return NULL;
            }
        }
 
        // save the tile image in memory encoded using PNG
        std::string png_metaTile = mapnik::save_to_string<mapnik::image_data_32>(buf.data(), png256?"png256":"png");
 
        // split the metatile up into 256*256 size tiles
 
        // grab the in memory PNG encoded image as a Magick::Image
        Magick::Image IM_metaTile;
        IM_metaTile.magick("PNG");
        Magick::Blob IM_blob(png_metaTile.data(), png_metaTile.size());
        IM_metaTile.read(IM_blob);
 
        metatile_count++;
        for (unsigned int mi = 0; mi < metaTiles; mi++) {
            for(unsigned int mj = 0; mj < metaTiles; mj++) {
                char tileName[FILENAME_MAX];
                snprintf(tileName, FILENAME_MAX, "%s/%d/%d/%d.png", outputDirName.c_str(), z, x + mi, y + mj);
                Magick::Image IM_tile = IM_metaTile;
                IM_tile.crop( Magick::Geometry(256, 256, 256*mi, 256*mj) );
                IM_tile.write(tileName);
                tile_count++;
            }
        }
    }
 
    printf("Thread %d ending, processed %d metatiles (%d tiles)\n", tid, metatile_count, tile_count);
 
    return NULL;
}
 
/* For a given OSM z/x/y tile, find the bounds of this tile in the Spherical
   Mercator projection which OSM uses. */
Bbox tileToMercBounds(int z, int x, int y) {
    // special case for zoom 0 which spans the while projection space
    if (z == 0) {
        Bbox bbox = {-MAX_X, -MAX_Y, MAX_X, MAX_Y};
        return bbox;
    }
 
    /* the reference number for the tile directly after the centre axis of the
       projected space */
    int tilePastCentreAxis = pow(2,z) / 2;
 
    /* if tile is beyond expected extents cap it
         x = 0 .. (2^z - 1)
         y = 0 .. (2^z - 1)   */
    assert((x >= 0) && (x < pow(2, z)));
    assert((y >= 0) && (y < pow(2, z)));
 
    /* flip all tiles left/above of the centre axis across to the right/bottom
       side for now in other words move everything to the bottom right quadrant */
    bool mirrored_x = false;
    bool mirrored_y = false;
 
    if (x < tilePastCentreAxis) {
        mirrored_x = true;
        x = (pow(2, z) - 1) - x;
    }
 
    if (y < tilePastCentreAxis) {
        mirrored_y = true;
        y = (pow(2, z) - 1) - y;
    }
 
    // 2^z / 2 <= x < 2^z
    assert ((tilePastCentreAxis <= x) && (x < (pow(2, z))));
    assert ((tilePastCentreAxis <= y) && (y < (pow(2, z))));
 
    /* we need to cast the numerator to double in case both the numerator and
       and denominator would have been integers, but true answer should have
       been a floating point number.
 
       we also need to cast the tile number to double so that when we multiply
       it by MAX_X, we don't overfloat the int type */
    double bound_left = ((double)(((double)x - tilePastCentreAxis) * (MAX_X))) / tilePastCentreAxis;
    double bound_right = ((double)(((double)(x + 1) - tilePastCentreAxis)  * (MAX_X))) / tilePastCentreAxis;
 
    double bound_top = (((double)(((double)y - tilePastCentreAxis) * (MAX_Y))) / tilePastCentreAxis) * -1;
    double bound_bottom = (((double)(((double)(y + 1) - tilePastCentreAxis) * (MAX_Y))) / tilePastCentreAxis) * -1;
 
    /* if we mirrored the tiles earlier, fix them up now that we have the merc
       bounds */
    if (mirrored_x) {
        double br = bound_right;
        double bl = bound_left;
 
        bound_left = -1 * br;
        bound_right = -1 * bl;
    }
 
    if (mirrored_y) {
        double bb = bound_bottom;
        double bt = bound_top;
 
        bound_bottom = -1 * bt;
        bound_top = -1 * bb;
    }
 
    // return the bounds of this tile in the Spherical Mercator projection
    struct Bbox bbox = {bound_left, bound_bottom, bound_right, bound_top};
    return bbox;
}
