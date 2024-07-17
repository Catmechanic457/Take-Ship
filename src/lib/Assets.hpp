#ifndef ASSETS_H
#define ASSETS_H

#include <SFML/Graphics.hpp>
#include <iostream>
#include <string.h>

#define ASSETS_PATH "../assets/"
#define TEXTURE_PATH "textures/"
#define FONT_PATH "fonts/"


struct ObjectFile {
    sf::Image texture, height_map, normal_map;
};

/**
 * \brief Load object files under a namespace
 * \param obj_ target object
 * \param namespace_ namespace or path to use
 * \return `true` if loading was successful
*/
bool load_object_file(ObjectFile& obj_, std::string namespace_) {
    // load all files
    // store the success into the array
    bool load_confirm[] = {
        obj_.texture.loadFromFile(ASSETS_PATH TEXTURE_PATH + namespace_ + "_texture.png"),
        obj_.height_map.loadFromFile(ASSETS_PATH TEXTURE_PATH + namespace_ + "_hmap.png"),
        obj_.normal_map.loadFromFile(ASSETS_PATH TEXTURE_PATH + namespace_ + "_normap.png")
    };
    // check all reads were a success
    for (bool b : load_confirm) {
        if (!b) {return false;}
    }
    return true;
}

namespace assets {
    namespace object_file {
        ObjectFile debug;
        ObjectFile tree;
        ObjectFile player_ship;
    }
}

/**
 * \brief Load all game assets
 * \return `true` if loading was successful
*/
bool load_assets() {
    using namespace assets;
    using namespace object_file;
    // load all assets
    // store the success into the array
    bool load_confirm[] = {
        load_object_file(debug, "debug"),
        load_object_file(tree, "level/tree"),
        load_object_file(tree, "mobiles/player_ship")
    };
    // check loading was a success
    for (bool b : load_confirm) {
        if (!b) {return false;}
    }
    return true;
}

#endif