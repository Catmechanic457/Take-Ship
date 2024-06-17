#ifndef NOISE_H
#define NOISE_H

#include "PerlinNoise.hpp"

#include <string.h>
#include <fstream>
#include <json/json.h>

#include <iostream>


namespace noise {
    class Noise : private siv::PerlinNoise {
        unsigned _seed;
        public:
        // constructor - takes an optional seed
        Noise(unsigned seed_) : siv::PerlinNoise(seed_), _seed(seed_) {}
        Noise() : Noise(0) {}
        /**
         * \brief Set the seed to be used
         * \param seed_ New seed
        */
        void set_seed(unsigned seed_) {_seed = seed_; reseed(seed_);}
        /**
         * \returns The current seed
        */
        unsigned get_seed() const {return _seed;}
        /**
         * \param x_ X position
         * \param y_ Y position
         * \param f_ Frequency to use
         * \param o_ Number of octaves to use
         * \return The amplitude of the noise at the given point
        */
        double value(double x_, double y_, double frequency_, unsigned octaves_) const {
            return octave2D_01(x_ * frequency_, y_ * frequency_, octaves_);
        }
    };
    struct Settings {
        double frequency;
        unsigned octaves;
        double threshold;
        constexpr Settings() : Settings(0.0,0,0.0) {}
        constexpr Settings(double frequency_, unsigned octaves_, double threshold_) : 
        frequency(frequency_), octaves(octaves_), threshold(threshold_) {}
        void load_from_path(std::string path_) {
            Json::Value file;
            std::ifstream values_file(path_, std::ifstream::binary);
            if (values_file.good()) {
                values_file >> file;
            }
            load_json(file["noise_parameters"]);
        }
        void load_json(Json::Value file_) {

            std::string format = file_["format"].asString();

            Json::Value data_ = file_["data"];

            // only "simple" format implemented
            // "format" key in place only for future-proofing

            if (format == "simple") {
                frequency = data_["frequency"].asDouble();
                octaves = data_["octaves"].asUInt();
                threshold = data_["threshold"].asDouble();
                return;
            }
            if (format == "path") {
                load_from_path(data_["path"].asString());
                return;
            }

            std::cout << "Invalid Format!" << std::endl;
        }
    };
}

#endif