#ifndef NOISE_H
#define NOISE_H

#include "PerlinNoise.hpp"
#include "Validate.hpp"

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
        /**
         * \brief Load noise parameters into settings
         * \param values_ JSON values to load
        */
        void load_json(Json::Value values_) {

            std::string format = values_["format"].asString();

            Json::Value data_ = values_["data"];

            // only "simple" format implemented
            // "format" key in place only for future-proofing

            if (format == "simple") {
                frequency = data_["frequency"].asDouble();
                octaves = data_["octaves"].asUInt();
                threshold = data_["threshold"].asDouble();

                // validate
                if (frequency <= 0) {
                    send_error_message("frequency", error_type::EXPECTED_POSITIVE, ERROR);
                    std::cout << std::endl;
                    throw;
                }
                if (threshold < 0 or threshold > 1) {
                    send_error_message("threshold", error_type::NOT_NORMALIZED, WARN);
                    std::cout << std::endl;
                }
                return;
            }

            send_error_message("generic_noise_parameters", error_type::FORMAT_ERROR, ERROR);
            std::cout << std::endl;
            throw;
        }
    };
}

#endif