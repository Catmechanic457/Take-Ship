#ifndef NOISE_H
#define NOISE_H

#include "PerlinNoise.hpp"

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
        constexpr Settings(double frequency_, unsigned octaves_, double threshold_) : 
        frequency(frequency_), octaves(octaves_), threshold(threshold_) {}
    };
}

#endif