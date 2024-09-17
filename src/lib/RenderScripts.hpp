#ifndef RENDER_H
#define RENDER_H

#include <SFML/Graphics.hpp>

#include "RelativeSpace.hpp"

#include <vector>
#include <thread>

namespace render {
    /**
     * \return The dot product of the two vectors
    */
    template<typename T>
    T dot_product(std::vector<T> v1_, std::vector<T> v2_) {
        T total;
        for (unsigned int i = 0; i < v1_.size(); i++) {
            total += v1_[i] * v2_[i];
        }
        return total;
    }
    template<typename T>
    /**
     * \brief Converts the given vector into a unit vector
     * \param v Vector to convert
     * \return The normalized vector
    */
    std::vector<T> normalize(std::vector<T> v_) {
        // find the magnitude of the vector
        T mag_sq = 0;
        for (auto i : v_) {
            mag_sq += i*i;
        }
        auto mag = sqrt(mag_sq);

        // divide each term by the magnitude
        for (auto& i : v_) {
            i = i / mag;
        }
        return v_;
    }
    /**
     * \brief Blend two colours according to a weight
     * \param a_ First color
     * \param b_ Second color
     * \param w_ Weight towards first color
     * \return The blended color
    */
    sf::Color blend_colors(sf::Color a_, sf::Color b_, double w_) {
        // define a lambda function to mix two values according to a weight
        auto mix = [w_](double a, double b) {return a*w_ + b*(1-w_);};
        return sf::Color(
            // mix each channel
            mix(a_.r, b_.r),
            mix(a_.g, b_.g),
            mix(a_.b, b_.b),
            a_.a
            );
    }
    /**
     * \return The darkened colour
    */
    sf::Color darken_color(sf::Color color_, double w_) {
        return blend_colors(color_, sf::Color::Black, 1.0 - w_);
    }


    struct RenderContext {
        std::vector<double> light_direction;
        double ambient_brightness;
    };

    struct PixelInfo {
        sf::Color color;
        bool in_shadow;
        bool is_fertile_land = false;
    };

    struct RenderableLevel {
        virtual PixelInfo pixel(unsigned x_, unsigned y_, unsigned layer_index_) = 0;
        virtual rs::Vector2<unsigned> render_domain() = 0;
        virtual unsigned layer_count() = 0;
        virtual RenderContext render_context() = 0;
    };

    class Scene : public sf::RenderWindow {
        
        // the image buffer stores the layer textures during creation
        // outside of scene creation, this pointer will be null
        sf::Image* layer_image_buffer = nullptr;

        // stores each layer's texture
        std::vector<sf::Texture> layer_textures;

        // stores extra infomation about pixels in the static layer
        std::vector<std::vector<std::vector<PixelInfo>>> pixel_info;

        rs::Vector2<double> camera; // origin of local position
    
        /**
         * \brief Render a single layer
         * \param index_ Layer index to render
        */
        void build_single_layer(unsigned index_) { 
            auto domain = level.render_domain();
            layer_image_buffer[index_].create(domain.x, domain.y);
            // for each row (y)
            for (unsigned y = 0; y < domain.y; y++) {
                // for each column (x)
                for (unsigned x = 0; x < domain.x; x++) {
                    layer_image_buffer[index_].setPixel(x, y, render_pixel(level.pixel(x, y, index_)));
                }
            }
            layer_textures[index_].loadFromImage(layer_image_buffer[index_]);
            //layer_textures[index_].setSmooth(true);
        }
        /**
         * \brief Update the `static_layer` image and apply it to the background texture
        */
        void build_static_layer() {
            // reset pixel info
            auto domain = level.render_domain();
            unsigned layer_count = level.layer_count();
            pixel_info = std::vector<std::vector<std::vector<PixelInfo>>>(
                domain.y,
                std::vector<std::vector<PixelInfo>>(
                    domain.x,
                    std::vector<PixelInfo>(layer_count)
                )
            );

            // allocate memory for the image buffer
            layer_image_buffer = new sf::Image[layer_count];

            // create textures
            layer_textures = std::vector<sf::Texture>(layer_count);

            // start a thread for each layer in the level
            std::vector<std::thread> render_threads;
            for (unsigned i = 0; i < layer_count; i++) {
                std::thread t(&build_single_layer, this, i);
                render_threads.push_back(std::move(t));
            }
            // wait for all threads to end
            for (auto& t : render_threads) {t.join();}

            // deallocate memory
            delete[] layer_image_buffer;
            layer_image_buffer = nullptr;
        }
        /**
         * \return The colour of the rendered pixel
        */
        sf::Color render_pixel(PixelInfo p_) {

            auto base = p_.color;

            // darken pixel if in shadow
            double shadow_darkness = 0.25;
            if (p_.in_shadow) {base = darken_color(base, shadow_darkness);}

            // round colour channels for stylistic effect

            int r = 16;

            base.r = (base.r / r) * r;
            base.g = (base.g / r) * r;
            base.b = (base.b / r) * r;
            return base;
        }

        public:
        RenderableLevel& level;

        Scene(RenderableLevel& level_) : sf::RenderWindow(), level(level_) {}

        // stores the address of all drawables
        std::vector<sf::Drawable*> drawables;
        /**
         * \brief Clear the drawables array
        */
        void clear_drawables() {
            for (auto* drawable : drawables) delete drawable;
            drawables.clear();
        }
        /**
         * \brief Removes the given drawable
         * \param drawable_ Drawable to remove
        */
        void remove_drawable(sf::Drawable* t_drawable_) {
            for (unsigned i = 0; i < drawables.size(); i++) {
                // check if addresses match
                if (t_drawable_ == drawables[i]) {
                    drawables.erase(drawables.begin()+i);
                    return;
                }
            }
            std::cout << "Deletion of " << t_drawable_ << " failed. Entity not found.";
            throw;
        }
        render::RenderContext get_render_context() {return level.render_context();}
        /**
         * \return The camera position
        */
        rs::Vector2<double> get_camera() {return camera;}
        /**
         * \brief Set the camera position
         * \param pos_ New position
        */
        void set_camera(rs::Vector2<double> pos_) {
            camera = pos_;
        }
        /**
         * \brief Set the camera position
         * \param x_
         * \param y_
        */
        void set_camera(int x_, int y_) {set_camera(rs::Vector2(x_, y_));}
        /**
         * \brief Build a new level using the current seed
        */
        void generate() {
            build_static_layer();
        }
        /**
         * \brief Render the drawables
        */
        void render_drawables() {
            for (auto& drawable : drawables) {
                draw(*drawable);
            }
        }
        /**
         * \brief Render the static layer
        */
        void render_static_layer() {
            // draw each layer texture to the screen
            for (auto& t : layer_textures) {
                sf::Sprite s(t);
                draw(s);
            }
        }
        /**
         * \brief Render the whole scene
        */
        void render() {
            render_static_layer();
            render_drawables();
        }
        /**
         * \return Pixel info for a specific point
        */
        std::vector<render::PixelInfo> get_pixel_info(unsigned x_, unsigned y_) {
            return pixel_info[y_][x_];
        }
    };
}

#endif