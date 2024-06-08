#ifndef GAMEOBJ_H
#define GAMEOBJ_H

#include <SFML/Graphics.hpp>
#include "Noise.hpp"
#include "RenderScripts.hpp"

const noise::Settings StageDefault(8.0,3,0.6);
const noise::Settings ForestDefault(0.75,1,0.5);

namespace game {
    class Stage : public noise::Noise {
        private:
        rs::Vector2<unsigned> _win;
        rs::Vector2<unsigned> _sp;

        noise::Settings _settings;

        unsigned _win_max = _win.x > _win.y ? _win.x : _win.y;

        //protected:
        public:
        /**
         * \param pos_ Point
         * \return The amplitude of the noise at the given point
        */
        double value(rs::Vector2<double> pos_) const {
            return value(pos_, _settings.frequency, _settings.octaves);
        }
        /**
         * \param pos_ Point
         * \param f_ Frequency to use
         * \param o_ Number of octaves to use
         * \return The amplitude of the noise at the given point
        */
        double value(rs::Vector2<double> pos_, double f_, unsigned o_) const {
            return noise::Noise::value(pos_.x, pos_.y, f_ / _win_max, o_);
        }

        public:

        Stage(rs::Vector2<unsigned> size_, noise::Settings settings_) : noise::Noise(), _win(size_), _sp(size_.x/2, size_.y/2), _settings(settings_) {}
        Stage(unsigned sizex_, unsigned sizey_, noise::Settings settings_) : Stage(rs::Vector2<unsigned>(sizex_, sizey_), settings_) {}

        /**
         * \brief Set the number of octaves of noise
         * \param o_ New octave value
        */
        void set_octaves(unsigned o_) {_settings.octaves = o_;}
        /**
         * \brief Set the frequency of the noise
         * \param f_ New frequency value
        */
        void set_frequency(double f_) {_settings.frequency = f_;}
        /**
         * \brief Set the threshold ("height" value) used for the map
         * \param t_ New threshold value
        */
        void set_threshold(double t_) {_settings.threshold = t_;}
        /**
         * \brief Set the generation settings
         * \param s_ New settings object
        */
        void set_settings(noise::Settings s_) {_settings = s_;}
        /**
         * \returns The generation settings as a non-const reference
        */
        noise::Settings& get_settings() {return _settings;}
        /**
         * \brief Check a point lies inside a collision area
         * \param pos_ The point to check
         * \return `true` if the point is inside a collision area
        */
        bool collision(rs::Vector2<double> pos_) const {return value(pos_) > _settings.threshold;}
        /**
         * \brief Check a point lies inside the stage confines
         * \param pos_ The point to check
         * \return `true` if the point is inside stage confines
        */
        bool in_bounds(rs::Vector2<double> pos_) const {return (pos_.x >= 0 and pos_.y >= 0 and pos_.x < _win.x and pos_.y < _win.y);}
        /**
         * \brief Check a point is a valid position.
         * Where it is in bounds and not colliding
         * \param pos_ The point to check
         * \return `true` if the point is valid
        */
        bool valid(rs::Vector2<double> pos_) const {return in_bounds(pos_) and !collision(pos_);}
        bool path_exists(rs::Vector2<unsigned> pos1_, rs::Vector2<unsigned> pos2_, unsigned long max_iterations_ = -1) const {

            // if either the start or goal point is invalid, return false
            if (!valid(pos1_) or !valid(pos2_)) {return false;}

            const unsigned cell_size = 10;

            rs::Vector2<unsigned> start_cell(pos1_.x/cell_size, pos1_.y/cell_size);
            rs::Vector2<unsigned> goal_cell(pos2_.x/cell_size, pos2_.y/cell_size);

            // stores the directions to move in [N,E,S,W]
            int neighbour_pattern[4][2] = {
                {1,0},
                {0,1},
                {-1,0},
                {0,-1}
            };

            struct Cell {
                rs::Vector2<unsigned> pos;
                double cost = __DBL_MAX__;
                Cell(rs::Vector2<unsigned> pos_) : pos(pos_) {}
                bool operator== (rs::Vector2<unsigned>pos_) {
                    return pos == pos_;
                }
            };
            
            // stores the cells that are available to move to
            std::vector<Cell> open;
            // stores cells that have already been evaluated
            std::vector<rs::Vector2<unsigned>> closed;

            // initially, only the start cell can be moved to
            open.push_back(Cell(start_cell));

            unsigned long iteration = 0;

            while (open.size() > 0 and iteration < max_iterations_) {

                iteration++;

                // get the cell with the lowest cost
                auto current_it = std::max_element(
                open.begin(), open.end(),
                [] (const Cell &a, const Cell &b)
                {
                    return a.cost > b.cost;
                }
                );

                auto current = *current_it;

                // if the current cell is the goal, a path has been found
                if (current == goal_cell) {return true;}

                // since we have moved to the current cell, remove it from the open list
                open.erase(current_it);
                // mark the current cell as evaluated
                closed.push_back(current.pos);

                for (auto neighbour : neighbour_pattern) {
                    auto neighbour_pos = rs::Vector2<unsigned>(current.pos.x + neighbour[0], current.pos.y + neighbour[1]);
                    if (!valid(rs::Vector2(neighbour_pos.x*cell_size, neighbour_pos.y*cell_size))) {continue;}
                    if (std::find(open.begin(), open.end(), neighbour_pos) == open.end()) {

                        if (std::find(closed.begin(), closed.end(), neighbour_pos) == closed.end()) {

                            // add each neighboring cell to the open list along with their cost
                            Cell new_cell(neighbour_pos);
                            auto dy = neighbour_pos.y - goal_cell.y;
                            auto dx = neighbour_pos.x - goal_cell.x;

                            new_cell.cost = dy + dx;

                            open.push_back(new_cell);
                        }
                    }
                }
            }
            // if there are no move valid moves or the iteration limit has been reached
            // mark the stage as impossible
            return false;
        }
        /**
         * \return The stage spawnpoint (located at the centre)
        */
        const rs::Vector2<unsigned>& spawn_point() const {return _sp;}
        /**
         * \return The stage size
        */
        const rs::Vector2<unsigned>& stage_size() const {return _win;}
    };

    class Island : public Stage {

        // define grass layer
        const double grass_layer_offset = 0.05;

        // define minimum steepness to draw a cliff
        const double cliff_threshold = 0.20;

        // height scale for visual details
        const double visual_height_multiplier = 100.0;

        const sf::Color grass = sf::Color(140, 200, 65);
        const sf::Color dirt = sf::Color(180, 130, 50);
        const sf::Color cliff = sf::Color(185, 180, 165);
        const sf::Color sand = sf::Color(255, 200, 120);
        const sf::Color water_shallow = sf::Color(55, 170, 165);
        const sf::Color water_deep = sf::Color(50, 70, 155);

        /**
         * \param x_
         * \param y_
         * \return The scaled height
        */
        double scaled_height(double x_, double y_) {
            return value(rs::Vector2(x_, y_)) * visual_height_multiplier;
        }
        /**
         * \brief Get the normal vector of the terrain at the given point
         * \param x_
         * \param y_
         * \return The normal as a normalized vector
        */
        std::vector<double> calculate_normal(double x_, double y_) {
            double h = 1.0; // the step amount (should be small)

            double base = scaled_height(x_, y_);

            double dx = scaled_height(x_ + h, y_) - base; // change in x
            double dy = scaled_height(x_, y_ + h) - base; // change in y

            std::vector<double> normal = {-dx, -dy, h}; // normal vector

            return render::normalize(normal); // return as unit vector
        }
        /**
         * \return `true` if the point is in shadow
        */
        bool in_shadow(double x_, double y_, double h_) {
            //return false;

            // trace starts at pixel
            double trace_x = x_;
            double trace_y = y_;
            double trace_h = h_;
            
            while (true) {

                double height = value(rs::Vector2(trace_x, trace_y));

                // if trace is below land, it's in shadow
                if (height > trace_h) {return true;}

                // if the trace reaches the max height, return
                if (trace_h > 1.0) {return false;}

                
                const double k = 5.0;

                // step based on distance to terrain
                double d = (trace_h - height); // distance to terrain below
                double step = 0.2 + (d * d * k * visual_height_multiplier);

                double dx = -light_direction[0] * step;
                double dy = -light_direction[1] * step;
                double dh = -(light_direction[2] / visual_height_multiplier) * step;

                // step towards light
                trace_x += dx; 
                trace_y += dy;
                trace_h += dh;
            }
        }

        public:

        // direction of incoming light (from the sun)
        const std::vector<double> light_direction = render::normalize(std::vector<double>{2.0,2.0, -1});

        // ambient brightness
        const double am_brightness = 1.0;

        Island(rs::Vector2<unsigned int> size_, noise::Settings terrain_settings_ = StageDefault) : Stage(size_, terrain_settings_) {}
        Island(unsigned int x_, unsigned int y_, noise::Settings terrain_settings_ = StageDefault) : Island(rs::Vector2<unsigned int>(x_, y_), terrain_settings_) {}

        render::PixelInfo land_pixel(unsigned x_, unsigned y_) {

            auto height = value(rs::Vector2(x_, y_));
            auto normal = calculate_normal(x_, y_);

            // A 3D vector pointing strait up
            std::vector<double> vertical = {0.0, 0.0, 1.0};


            // base starts as grass or sand
            sf::Color base_land = (height - get_settings().threshold) > grass_layer_offset ? grass : sand;

            // find the steepness
            double land_steepness = 1 - render::dot_product(normal, vertical);

            // draw cliff if the land is steep
            if (cliff_threshold < land_steepness) {
                // map the steepness to be from `0.0` to `1.0`
                double cliff_steepness = (land_steepness - cliff_threshold) / (1 - cliff_threshold);

                // base becomes mix of gray and brown
                base_land = render::blend_colors(cliff, dirt, cbrt(cliff_steepness));
            }

            // darken areas that receive light at an angle
            // dot is negated as normal opposes light direction
            auto dot = -render::dot_product(normal, light_direction);

            dot = dot > 0.0 ? dot : 0; // Assume dot is 0 if dot is negative

            double brightness = (dot + am_brightness)/(1.0 + am_brightness);

            base_land = render::darken_color(base_land, 1 - brightness);

            render::PixelInfo pixel;
            pixel.color = base_land;
            pixel.in_shadow = in_shadow(x_, y_, height);

            // is fertile land on grass
            pixel.is_fertile_land = height - get_settings().threshold > grass_layer_offset and cliff_threshold > land_steepness;

            return pixel;
        }

        render::PixelInfo water_pixel(unsigned x_, unsigned y_) {

            render::PixelInfo pixel;
            
            // find the terrain height
            double height = value(rs::Vector2(x_, y_));

            // water level is threshold
            double water_level = get_settings().threshold;

            // if the terrain is above the water level
            // return a transparent color (no water)
            if (height > water_level) {
                pixel.color = sf::Color::Transparent;
                pixel.in_shadow = false;
                return pixel;
            }

            // water depth is portion of height below the water level
            double water_depth = (water_level - height) / water_level;

            // blend the deep and light colors based on the depth
            sf::Color water_colour = render::blend_colors(water_deep, water_shallow, cbrt(water_depth));

            // make shallower areas more transparent
            const double offset = 2.0;
            water_colour.a = static_cast<uint8_t>(((cbrt(water_depth) + offset) / (1.0 + offset)) * 255);
            
            pixel.color = water_colour;
            pixel.in_shadow = in_shadow(x_, y_, water_level);

            return pixel;
        }
    };

    class RenderableIsland : public Island, public render::RenderableLevel {
        public:
        RenderableIsland(rs::Vector2<unsigned int> size_, noise::Settings terrain_settings_ = StageDefault) : Island(size_, terrain_settings_), RenderableLevel() {}
        RenderableIsland(unsigned int x_, unsigned int y_, noise::Settings terrain_settings_ = StageDefault) : Island(rs::Vector2<unsigned int>(x_, y_), terrain_settings_) {}
        render::PixelInfo pixel(unsigned x_, unsigned y_, unsigned layer_index_) override {
            switch (layer_index_) {
            case 0: return land_pixel(x_, y_);
            case 1: return water_pixel(x_, y_);
            default: throw;
            }
        }
        rs::Vector2<unsigned> render_domain() override {return stage_size();}
        unsigned layer_count() override {return 2U;}
        render::RenderContext render_context() override {
            render::RenderContext context;
            context.ambient_brightness = am_brightness;
            context.light_direction = light_direction;
            return context;
        }
    };

    class BasicEntity : public sf::Drawable {
        // texture to use when drawing
        sf::Texture texture_source;
        // reference to the parent scene
        render::Scene& scene;
        // sources
        ObjectFile files;

        public:

        // global position of the entity
        rs::Position position;


        BasicEntity(render::Scene& parent_scene_) : scene(parent_scene_) {}
        BasicEntity(render::Scene& parent_scene_, ObjectFile files_) : scene(parent_scene_), files(files_) {}
        /**
         * \brief Render a new texture and apply it to the sprite
         * \param context_ Render context
        */
        void update_texture(render::RenderContext context_) {
            // create copy of base texture
            sf::Image texture = files.texture;
            // convert deg to rad for trig
            auto r = position.rotation * M_PI/ 180.0;
            // get components of light dir
            auto x1 = context_.light_direction[0];
            auto y1 = context_.light_direction[0];
            auto z1 = context_.light_direction[0];
            std::vector<double> light_direction = {
                (cos(r)*x1) - (sin(r)*y1),
                (sin(r)*x1) - (cos(r)*y1),
                z1
            };

            for (unsigned i = 0; i < files.texture.getSize().x; i++) {
                for (unsigned j = 0; j < files.texture.getSize().y; j++) {
                // for every pixel:

                // get the pixel of the normal map
                sf::Color n = files.normal_map.getPixel(i, j);

                // convert to a normalized vector
                std::vector<double> normal = {n.r/255.0, n.g/255.0, n.b/255.0};
                // get the dot product from incoming light
                auto dot = render::dot_product(normal, light_direction);
                    dot = dot > 0.0 ? dot : 0; // Assume dot is 0 if dot is negative

                    double brightness = (dot + context_.ambient_brightness)/(1.0 + context_.ambient_brightness);

                    // darken pixel
                    texture.setPixel(i, j, render::darken_color(texture.getPixel(i, j), 1.0 - brightness));

                }
            }

            // assign the new texture to the sprite
            texture_source.loadFromImage(texture);
        }
        void draw(sf::RenderTarget& target, sf::RenderStates states) const override {
            // sprite object to position and draw the texture
            sf::Sprite stamp(texture_source);
            // get the global position and find its relative position
            auto global_pos = position.position;
            auto local_pos = scene.relative_from_global(global_pos.x, global_pos.y);
            // draw in the correct position
            stamp.setOrigin(texture_source.getSize().x/2.0,texture_source.getSize().y/2.0);
            stamp.setPosition(local_pos.x + scene.getSize().x/2, local_pos.y + scene.getSize().y/2);
            stamp.setRotation(position.rotation);
            target.draw(stamp);
        }
        /**
         * \brief Set object files
         * \param objf_ Reference to object files
        */
        void set_object_file(ObjectFile& objf_) {
            files = objf_;
        }
    };

    class ForestGenerator : public noise::Noise {
        // reference to parent scene
        render::Scene& scene;
        // settings to use
        noise::Settings settings;

        public:
        
        ForestGenerator(render::Scene& scene_, noise::Settings settings_ = ForestDefault) :
        scene(scene_), settings(settings_) {}
        /**
         * \brief Generate a forest and append the entities to the scene
         * \param max_count_ The maximum possible trees to generate. The fraction of trees that
         * will spawn is approximately half the threshold value multiplied by the land coverage.
        */
        void generate(unsigned max_count_ = 500) {
            // seed the noise
            set_seed(rand());
            // get the size of the scene
            auto domain = scene.getSize();
            for (unsigned i = 0; i < max_count_; i++) {
                // chose a random coordinate
                // TODO - Currently deterministic, needs to be reseeded to some random input
                unsigned x = rand() % domain.x;
                unsigned y = rand() % domain.y;

                // check if any layer in the pixel is fertile
                for (auto p : scene.get_pixel_info(x, y)) {
                    if (not p.is_fertile_land) {continue;}
                    // is fertile
                    // get the noise level
                    double n = value(x, y, settings.frequency, settings.octaves);
                    // rescale based on threshold
                    n = (n - settings.threshold) / (1.0 - settings.threshold);
                    if (n <= 0) break;
                    // get a random 3 s.f float
                    double r = rand() % 1000;

                    // place a tree if the noise is greater than the random number
                    if (n * 1000 > r) {
                        // allocate new memory
                        // TODO memory needs to be freed when the object is destroyed
                        BasicEntity* tree = new BasicEntity(scene, assets::object_file::tree);
                        tree->position.position = rs::Vector2<double>(x,y);
                        tree->position.rotation = rand() % 360;
                        tree->update_texture(scene.get_render_context());
                        scene.drawables.push_back(tree);
                        break;
                    }
                } 
            }
        }
    };
}

#endif