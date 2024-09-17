#ifndef GAMEOBJ_H
#define GAMEOBJ_H

#include <json/json.h>
#include <random>

#include <SFML/Graphics.hpp>
#include "Noise.hpp"
#include "RenderScripts.hpp"

#define EULERS_NUM 2.718281828459

const noise::Settings StageDefault(8.0,3,0.6);
const noise::Settings ForestDefault(0.75,1,0.5);

namespace game {
    class Stage : public noise::Noise {
        private:
        rs::Vector2<unsigned> _win;
        
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

        Stage(rs::Vector2<unsigned> size_, noise::Settings settings_) : noise::Noise(), _win(size_), _settings(settings_) {}
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
         * \return The stage size
        */
        const rs::Vector2<unsigned>& get_stage_size() const {return _win;}
        /**
         * \brief Set the stage size
        */
        void set_stage_size(rs::Vector2<unsigned> size_) {_win = size_;}
        /**
         * \brief Set the stage size
        */
        void set_stage_size(unsigned x_, unsigned y_) {_win.x = x_; _win.y = y_;}
    };

    class Island : public Stage {

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

        // define grass layer
        double grass_layer_offset = 0.05;

        // define minimum steepness to draw a cliff
        double cliff_threshold = 0.20;

        // height scale for visual details
        double visual_height_multiplier = 100.0;

        sf::Color grass = sf::Color(140, 200, 65);
        sf::Color dirt = sf::Color(180, 130, 50);
        sf::Color cliff = sf::Color(185, 180, 165);
        sf::Color sand = sf::Color(255, 200, 120);
        sf::Color water_shallow = sf::Color(55, 170, 165);
        sf::Color water_deep = sf::Color(50, 70, 155);

        // direction of incoming light (from the sun)
        std::vector<double> light_direction = render::normalize(std::vector<double>{2.0,2.0, -1});

        // ambient brightness
        double am_brightness = 1.0;

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
        rs::Vector2<unsigned> render_domain() override {return get_stage_size();}
        unsigned layer_count() override {return 2U;}
        render::RenderContext render_context() override {
            render::RenderContext context;
            context.ambient_brightness = am_brightness;
            context.light_direction = light_direction;
            return context;
        }
    };

    class BasicEntity : public sf::Drawable {
        protected:
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
            auto r = position.rotation * 3.14159 / 180.0;
            // get components of light dir
            auto x1 = context_.light_direction[0];
            auto y1 = context_.light_direction[1];
            auto z1 = context_.light_direction[2];
            std::vector<double> light_direction = {
                (cos(-r)*x1) - (sin(-r)*y1),
                (cos(-r)*y1) + (sin(-r)*x1),
                z1
            };

            for (unsigned i = 0; i < files.texture.getSize().x; i++) {
                for (unsigned j = 0; j < files.texture.getSize().y; j++) {
                // for every pixel:

                    // get the pixel of the normal map
                    sf::Color n = files.normal_map.getPixel(i, j);

                    // convert to a normalized vector
                    std::vector<double> normal = render::normalize<double>({n.r - 127.0, n.g - 127.0, -(n.b - 127.0)});
                    // get the dot product from incoming light
                    auto dot = render::dot_product(normal, light_direction);

                    dot = dot > 0.0 ? dot : 0; // Assume dot is 0 if dot is negative
                    // dot must be clamped. ideally, this should be avoided if possible
                    dot = dot < 1.0 ? dot : 1; // Assume dot is 1 if greater than 1 (due to rounding errors)

                    double brightness = (dot + context_.ambient_brightness)/(1.0 + context_.ambient_brightness);

                    // darken pixel
                    texture.setPixel(i, j, render::darken_color(texture.getPixel(i, j), 1.0 - brightness));

                }
            }

            // assign the new texture to the sprite
            texture_source.loadFromImage(texture); // When player killed, "Unknown Stopping event" fault
        }
        void draw(sf::RenderTarget& target, sf::RenderStates states) const override {
            // set sprite object to position and draw the texture
            sf::Sprite stamp(texture_source);
            auto global_pos = position.position;
            // draw in the correct position
            stamp.setOrigin(texture_source.getSize().x/2.0,texture_source.getSize().y/2.0);
            stamp.setPosition(global_pos.x, global_pos.y);
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

    class ComplexEntity : public BasicEntity {
        // region around centre of entity where it can collide
        const rs::Vector2<double> hit_box;

        const double max_health;
        const double armor;
        const double mass;
        // how much the entity slows down over time, based on its current velocity
        const double drag_coefficient;
        // determines how the entity behaves in collisions
        const double elasticity;

        public:

        double health;
        bool can_collide = true;
        rs::Vector2<double> velocity;
        // the resultant force currently applied to spire
        // excludes force from drag
        rs::Vector2<double> resultant_force;

        ComplexEntity(render::Scene& parent_scene_, rs::Vector2<double> hit_box_, double base_health_ = 100.0, double armor_ = 0.0, double mass_ = 100.0, double drag_coeff_ = 0.2, double elasticity_ = 1.0) :
        BasicEntity(parent_scene_), hit_box(hit_box_), max_health(base_health_), armor(armor_), mass(mass_), drag_coefficient(drag_coeff_), elasticity(elasticity_), health(max_health), velocity(0.0,0.0), resultant_force(0.0,0.0)
        {}

        ComplexEntity(render::Scene& parent_scene_, ObjectFile files_, rs::Vector2<double> hit_box_, double base_health_ = 100.0, double armor_ = 0.0, double mass_ = 100.0, double drag_coeff_ = 0.2, double elasticity_ = 1.0) :
        ComplexEntity(parent_scene_, hit_box_, base_health_, armor_, mass_, drag_coeff_, elasticity_)
        {
            set_object_file(files_);
        }
        /**
         * \return The entity's hit-box
        */
        rs::Vector2<double> get_hitbox() {return hit_box;}
        /**
         * \return The maximum health
        */
        double get_max_health() {return max_health;}
        /**
         * \return The mass of the entity
        */
        double get_mass() {return mass;}
        /**
         * \return The entity's elasticity in collisions
        */
        double get_elasticity() {return elasticity;}
        /**
         * \brief Apply damage to the entity
         * \param damage Damage to receive
         * \param piercing Ignore armour if `true`
        */
        void apply_damage(double damage, bool piercing = false) {
            double dmg = piercing ? damage : damage / armor;
            health -= dmg;
        }
        /**
         * \return `true` if the health is below 0
        */
        bool is_killed() {return health <= 0.0;}
        /**
         * \brief Rotate the entity based on how it's currently moving
        */
        void goto_natural_rotation() {
            // weight towards velocity
            double k = 10;

            // add the velocity and acceleration vectors
            double i = resultant_force.x / mass + k * velocity.x;
            double j = resultant_force.y / mass + k * velocity.y;

            // use trig to find angle
            double angle = atan2(j, i) / 3.14159 * 180.0;

            // ensure angle in range 0 to 350
            while (angle < 0) {
                angle += 360.0;
            }
            while (angle >= 360) {
                angle -= 360.0;
            }
            
            // set rotation
            position.rotation = angle;
        }
        /**
         * \brief Calculates the sprites new position
         * \param t Time in seconds that has elapsed
        */
        void move_over_time(double t) {
            
            // equation for final velocity
            auto v = [&] (double u, double a_0) {
                // let drag coeff = `k`
                auto& k = drag_coefficient;
                return ((a_0 - (u * k)) * pow(EULERS_NUM, -k * t) - a_0) / -k;
            };

            // equation for displacement
            auto s = [&] (double u, double a_0) {
                // let drag coeff = `k`
                auto& k = drag_coefficient;
                return ((a_0 - (u * k)) * (pow(EULERS_NUM, -k * t) - 1) + (a_0 * t * k)) / (k*k);
            };

            // `i` represents x direction
            // `j` represents y direction

            // find acceleration from resultant force
            double acc_i = resultant_force.x / mass;
            double acc_j = resultant_force.y / mass;

            // calculate final velocity
            double v_i = v(velocity.x, acc_i);
            double v_j = v(velocity.y, acc_j);

            // calculate displacement
            double dx = s(velocity.x, acc_i);
            double dy = s(velocity.y, acc_j);

            // set entity velocity to final velocity
            velocity.x = v_i;
            velocity.y = v_j;

            // add displacement to entity position
            position.position.x += dx;
            position.position.y += dy;
        }
        /**
         * \brief Check if the given hit-box intersects the entity's hit-box
         * \param o_ Position of hit-box to check
         * \param hit_box_ Size of hit-box
        */
        bool collides_with(rs::Vector2<double> o_, rs::Vector2<double> hit_box_ = rs::Vector2<double>(0,0)) {
            rs::Vector2<double> p1(position.position.x + (hit_box.x/2.0), position.position.y + (hit_box.y/2.0));
            rs::Vector2<double> p2(position.position.x - (hit_box.x/2.0), position.position.y - (hit_box.y/2.0));
            rs::Vector2<double> p3(p1.x, p2.y);
            rs::Vector2<double> p4(p2.x, p1.y);

            rs::Vector2<double> p5(o_.x + (hit_box_.x/2.0), o_.y + (hit_box_.y/2.0));
            rs::Vector2<double> p6(o_.x - (hit_box_.x/2.0), o_.y - (hit_box_.y/2.0));
            rs::Vector2<double> p7(p5.x, p6.y);
            rs::Vector2<double> p8(p6.x, p5.y);

            auto check_intersect = [] (rs::Vector2<double> p1_, rs::Vector2<double> p2_, rs::Vector2<double> cp_) {
                return p1_.x >= cp_.x and cp_.x >= p2_.x and p1_.y >= cp_.y and cp_.y >= p2_.y;
            };

            // check if any point in the check hit-box intersects this hit-box
            if (check_intersect(p1, p2, p5)) {return true;}
            if (check_intersect(p1, p2, p6)) {return true;}
            if (check_intersect(p1, p2, p7)) {return true;}
            if (check_intersect(p1, p2, p8)) {return true;}

            // check if any point in this hit-box intersects the check hit-box
            if (check_intersect(p5, p6, p1)) {return true;}
            if (check_intersect(p5, p6, p2)) {return true;}
            if (check_intersect(p5, p6, p3)) {return true;}
            if (check_intersect(p5, p6, p4)) {return true;}

            // if all checks fail, return false
            return false;
        }
        virtual void kill() {delete this;}
    };
    struct Weapon {
        // cooldown between weapon uses in seconds
        double base_cooldown;
        // remaining cooldown in seconds
        double current_cooldown = 0.0;
        // the velocity projectiles are launched
        double launch_velocity;
        // how long the projectile stays in the air
        double flight_time = __DBL_MAX__;
        // relative placement of weapon in relation to entities centre
        rs::Position origin;
    };

    class HostileEntity : public ComplexEntity {
        public:
        // list of all weapons belonging to the entity
        std::vector<Weapon> weapons;

        HostileEntity(render::Scene& parent_scene_, rs::Vector2<double> hit_box_, double base_health_ = 100.0, double armor_ = 0.0, double mass_ = 100.0, double drag_coeff_ = 0.2, double elasticity_ = 1.0) :
        ComplexEntity(parent_scene_, hit_box_, base_health_, armor_, mass_, drag_coeff_, elasticity_)
        {}
        HostileEntity(render::Scene& parent_scene_, ObjectFile files_, rs::Vector2<double> hit_box_, double base_health_ = 100.0, double armor_ = 0.0, double mass_ = 100.0, double drag_coeff_ = 0.2, double elasticity_ = 1.0) :
        ComplexEntity(parent_scene_, files_, hit_box_, base_health_, armor_, mass_, drag_coeff_, elasticity_)
        {}
        /**
         * \return `true` if weapon successfully fired
        */
        bool try_fire(unsigned index_) {
            // get the memory address of the weapon
            auto& weapon = weapons[index_];
            // check it's not on cooldown
            if (weapon.current_cooldown <= 0.0) {
                // reset cooldown
                weapon.current_cooldown = weapon.base_cooldown;
                return true;
            }
            return false;
        }
        rs::Position weapon_position(unsigned index_) {
            auto weapon = weapons[index_];
            auto r = (position.rotation + weapon.origin.rotation) * 3.14159 / 180.0;
            auto x = weapon.origin.position.x;
            auto y = weapon.origin.position.y;
            double i = (cos(r)*x) - (sin(r)*y);
            double j = (cos(r)*y) + (sin(r)*x);
            rs::Position pos;
            pos.position.x = position.position.x + i;
            pos.position.y = position.position.y + j;
            pos.rotation = position.rotation + weapon.origin.rotation;

            return pos;
        }
    };

    class Projectile : public ComplexEntity {
        public:
        Projectile(render::Scene& parent_scene_, rs::Vector2<double> hit_box_, double base_health_ = 1.0, double armor_ = 0.0, double mass_ = 100.0, double drag_coeff_ = 0.05, double elasticity_ = 0.4) :
        ComplexEntity(parent_scene_, hit_box_, base_health_, armor_, mass_, drag_coeff_, elasticity_)
        {}
        Projectile(render::Scene& parent_scene_, ObjectFile files_, rs::Vector2<double> hit_box_, double base_health_ = 1.0, double armor_ = 0.0, double mass_ = 100.0, double drag_coeff_ = 0.05, double elasticity_ = 0.4) :
        ComplexEntity(parent_scene_, files_, hit_box_, base_health_, armor_, mass_, drag_coeff_, elasticity_)
        {}
        double r_collision_grace = 0.25;
        double r_flight_time;
        void kill() override {r_flight_time = 0;}
    };

    struct HostileController {
        virtual void action(HostileEntity& entity) = 0;
    };
    struct ControlledHostile {
        HostileEntity* entity;
        HostileController* controller;
    };

    class ForestGenerator : public noise::Noise {
        // reference to parent scene
        render::Scene& scene;
        //game::Stage& stage; // ?? idk
        // settings to use
        noise::Settings settings;

        public:
        
        ForestGenerator(render::Scene& scene_, noise::Settings settings_ = ForestDefault) :
        scene(scene_), settings(settings_) {}
        /**
         * \brief Set the generation settings
         * \param s_ New settings object
        */
        void set_settings(noise::Settings s_) {settings = s_;}
        /**
         * \returns The generation settings as a non-const reference
        */
        noise::Settings& get_settings() {return settings;}
        /**
         * \brief Generate a forest and append the entities to the scene
         * \param max_count_ The maximum possible trees to generate. The fraction of trees that
         * will spawn is approximately half the threshold value multiplied by the land coverage.
        */
        void generate(unsigned max_count_ = 500) {
            // seed the noise
            set_seed(rand());
            // get the size of the scene
            auto domain = scene.level.render_domain();
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
    class CollisionHandler {
        // Determines how much damage to apply from a collision
        static constexpr double collision_damage_constant = 0.01;
        public:
        CollisionHandler() = delete; // No instances should be made
        /**
         * \brief Handle a *known* collision between two entities
         * \param e1_ Entity 1
         * \param e1_ Entity 2
        */
        static void handle_collision(ComplexEntity& e1_, ComplexEntity& e2_) {
            // momentum of e1
            double m1x = e1_.velocity.x * e1_.get_mass();
            double m1y = e1_.velocity.y * e1_.get_mass();

            // momentum of e2
            double m2x = e2_.velocity.x * e2_.get_mass();
            double m2y = e2_.velocity.y * e2_.get_mass();

            // when elasticity = 0, k = 0.5
            
            double k = (e1_.get_elasticity() + e2_.get_elasticity() + 2) / 4.0;

            // find new momentum for each entity
            double m3x = m2x * k + m1x * (1-k);
            double m3y = m2y * k + m1y * (1-k);

            double m4x = m1x * k + m2x * (1-k);
            double m4y = m1y * k + m2y * (1-k);

            // damage is proportional to change in momentum
            double damage1 = sqrt((m3x - m1x)*(m3x - m1x) + (m3y - m1y)*(m3y - m1y)) * collision_damage_constant;
            double damage2 = sqrt((m4x - m2x)*(m4x - m2x) + (m4y - m2y)*(m4y - m2y)) * collision_damage_constant;

            e1_.apply_damage(damage1);
            e2_.apply_damage(damage2);

            // find new velocity
            e1_.velocity.x = m3x / e1_.get_mass();
            e1_.velocity.y = m3y / e1_.get_mass();

            e2_.velocity.x = m4x / e2_.get_mass();
            e2_.velocity.y = m4y / e2_.get_mass();
        }
    };
    class GameHandler {

        render::Scene& scene;
        game::Stage& stage;
        HostileEntity* player_ship;
        Projectile cannonball;

        // stores a reference to every mobile on the scene
        std::vector<ComplexEntity*> mobiles;

        // stores a reference to every enemy on the scene
        std::vector<ComplexEntity*> enemies;

        // stores a reference to every projectile on the scene
        std::vector<Projectile*> projectiles;

        // stores the controller object for each controlled mobile
        std::vector<ControlledHostile*> controlled_hostiles;
        /**
         * \brief Add an entity to the scene
        */
        void add_entity(ComplexEntity* e_) {
            mobiles.push_back(e_);
            scene.drawables.push_back(e_);
        }
        /**
         * \brief Remove the given mobile from the mobile array. Does not delete the pointer
         * \param t_mobile_ Mobile to remove
        */
        void remove_mobile(ComplexEntity* t_mobile_) {
            for (unsigned i = 0; i < mobiles.size(); i++) {
                // check if addresses match
                if (t_mobile_ == mobiles[i]) {
                    mobiles.erase(mobiles.begin()+i);
                    return;
                }
            }
            // throw if not found
            throw;
        }
        /**
         * \brief Delete all mobiles
        */
        void delete_mobiles() {
            for (auto* mobile : mobiles) delete mobile; // free memory
            mobiles.clear(); // clear the array
        }
        /**
         * \brief Remove the given mobile from the mobile array and drawables
         * \param t_mobile_ Mobile to kill
        */
        void kill_mobile(ComplexEntity* t_mobile_) {
            std::cout << "Deleting entity `" << t_mobile_ << "`\n";

            t_mobile_->health = 0.0;
            scene.remove_drawable(t_mobile_);
            remove_mobile(t_mobile_);
            t_mobile_->kill(); // mobile may be left undeleted

            std::cout << "Mobile array size is now " << mobiles.size() << '\n';
        }
        /**
         * \brief Detect and handle any collisions between mobiles
        */
        void handle_collisions() {
            // for all mobiles
            // note - last mobile is skipped as all pairs will already be exhausted
            for (unsigned i = 0; i < mobiles.size() - 1; i++) {
                if (mobiles[i]->can_collide == false) {continue;}
                // check remaining pairs
                for (auto j = i + 1; j < mobiles.size(); j++) {
                    if (mobiles[j]->can_collide == false) {continue;}
                    // check for collision
                    if (mobiles[i]->collides_with(mobiles[j]->position.position, mobiles[j]->get_hitbox())) {
                        // handle collision
                        CollisionHandler::handle_collision(*mobiles[i], *mobiles[j]);
                    }
                }
            }
        }
        /**
         * \brief Move all mobiles to their next positions
         * \param t_ Time elapsed
        */
        void move_mobiles(double t_) {
            for (auto mobile : mobiles) {
                mobile->move_over_time(t_);
                mobile->goto_natural_rotation();
            }
        }
        /**
         * \brief Remove mobiles with a health that is < 0
        */
        void remove_killed_mobiles() {
            for (auto mobile : mobiles) {
                if (mobile->is_killed()) kill_mobile(mobile);
            }
        }
        /**
         * \brief Update the textures of all mobiles
        */
        void update_textures() {
            for (auto mobile : mobiles) {mobile->update_texture(scene.get_render_context());}
        }
        /**
         * \brief Initialises a new player ship
         * \return A pointer to the new player
        */
        HostileEntity* create_player() {

            // create a new player
            player_ship = new HostileEntity(scene, assets::object_file::player_ship, rs::Vector2(32.0,32.0), 50, 3, 75, 0.4, 0.8);

            // define left and right weapons
            Weapon player_cannon_left, player_cannon_right;

            // use left cannon to define properties
            player_cannon_left.base_cooldown = 1.0; // change to 1.0
            player_cannon_left.launch_velocity = 100.0;
            player_cannon_left.flight_time = 4.0;
            player_cannon_left.origin.position.x = 0.0;
            player_cannon_left.origin.position.y = 10.0; // change
            player_cannon_left.origin.rotation = 180.0;

            // modify left side properties for right side
            player_cannon_right = player_cannon_left;
            player_cannon_right.origin.position.y *= 1; // flip y-axis
            player_cannon_right.origin.rotation -= 180.0; // rotate 180

            // add weapons to the player ship
            player_ship->weapons.push_back(player_cannon_left);
            player_ship->weapons.push_back(player_cannon_right);

            // return the player object
            return player_ship;
        }
        /**
         * \brief Resets all entities for a new round
        */
        void reset_entities() {
            
            // clear all old entities from memory
            delete_mobiles();
            scene.clear_drawables();

            // delete the old player
            //delete player_ship; seg fault

            // create a new player
            create_player();
            add_entity(player_ship);

            // TODO - create the controller objects

            // create 3 tanks
            for (unsigned i = 0; i < 3; i++) {
                enemies.push_back(new HostileEntity(scene, assets::object_file::tank, rs::Vector2(32.0,32.0), 150, 4, 145, 0.4, 0.4));
            }

            // create 3 snipers
            for (unsigned i = 0; i < 3; i++) {
                enemies.push_back(new HostileEntity(scene, assets::object_file::sniper, rs::Vector2(32.0,32.0), 75, 1, 80, 0.4, 0.8));
            }

            // push all enemies onto the mobile array
            for (auto* e : enemies) {
                add_entity(e);
            }
        }
        /**
         * \return `true` if the given position is a valid spawn for a ship entity
        */
        bool valid_spawn(double x_, double y_) {
            int check_pattern[9][2] = {
                {1,1}, {1,0}, {1,-1},
                {0,1}, {0,0}, {0,-1},
                {-1,1}, {-1,0}, {-1,-1}
            };

            double scale = 32.0;

            // if any point in the check pattern is invalid, spawn is invalid
            for (auto p : check_pattern) {
                if (!stage.valid(rs::Vector2(x_ + (p[0] * scale), y_ + (p[1] * scale)))) {
                    return false;
                }
            }

            return true;
        }
        /**
         * \brief Spawn the player and enemy ships
        */
        void spawn_ships() {

            // midpoint of level
            auto mx = scene.level.render_domain().x/2;
            auto my = scene.level.render_domain().y/2;

            // spawn player at centre

            player_ship->position.position.x = mx;
            player_ship->position.position.y = my;
            player_ship->position.rotation = 0;

            // TODO - Use a better method for random coordinates
            
            // "minimum player separation" - the closest distance an enemy can spawn from the player
            double mps = 180.0;

            // "minimum enemy separation" - the closest distance an enemy can spawn from another enemy
            double mes = 32.0;

            for (unsigned i = 0; i < enemies.size(); i++) {
                while (true) {
                    
                    try_spawn_enemy:

                    int x = rand() % scene.level.render_domain().x;
                    int y = rand() % scene.level.render_domain().y;

                    auto dx = x - scene.level.render_domain().x/2;
                    auto dy = y - scene.level.render_domain().y/2;

                    // check mps
                    if ((dx*dx + dy*dy) < mps*mps) {goto try_spawn_enemy;}

                    // check mes
                    for (unsigned j = 0; j < i; j++) {
                        // for all already placed enemies
                        auto dx = x - enemies[j]->position.position.x;
                        auto dy = y - enemies[j]->position.position.y;

                        if ((dx*dx + dy*dy) < mes*mes) {goto try_spawn_enemy;}
                    }

                    // check valid
                    if (!valid_spawn(x,y)) {goto try_spawn_enemy;}

                    // ensure there exists a valid path between the enemy and player
                    // this ensures that the enemy is able to reach the player
                    // this is an expensive process. optimisation steps should start here
                    if (!stage.path_exists(rs::Vector2(mx,my), rs::Vector2(x,y))) {goto try_spawn_enemy;}

                    enemies[i]->position.position.x = x;
                    enemies[i]->position.position.y = y;
                    enemies[i]->position.rotation = 0; // face centre maybe?

                    break;
                }
            }
        }
        /**
         * \brief Give the projectile velocity and add it to the projectile array
        */
        void launch_projectile(HostileEntity* parent_, unsigned weapon_index_, Projectile* projectile_) {
            auto& weapon = parent_->weapons[weapon_index_];

            double x = parent_->weapon_position(weapon_index_).position.x;
            double y = parent_->weapon_position(weapon_index_).position.y;

            double r = (parent_->position.rotation + weapon.origin.rotation) * 3.14159 / 180.0;
            double v = weapon.launch_velocity;

            double i = parent_->velocity.x;
            double j = parent_->velocity.y;

            projectile_->position.position.x = x;
            projectile_->position.position.y = y;
            projectile_->velocity.x = i + (sin(-r) * v);
            projectile_->velocity.y = j + (cos(-r) * v);

            projectile_->r_flight_time = weapon.flight_time;
            projectile_->r_collision_grace = 0.25;
            projectile_->can_collide = false;
            projectiles.push_back(projectile_);
        }
        /**
         * \brief Update the flight time and collision grace of projectiles
         * \param t_ time elapsed
        */
        void step_projectiles(double t_) {
            // this could be optimized by sorting by least remaining time, etc.
            for (unsigned i = 0; i < projectiles.size(); i++) {
                auto* e = projectiles[i];

                e->r_collision_grace -= t_;
                e->can_collide = e->r_collision_grace < 0.0;
                e->r_flight_time -= t_;

                if (e->r_flight_time < 0.0) {

                    // kill only if not already dead
                    if (e->health > 0.0) {kill_mobile(e);};

                    projectiles.erase(projectiles.begin()+i);
                    delete e; // we can be sure the projectile is no longer needed
                    
                    i--; // reduce i by one since array size has changed
                }
            }
        }
        /**
         * \brief Updates weapon cooldowns
         * \param t_ time elapsed
        */
        void update_weapon_cooldowns(HostileEntity* e_, double t_) {
            for (auto& w : e_->weapons) {w.current_cooldown -= t_;}
        }
        /**
         * \brief React to player input
        */
        void player_input() {
            // force of player input
            double motor_force = 3000.0;

            // mouse input
            auto mouse_pos = scene.mapPixelToCoords(sf::Mouse::getPosition(scene));
            mouse_pos.x -= player_ship->position.position.x;
            mouse_pos.y -= player_ship->position.position.y;

            // force direction as unit vector
            auto fuv = render::normalize<double>({mouse_pos.x, mouse_pos.y});
            player_ship->resultant_force.x = fuv[0] * motor_force;
            player_ship->resultant_force.y = fuv[1] * motor_force;

            // left and right cannons
            if (sf::Mouse::isButtonPressed(sf::Mouse::Left)) {
                std::cout << "Try fire: Left\n";
                if (player_ship->try_fire(0)) {
                    // allocate memory for a new cannonball
                    auto* new_cannonball = new Projectile(cannonball);
                    add_entity(new_cannonball);
                    launch_projectile(player_ship, 0, new_cannonball);
                }
            }
            if (sf::Mouse::isButtonPressed(sf::Mouse::Right)) {
                std::cout << "Try fire: Right\n";
                if (player_ship->try_fire(1)) {
                    // allocate memory for a new cannonball
                    auto* new_cannonball = new Projectile(cannonball);
                    add_entity(new_cannonball);
                    launch_projectile(player_ship, 1, new_cannonball);
                }
            }
        }
        public:
        GameHandler(render::Scene& scene_, game::Stage& stage_) :
        scene(scene_),
        stage(stage_),
        cannonball(scene_, rs::Vector2(4,4))
        {
            load_assets();
            cannonball.set_object_file(assets::object_file::cannonball);
        }
        /**
         * \brief Prepare for new round
        */
        void initialise_round() {
            reset_entities();
            spawn_ships();

            std::cout << "Created " << mobiles.size() << " mobile entities.\n";
        }
        /**
         * \brief Step all continuous processes
         * \param time_elapsed Time since last frame
        */
        void step(double time_elapsed_) {

            scene.set_camera(player_ship->position.position);
            
            player_input();

            handle_collisions();
            move_mobiles(time_elapsed_);
            update_textures();
            remove_killed_mobiles();

            step_projectiles(time_elapsed_);
            update_weapon_cooldowns(player_ship, time_elapsed_);
        }
    };
}
/**
 * \brief Load stage parameters
 * \param values_ JSON values to load
*/
void load_stage_json(const Json::Value& data_, game::Stage& stage_) {
    noise::Settings settings;
    settings.load_json(data_["noise_parameters"]);
    stage_.set_settings(settings);
    
    Json::Value size = data_["size"];
    stage_.set_stage_size(size["x"].asUInt(), size["y"].asUInt());
}
/**
 * \brief Load island parameters. Stage parameters must be loaded separately
 * \param values_ JSON values to load
*/
void load_island_json(const Json::Value& data_, game::Island& island_) {

    island_.grass_layer_offset = data_["grass_layer_offset"].asDouble();
    island_.cliff_threshold = data_["cliff_threshold"].asDouble();
    island_.visual_height_multiplier = data_["visual_height_multiplier"].asDouble();

    island_.am_brightness = data_["ambient_brightness"].asDouble();

    Json::Value light_dir = data_["light_direction"];
    island_.light_direction[0] = light_dir["x"].asDouble();
    island_.light_direction[1] = light_dir["y"].asDouble();
    island_.light_direction[2] = -std::abs(light_dir["z"].asDouble()); // force `z` to be negative

    island_.light_direction = render::normalize(island_.light_direction);

    auto read_color = [] (Json::Value c) {
        sf::Color color;
        color.r = c[0].asUInt();
        color.g = c[1].asUInt();
        color.b = c[2].asUInt();
        color.a = c.size() == 4 ? c[3].asUInt() : 255;
        return color;
    };

    Json::Value colors = data_["colors"];

    island_.grass = read_color(colors["grass"]);
    island_.dirt = read_color(colors["dirt"]);
    island_.cliff = read_color(colors["cliff"]);
    island_.sand = read_color(colors["sand"]);
    island_.water_shallow = read_color(colors["water_shallow"]);
    island_.water_deep = read_color(colors["water_deep"]);
}
/**
 * \brief Load forest parameters
 * \param values_ JSON values to load
*/
void load_forest_json(const Json::Value& data_, game::ForestGenerator& forest_) {
    noise::Settings settings;
    settings.load_json(data_["noise_parameters"]);
    forest_.set_settings(settings);
}

#endif