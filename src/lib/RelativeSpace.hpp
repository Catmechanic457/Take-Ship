#ifndef RS_H
#define RS_H

namespace rs {
    /**
     * \brief Represents a point in 2D space
    */
    template <typename T>
    struct Vector2 {
        T x;
        T y;
        Vector2() {}
        template <typename G>
        Vector2(Vector2<G> vect_) : x(vect_.x), y(vect_.y) {}
        Vector2(T x_, T y_) : x(x_), y(y_) {}
        bool operator== (Vector2<T> v_) {
            return x == v_.x and y == v_.y;
        }
        void from_bearing(double d_, double r_) {
            x = cos(r_) * d_;
            y = sin(r_) * d_;
        }
    };
    /**
     * \brief Represents a point and rotation in 2D space
    */
    struct Position {
        Vector2<double> position;
        double rotation;
        Position() {}
        Position(Vector2<double> pos_, double rot_) : position(pos_), rotation(rot_) {}
        Position(double x_, double y_, double rot_) : position(x_, y_), rotation(rot_) {}
    };
}

#endif