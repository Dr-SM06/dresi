#ifndef INDEX_SHRINKINGCONESEGMENTATION_HPP_
#define INDEX_SHRINKINGCONESEGMENTATION_HPP_

#include <vector>
#include <iostream>
#include <cmath>
#include <stdexcept>
#include "Segment.hpp"

template<typename KeyType, typename PosType>
class ShrinkingConeSegmentation {
private:
    struct Point {
        KeyType x;
        double y;
    };

    struct Slope {
        double dx;
        double dy;

        bool operator<(const Slope& p) const {
            return dy * p.dx < dx * p.dy;
        }

        bool operator>(const Slope& p) const {
            return dy * p.dx > dx * p.dy;
        }

        explicit operator double() const {
            return dy / dx;
        }
    };

    const double error;
    Point first_point;
    Point last_point;
    Slope lower_slope = {1, 0};
    Slope upper_slope = {0, 1};
    size_t points_in_segment = 0;

public:
    explicit ShrinkingConeSegmentation(double error) : error(error) {
        if (error < 0) {
            throw std::invalid_argument("Error can't be less than zero");
        }
    }

    bool add_point(const KeyType& x, double y) {
        Point current_point{x, y};
        Point p1{x, y + error};
        Point p2{x, y - error};

        if (points_in_segment == 0) {
            first_point = current_point;
            last_point = current_point;
            lower_slope = {1, 0};
            upper_slope = {0, 1};
            ++points_in_segment;
            return true;
        }

        if (points_in_segment == 1) {
            lower_slope = {x - first_point.x, y - first_point.y - error};
            upper_slope = {x - first_point.x, y - first_point.y + error};
            ++points_in_segment;
            last_point = current_point;
            return true;
        }

        Slope slope = {x - first_point.x, y - first_point.y};

        bool outside_lower_slope = slope < lower_slope;
        bool outside_upper_slope = slope > upper_slope;

        if (outside_lower_slope || outside_upper_slope) {
            points_in_segment = 0;
            return false;
        }

        if (p1.x - first_point.x < upper_slope.dx) {
            upper_slope = {x - first_point.x, y - first_point.y + error};
        }

        if (p2.x - first_point.x > lower_slope.dx) {
            lower_slope = {x - first_point.x, y - first_point.y - error};
        }

        last_point = current_point;
        ++points_in_segment;
        return true;
    }

    Segment<KeyType, PosType> get_segment(std::vector<std::pair<KeyType, PosType>>& keys, uint64_t buf_size) {
        if (points_in_segment == 1) {
            return Segment<KeyType, PosType>(first_point.x, first_point.y, last_point.x, 0.0, keys, buf_size);
        }

        double u_slope = static_cast<double>(upper_slope);
        double l_slope = static_cast<double>(lower_slope);
        double slope = (u_slope + l_slope) / 2;
        return Segment<KeyType, PosType>(first_point.x, first_point.y, last_point.x, slope, keys, buf_size);
    }

    size_t partition_data(const std::vector<std::pair<KeyType, PosType>>& data, uint64_t buf_size, std::vector<Segment<KeyType, PosType>>& segments) {
        if (data.empty()) {
            return 0;
        }

        std::vector<std::pair<KeyType, PosType>> keys;
        size_t num_segments = 0;
        size_t start = 0;

        // Start the segmentation
        for (const auto& kv : data) {
            if (!add_point(kv.first, kv.second)) {
                segments.push_back(get_segment(keys, buf_size));
                start = keys.size();
                ++num_segments;
                keys.clear(); // clear the buffer for the next segment
            }
            keys.push_back(kv);
        }

        // Add the last segment if any points were collected
        if (!keys.empty()) {
            segments.push_back(get_segment(keys, buf_size));
            ++num_segments;
        }

        return num_segments;
    }
};

#endif /* INDEX_SHRINKINGCONESEGMENTATION_HPP_ */
