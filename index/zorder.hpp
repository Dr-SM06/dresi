/*
 * zorder.hpp
 *
 *  Created on: 16-Nov-2024
 *      Author: akhilesh
 */

#ifndef INDEX_ZORDER_HPP_
#define INDEX_ZORDER_HPP_

#include <cstdint>
#include <cmath>

class zorder {
public:
    static uint64_t computeZValue(double lat, double lon) {
        // Normalize latitude and longitude to [0, 1] range
        auto normalize = [](double value, double min, double max) {
            return (value - min) / (max - min);
        };

        uint32_t x = static_cast<uint32_t>(normalize(lat, -90, 90) * 1e6);
        uint32_t y = static_cast<uint32_t>(normalize(lon, -180, 180) * 1e6);

        return interleaveBits(x, y);
    }
	static std::pair<double, double> decodeZValue(uint64_t z) {
		// Extract x and y by deinterleaving bits
		uint32_t x = deinterleaveBits(z, false);
		uint32_t y = deinterleaveBits(z, true);

		// Denormalize back to latitude and longitude
		auto denormalize = [](uint32_t value, double min, double max) {
			//return min + (static_cast<double>(value) * (max - min)) ;/// 1e6;
			//return min + (static_cast<double>(value) / (max - min)) * 1e6;
			return min + (static_cast<double>(value) / 1e6) * (max - min);
		};

		double lat = denormalize(x, -90, 90);
		double lon = denormalize(y, -180, 180);

		return {lat, lon};
	}

private:
    static uint64_t interleaveBits(uint32_t x, uint32_t y) {
        uint64_t z = 0;
        for (int i = 0; i < 32; ++i) {
            //z |= ((x & (1 << i)) << i) | ((y & (1 << i)) << (i + 1));
        	z |= ((static_cast<uint64_t>(x) & (1ULL << i)) << i) |
        	             ((static_cast<uint64_t>(y) & (1ULL << i)) << (i + 1));
        }
        return z;
    }
	static uint32_t deinterleaveBits(uint64_t z, bool extractY) {
		uint32_t value = 0;
		for (int i = 0; i < 32; ++i) {
			if (extractY) {
				//value |= ((z & (1ULL << (2 * i + 1))) >> (i + 1));
				if ((z & (1ULL << (2 * i + 1)))) {
				                value |= (1 << i);
				            }
			} else {
				//value |= ((z & (1ULL << (2 * i))) >> i);
				if ((z & (1ULL << (2 * i)))) {
				                value |= (1 << i);
				            }
			}
		}
		return value;
	}
};
/////////////////////////////////////////////////////////
//#include <stdexcept>
//#include <utility>
//
//class ZOrder {
//public:
//    static uint64_t computeZValue(double lat, double lon) {
//        if (lat < -90.0 || lat > 90.0 || lon < -180.0 || lon > 180.0) {
//            throw std::invalid_argument("Latitude or longitude out of bounds");
//        }
//
//        auto normalize = [](double value, double min, double max) {
//            return (value - min) / (max - min);
//        };
//
//        uint32_t x = static_cast<uint32_t>(normalize(lat, -90.0, 90.0) * UINT32_MAX);
//        uint32_t y = static_cast<uint32_t>(normalize(lon, -180.0, 180.0) * UINT32_MAX);
//
//        return interleaveBits(x, y);
//    }
//
//    static std::pair<double, double> decodeZValue(uint64_t z) {
//        uint32_t x = deinterleaveBits(z, false);
//        uint32_t y = deinterleaveBits(z, true);
//
//        auto denormalize = [](uint32_t value, double min, double max) {
//            return min + (static_cast<double>(value) / UINT32_MAX) * (max - min);
//        };
//
//        double lat = denormalize(x, -90.0, 90.0);
//        double lon = denormalize(y, -180.0, 180.0);
//
//        return {lat, lon};
//    }
//
//private:
//    static uint64_t interleaveBits(uint32_t x, uint32_t y) {
//        uint64_t z = 0;
//        for (int i = 0; i < 32; ++i) {
//            z |= ((x & (1ULL << i)) << i) | ((y & (1ULL << i)) << (i + 1));
//        }
//        return z;
//    }
//
//    static uint32_t deinterleaveBits(uint64_t z, bool extractY) {
//        uint32_t value = 0;
//        for (int i = 0; i < 32; ++i) {
//            if (extractY) {
//                value |= ((z & (1ULL << (2 * i + 1))) >> (i + 1));
//            } else {
//                value |= ((z & (1ULL << (2 * i))) >> i);
//            }
//        }
//        return value;
//    }
//};



#endif /* INDEX_ZORDER_HPP_ */
