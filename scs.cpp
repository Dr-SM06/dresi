/*
 * scs.cpp
 *
 *  Created on: 16-Nov-2024
 *      Author: akhilesh
 */
#include "Index/ShrinkingConeSegmentation.hpp"
#include "Index/Segment.hpp"
#include <iostream>
#include <vector>

//int main() {
//    // Sample data points (x, y)
////    std::vector<std::pair<int, double>> data = {{1, 2.0}, {2, 3.0}, {3, 4.0}, {4, 7.0}, {5, 10.0}};
////
////    // Create the segmentation model with an error bound of 1.0
////    ShrinkingConeSegmentation<int, size_t> segmentation(1.0);
////
////    // Vector to store the resulting segments
////    std::vector<Segment<int, size_t>> segments;
////
////    // Partition the data into segments
////    size_t num_segments = segmentation.partition_data(data, 10, segments);
////
////    // Output the resulting segments
////    std::cout << "Number of segments: " << num_segments << std::endl;
////    for (const auto& segment : segments) {
////        std::cout << "Segment: Start Key = " << segment.startKey
////                  << ", End Key = " << segment.endKey
////                  << ", Slope = " << segment.slope << std::endl;
////    }
//
//    return 0;
//}




