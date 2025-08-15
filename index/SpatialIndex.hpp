#ifndef INDEX_SPATIALINDEX_HPP_
#define INDEX_SPATIALINDEX_HPP_

#include <vector>
#include<unordered_set>									 ///////Not ORIGINAL IMPLEMENTATION v1///////////////////////
#include<queue>											//SEGMENT BUFFER OF SEGMENT STORED INSIDE BTREE WILL BE//
#include <iostream>									   ////USED FOR INSERTION AND LOOKUP DIRECTLY///////////////
#include <cassert>
#include "../include/stx/btree.h"
#include "Segment.hpp"
#include "ShrinkingConeSegmentation.hpp"
#include "zorder.hpp"

#include <chrono>

#define ADD_ERR(x, error, size) ((x) + (error) >= (size) ? (size) : (x) + (error))
#define SUB_ERR(x, error) ((x) <= (error) ? 0 : ((x) - (error)))

// Function to compute Euclidean distance
double computeDistance(const std::pair<double, double> p1, const std::pair<double, double> p2) {
    return std::sqrt((p1.first - p2.first) * (p1.first - p2.first) + (p1.second - p2.second) * (p1.second - p2.second));
}

// Assuming a key-value pair is represented as a pair<KeyType, PosType>
template <typename KeyType, typename PosType>
int binarySearchForKey(const std::vector<std::pair<KeyType, PosType>>& vec, const KeyType& key, int start_index, int end_index) {
    auto it = std::lower_bound(vec.begin() + start_index, vec.begin() + end_index-1 , std::make_pair(key, PosType()), [](const std::pair<KeyType, PosType>& a, const std::pair<KeyType, PosType>& b) {
        return a.first < b.first;
    });
//	auto sit=vec.begin();
//	auto it = std::lower_bound(vec.begin()+start_index, vec[end_index-1] , std::make_pair<KeyType, PosType>, [](const std::pair<KeyType, PosType>& a, const std::pair<KeyType, PosType>& b) {
//	        return a.first < b.first;
//	    });

    if (it != vec.end() && it->first == key) {
	//if (it->first == key){
        return std::distance(vec.begin(), it);

    } else {
        return -1; // Not found
    }
}

template <typename KeyType, typename PosType, uint64_t Error = 64, uint64_t BufferSize = 32, typename Floating = long double>
class SpatialIndex
{
	static_assert(Error > 0);
private:
    // Btree structure: key = Z-value (start of the segment), value = Segment
	std::vector<KeyType> SI_keys;
	std::vector<Segment<KeyType, PosType>> segments;
    //stx::btree<KeyType, Segment<KeyType, PosType>> btree;
	stx::btree<KeyType,
    Segment<KeyType, PosType>,
    std::pair<KeyType, Segment<KeyType, PosType>>,
    std::greater<KeyType>,
    stx::btree_default_map_traits<KeyType, Segment<KeyType, PosType>>,
    false,
    std::allocator<std::pair<KeyType, Segment<KeyType, PosType>>>,
    false> btree;
	std::unordered_set<KeyType> deletedKeys;

    // Buffer and error parameters
    double errorBound;
    size_t bufferSize;
    size_t n;

public:
    // Constants for error handling and buffer sizes
    static constexpr uint64_t error_value = Error;
    static constexpr uint64_t seg_error = Error - BufferSize;
    static constexpr uint64_t buffer_size = BufferSize;

    // Define the iterator type using the btree iterator
    using iterator = typename stx::btree<KeyType, Segment<KeyType, PosType>>::const_iterator;

    SpatialIndex(double errorBound, size_t bufferSize)
        : errorBound(errorBound), bufferSize(bufferSize),btree(),segments() {}

    // function to build index for stoing data. Create an object of ShrinkingConeSegmentation and partition data
    void buildIndex(const std::vector<std::pair<double, double>>& points)
    {
        // Convert points to Z-values (keys for the btree)
    	auto start = std::chrono::high_resolution_clock::now();
    	size_t i=0;
    	n=points.size();
       std::vector<std::pair<KeyType, PosType>> index_keys;
        uint64_t zvalue_key;
        for (const auto& point : points)
        {
        	//auto key=zorder::computeZValue(point.first, point.second);

        	zvalue_key=zorder::computeZValue(point.first, point.second);
        	SI_keys.push_back(zvalue_key);
        }

        std::sort(SI_keys.begin(),SI_keys.end());

        for(const auto k:SI_keys){
        	index_keys.push_back(std::make_pair(k,i));
        	i++;
        }

        // Create an instance of ShrinkingConeSegmentation
        ShrinkingConeSegmentation<KeyType, PosType/*, Error, BufferSize*/> segmentation(errorBound/*, bufferSize*/);

        // Vector to hold the segments
        //std::vector<Segment<KeyType, PosType>> segments;

        // Call partitionData to partition the keys and populate the segments vector
        //segmentation.partition_data(keys, bufferSize, segments);
        size_t num_segments=segmentation.partition_data(index_keys, bufferSize, segments);
        //additional code according to fitting_tree.h file to make lowerbound work
        std::vector<std::pair<KeyType, Segment<KeyType, PosType>>> formatted_segments;
		formatted_segments.reserve(num_segments);
		for (auto it = segments.rbegin(); it != segments.rend(); ++it) {
			formatted_segments.emplace_back(it->get_start_key(), *it);
		}
		std::cout<<"Bulk loading the index segments"<<std::endl;
		btree.bulk_load(formatted_segments.begin(), formatted_segments.end());
		//btree.flush();

        // Insert the segments into the btree
//        for (const auto& segment : segments)
//        {
//            btree.insert(segment.get_start_key(), segment);
//        }
		auto end = std::chrono::high_resolution_clock::now();
		//std::chrono::duration<double> duration = end - start;
		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
		std::cout << "index building time: " << duration.count() << " microseconds" << std::endl;

		//// Print the elapsed time in milliseconds
	    //std::cout << "Function execution time: " << duration.count() << " milliseconds" << std::endl;
        std::cout<<"index built successfully"<<std::endl;
        std::cout<<"DReSI size:"<<index_size()<<"B"<<std::endl;
    }

    inline size_t index_size(){
    	size_t segment_size=this->segments.size()*sizeof(this->segments[0]);
    	size_t tree_size=sizeof(this->btree);
    	//size_t insertBuff_size=sizeof(this->);
    	size_t deleteBuff_size=sizeof(this->deletedKeys);
    	return segment_size+tree_size+deleteBuff_size+sizeof(SI_keys);
    }


    	void lookupIndex(const std::vector<std::pair<double, double>>& points){
    		auto start = std::chrono::high_resolution_clock::now();
    		int count=0;
    		for(auto p:points){
    			PosType loc=lookup(p);
    			if(loc!=-1){
    				//std::cout<<"record successfully found: "<<count++<<std::endl;
    				count++;
    			}
    			else{
    				//std::cout<<"record not found"<<std::endl;
    				//count++;
    			}
    		}
			auto end = std::chrono::high_resolution_clock::now();

			//std::chrono::duration<double> duration = end - start;
			auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
					end - start);
			std::cout<<"working till this point, duration is"<<duration.count()<<std::endl;
			std::cout << "time taken in 100 lookups in microseconds:"
							<< duration.count()<<"Average latency of lookup query in microseconds"<<duration.count()/points.size()<<std::endl;//<<"lookup throughput is"<<(100/duration.count()*1000000) << endl;
    	}
    // Modified lookup function to handle deleted keys
        PosType lookup(const std::pair<double, double>& queryPoint) const {
        	if (n == 0){
        		std::cout<<"Index is empty"<<std::endl;
        		return -1;
        	}
            KeyType queryKey = zorder::computeZValue(queryPoint.first, queryPoint.second);

            // Check if the key is marked as deleted
            if (deletedKeys.find(queryKey) != deletedKeys.end()) {
                std::cout << "Key is deleted: " << queryKey << std::endl;
                return -1;
            }

            // Perform the usual lookup
            auto it = btree.lower_bound(queryKey);
            //auto it = btree.find(queryKey);
            if (it == btree.end()) {
                return -1;
            }
            //size_t found=-1;
            const Segment<KeyType, PosType> segment=it.data();
//            for (const Segment<KeyType, PosType> s : segments) {
//            	if (s.get_start_key() == segment.get_start_key()) {
//            		if (s.get_buffer().count(queryKey)>0) {
//            			found=s.get_buffer().at(queryKey);
//            		}
//            	}
//            }
            if(segment.get_buffer().count(queryKey)>0){
            	return segment.get_buffer().at(queryKey);
            }
//            if(found!=-1){
//            	return found;
//            }
            else{
			// Find the position of the key in the segment
            	KeyType startKey = it.data().get_start_key();
            	auto [slope, intercept] = it.data().get_slope_intercept();
            	auto pos = (queryKey - startKey) * slope + intercept;

            	uint64_t hi = ADD_ERR(pos, Error, n);
            	uint64_t lo = SUB_ERR(pos, Error);
//            	return binarySearchForKey(it.data().get_segment_keys(), queryKey,
//					0, it.data().get_segment_keys().size());
            	if(lo<hi)
            		return binarySearchForKey(it.data().get_segment_keys(), queryKey,
            						lo, hi);
            	else
            		return -1;
            }
            std::cout<<"lookup successful"<<std::endl;
        }
    ///need to implement insert, update and delete functions and also to modify lookup to consider
    //that there can be deleted items as well as updated items, also segment buffer should also be checked

    bool insertData(const std::pair<double, double>& dataPoint){
    	// Compute the Z-order value of the query point
    	//auto start = std::chrono::high_resolution_clock::now();
    	    	    KeyType key = zorder::computeZValue(dataPoint.first, dataPoint.second);
    	    	   // std::cout<<"to be inserted key: "<<key<<std::endl;
    	    	    return insertKey(key);

    }
    bool insertKey(const KeyType& key) {
        /// Step 1: Check if the key is marked as deleted
        if (deletedKeys.find(key) != deletedKeys.end()) {
            // Remove the key from the deletedKeys set
            deletedKeys.erase(key);

            // No need to proceed further, as the key is already in the index
           // std::cout << "Key " << key << " was marked as deleted but is already in the index. Marker removed." << std::endl;
            return true;
        }
        // Step 2: Find the appropriate segment using the btree
        auto segmentIt = btree.lower_bound(key);
		if (segmentIt == btree.end()) {
			std::cerr << "No appropriate segment found for key: " << key
					<< std::endl;
			return false;
		}



        //index_keys.emplace_back(key, pos);

        // Step 3: Check if the key is already in the index
        // Reference to the found segment
		Segment<KeyType, PosType> &segment = segmentIt->second;
		if (segment.containsKey(key)) {
			//std::cout << "Key " << key
					//<< " is already present in the index. No insertion performed."
					//<< std::endl;
			return true;
		}



		// Step 4a: Add the key to the end of index_keys
    	SI_keys.push_back(key);
		n++;
        // Step 4: Handle the two cases
        if (segment.insertToBuffer(key,n-1 )) {
            // Case 1: The segment's buffer is not full
            //std::cout << "Key " << key << " added to buffer of segment starting at " << segment.get_start_key() << std::endl;
//            std::vector<std::pair<KeyType, Segment<KeyType, PosType>>> formatted_segments;
//            formatted_segments.emplace_back(segment.get_start_key(), segment);
//            btree.erase(segment.get_start_key());
//            btree.insert(formatted_segments.begin(),formatted_segments.end());

        } else {
            // Case 2: The segment's buffer is full, perform segmentation
            std::cout << "Segment buffer full. Performing segmentation..." << std::endl;

            // Merge keys and buffer and include the new key
            std::vector<std::pair<KeyType, PosType>> mergedKeys = segment.mergeSegmentKeysAndBuffer(key, n-1);
            // Use ShrinkingConeSegmentation to split the merged keys into new segments
            //use segment error instead of errorbound this time
            //More formally, given a specified error of error , we transparently set the error threshold for the
            //segmentation process to (error − buff ). This ensures that a
            //lookup operation will satisfy the specified error even if the element is located in the buffer.
            ShrinkingConeSegmentation<KeyType, PosType> segmentation(errorBound);
            std::vector<Segment<KeyType, PosType>> newSegments;
            std::vector<std::pair<KeyType, Segment<KeyType, PosType>>> formatted_segments;
            size_t num_segments=segmentation.partition_data(mergedKeys, bufferSize, newSegments);
            formatted_segments.reserve(num_segments);
			for (auto it = newSegments.rbegin(); it != newSegments.rend(); ++it) {
				formatted_segments.emplace_back(it->get_start_key(), *it);
			}

            // Remove the old segment from the btree
            btree.erase(segment.get_start_key());
            btree.insert(formatted_segments.begin(),formatted_segments.end());
//            for (auto& newSegment : formatted_segments) {
//                btree.insert(newSegment.get_start_key(), newSegment);
//                std::cout << "New segment inserted with start key: " << newSegment.get_start_key() << std::endl;
//            }
        }

        return true;
    }


    // Delete function
        bool deleteKey(const std::pair<double, double>& dataPoint) {
            // Add the key to the deleted set
        	KeyType key = zorder::computeZValue(dataPoint.first, dataPoint.second);
        	//std::cout<<"to be deleted key: "<<key<<std::endl;
            deletedKeys.insert(key);
            //std::cout<<"deleted successfully "<<std::endl;
            return true;
        }
        void performUpdateQuery(std::vector<std::pair<std::pair<double, double>,std::pair<double, double>>> queries){
        	auto start = std::chrono::high_resolution_clock::now();
        	for(auto q:queries){
        		updateKey(q.first,q.second);
        	}
			auto end = std::chrono::high_resolution_clock::now();
			auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
					end - start);
			std::cout << "time taken in 50 updates in microseconds:" << duration.count()
					<< "Updation throughput "
					<< (queries.size() * 1000000) / duration.count() << std::endl;
        }
        // Update function
        bool updateKey(const std::pair<double, double>& oldDataPoint,const std::pair<double, double>& newDataPoint) {
            // Check if the oldKey exists in the structure
        	KeyType oldKey = zorder::computeZValue(oldDataPoint.first, oldDataPoint.second);
        	KeyType newKey = zorder::computeZValue(newDataPoint.first, newDataPoint.second);
            if (deletedKeys.count(oldKey) > 0) {
                std::cerr << "Error: Old key is already marked as deleted and hence doesn't exist" << std::endl;
                return false;
            }

            // Mark the old key as deleted
            deleteKey(oldDataPoint);

            // Insert the new key
            return insertKey(newKey);
         //   std::cout<<"updated successfully"<<std::endl;
        }

        void performRangeQuery(std::vector<std::pair<std::pair<double, double>,std::pair<double, double>>> queries){
        	auto start = std::chrono::high_resolution_clock::now();
        	for (auto q : queries) {
				rangeQuery(q.first, q.second);
			}
			auto end = std::chrono::high_resolution_clock::now();
			auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
					end - start);
			std::cout << "time taken in 50 range query in microseconds:"
					<< duration.count()
					<< "Average latency of range query in microseconds"
					<< duration.count() / queries.size() << std::endl;
        }
//        bool static comp(const std::pair<KeyType,PosType>&pp1,const std::pair<KeyType,PosType>&pp2){
//        	return pp1.first<pp2.first;
//        }
    //Another requirement is to implement range query, and knn query from the modified z-order curve based spatial data
    //transformation
        std::vector<std::pair<double, double>> rangeQuery(
            const std::pair<double, double>& lowerLeft,
            const std::pair<double, double>& upperRight) const {

            if (n == 0) {
                std::cout << "Index is empty" << std::endl;
                return {};
            }

            // Compute Z-order values for bounding box corners
            KeyType zMin = zorder::computeZValue(lowerLeft.first, lowerLeft.second);
            KeyType zMax = zorder::computeZValue(upperRight.first, upperRight.second);

            // Step 1: Retrieve candidate segments using the Z-order range
            // This loop iterates through segments in the B-tree whose starting key is <= zMax.
            // Depending on how segments are structured (e.g., if a segment can span a large Z-range),
            // you might need a more sophisticated check here to ensure the segment's Z-range
            // overlaps with [zMin, zMax]. Assuming the B-tree key is the segment's min Z-value,
            // this initial scan is a reasonable starting point.
            std::vector<Segment<KeyType, PosType>> candidateSegments;
            for (auto it = btree.lower_bound(zMin); it != btree.end() && it.key() <= zMax; ++it) {
                 // Optional: Add a check here if the segment's max Z-value is >= zMin
                 // if (it.data().get_max_z() >= zMin) {
                     candidateSegments.push_back(it.data());
                 // }
            }

            std::vector<std::pair<double, double>> results;

            // Step 2: Post-process each segment to refine results
            auto comp = [](const std::pair<KeyType, PosType>& a, const std::pair<KeyType, PosType>& b) {
                return a.first < b.first;
            };

            for(const auto& segment : candidateSegments){ // Use const reference to avoid copying segments
                const auto& segmentKeys = segment.get_segment_keys();

                // Find the first key in the segment >= zMin
                auto sit = std::lower_bound(segmentKeys.begin(), segmentKeys.end(), std::make_pair(zMin, 0), comp);

                // Iterate from sit, but STOP when the Z-key EXCEEDS zMax
                for (auto sk_it = sit; sk_it != segmentKeys.end(); ++sk_it) {
                    // If the current Z-key is greater than zMax, we can stop processing this segment
                    if (sk_it->first > zMax) {
                        break;
                    }

                    // The point's Z-value is within [zMin, zMax], now perform the rectangular check
                    auto point = zorder::decodeZValue(sk_it->first);

                    // Check if the point is within the original query rectangle
                    if (point.first >= lowerLeft.first
                        && point.first <= upperRight.first
                        && point.second >= lowerLeft.second
                        && point.second <= upperRight.second) {
                        results.emplace_back(point);
                    }
                }

                // Processing the buffer is still commented out as per your original code.
                // If you re-enable this, consider optimizing buffer access/search as well.
                // const auto& buffer = segment.get_buffer();
                // ... buffer processing logic ...
            }

            return results;
        }

        std::vector<std::pair<double, double>> kNNQuery(
            const std::pair<double, double>& queryPoint,
            int k,
            double alphaX,
            double alphaY) const {

            using DistPointPair = std::pair<double, std::pair<double, double>>;
            auto compare = [](const DistPointPair& a, const DistPointPair& b) {
                return a.first < b.first; // Min-heap for distances
            };
            std::priority_queue<DistPointPair, std::vector<DistPointPair>, decltype(compare)> knnHeap(compare);

            std::unordered_set<KeyType> visitedSegments;

            // Count total points in the dataset
//            size_t totalPoints = 0;
//            for (const auto& segment : segments) {
//                totalPoints += segment.get_segment_keys().size();
//            }
//
//            // Handle case where dataset is empty
//            if (totalPoints == 0) {
//                std::cout << "Dataset is empty. Returning an empty result." << std::endl;
//                return {};
//            }

            // Handle case where total points are less than k
            //if (totalPoints < static_cast<size_t>(k)) {
            if (n < static_cast<size_t>(k)) {
                std::cout << "Insufficient points in the dataset. Returning all available points." << std::endl;

                std::vector<std::pair<double, double>> allPoints;
                for (const auto& segment : segments) {
                    const auto& points = segment.get_segment_keys();
                    for (const auto& point : points) {
                        if (deletedKeys.count(point.first) == 0) {
                            allPoints.push_back(zorder::decodeZValue(point.first));
                        }
                    }
                }
//                for(const auto& k:btree.lower_bound(zorder::decodeZValue(queryPoint))){
//
//                }
                return allPoints;
            }

            // Initialize search window
//            double width = alphaX * std::max(1.0, k / static_cast<double>(std::max(segments.size(), size_t(1))));
//            double height = alphaY * std::max(1.0, k / static_cast<double>(std::max(segments.size(), size_t(1))));
            int width = alphaX * std::max(2.0, k / static_cast<double>(std::max(segments.size(), size_t(1))));
            int height = alphaY * std::max(2.0, k / static_cast<double>(std::max(segments.size(), size_t(1))));
            KeyType queryKey = zorder::computeZValue(queryPoint.first, queryPoint.second);

            size_t processedPoints = 0; // Track processed points
            const double maxWindowSize = 1e6; // Safeguard for infinite expansion

            while (true) {
                // Print debug info
               /* std::cout << "Loop Iteration: width = " << width
                          << ", height = " << height
                          << ", knnHeap.size() = " << knnHeap.size() << std::endl;*/

                // Termination safeguard
                if (width > maxWindowSize || height > maxWindowSize) {
                	goto label;
                    std::cout << "Search window exceeded max size. Exiting loop." << std::endl;


                }

                // Determine Z-order range for the current window
                KeyType zMin = zorder::computeZValue(queryPoint.first - width, queryPoint.second - height);
                KeyType zMax = zorder::computeZValue(queryPoint.first + width, queryPoint.second + height);

                // Retrieve candidate segments
                for (auto it = btree.lower_bound(zMin); it != btree.end() && it.data().get_start_key() <= zMax; ++it) {
                    KeyType segmentStartKey = it.data().get_start_key();
                    if (visitedSegments.count(segmentStartKey) > 0) continue;
                    visitedSegments.insert(segmentStartKey);

                    const auto& segment = it.data();
                    const auto& points = segment.get_segment_keys();
                    for (const auto& point : points) {
                        if (deletedKeys.count(point.first) > 0) continue;

                        std::pair<double, double> qpoint = zorder::decodeZValue(point.first);
                        double dist = computeDistance(queryPoint, qpoint);

                        if (knnHeap.size() < k) {
                            knnHeap.emplace(dist, qpoint);
                        } else if (dist < knnHeap.top().first) {
                            knnHeap.pop();
                            knnHeap.emplace(dist, qpoint);
                        }
                        processedPoints++;
                    }
                }

                // Termination conditions
               // if (knnHeap.size() < k && processedPoints == totalPoints) {
                if (knnHeap.size() < k && processedPoints == n) {
                    //std::cout << "All points processed, but heap size < k. Exiting loop." << std::endl;
                    break;
                } else if (knnHeap.size() == k && knnHeap.top().first > width * width + height * height) {
                   // std::cout << "Termination condition met. Breaking loop. Number of neighbours found"<<knnHeap.size() << std::endl;
                    break;
                } else {
                    width *= 2;
                    height *= 2;
                }
            }

            label:
            // Extract results from the heap
            std::vector<std::pair<double, double>> result;
            while (!knnHeap.empty()) {
                result.push_back(knnHeap.top().second);
                knnHeap.pop();
            }
            //std::cout << "KNN query executed successfully" << std::endl;
            return result;
        }





};

#endif /* INDEX_SPATIALINDEX_HPP_ */
