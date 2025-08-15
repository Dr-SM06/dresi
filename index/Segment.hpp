/*
 * Segment.hpp
 *
 *  Created on: 16-Nov-2024
 *      Author: akhilesh
 */

#ifndef INDEX_SEGMENT_HPP_
#define INDEX_SEGMENT_HPP_

#include <vector>
#include <utility>
#include <cstdint>
#include <iostream>
#include <map>



// Type definitions for consistency
//using KeyType = uint64_t; // Z-values
//using PosType = size_t;   // Position index
template<typename KeyType,typename PosType >
class Segment {
public:
	//struct DataItem;
private:
    KeyType startKey;			// The smallest key in the segment
    PosType startPos;			// The position of the smallest key
    double slope;				// The slope of the segment
    KeyType endKey;				// The largest key in the segment
    std::vector<std::pair<KeyType,PosType>> keys; // Stores all the key value pairs present in the segment
    //std::vector<std::pair<KeyType, PosType>> buffer;
    std::map<KeyType,
                 PosType,
                 std::less<KeyType>> buffer;               // A buffer maintained in sorted order for inserts
    uint64_t buffer_size;     // Current Buffer size
    uint64_t max_buffer_size; // Maximum Buffer size allowed
public:
    Segment() = default;
    // Constructor
    Segment(KeyType startKey, PosType startPos, double slope, KeyType endKey, size_t bufferSize)
        : startKey(startKey), startPos(startPos), endKey(endKey), slope(slope), buffer(), buffer_size(0), max_buffer_size(bufferSize){
    	//:startKey(startKey), startPos(startPos), endKey(endKey), slope(slope), buffer_size(0), max_buffer_size(bufferSize){
    	//buffer.reserve(bufferSize);
    }
    Segment(KeyType startKey, PosType startPos, KeyType endKey, double slope,
                std::vector<std::pair<KeyType, PosType>>& seg_keys, uint64_t bufferSize)
            : startKey(startKey), startPos(startPos), endKey(endKey), slope(slope), buffer(),buffer_size(0), max_buffer_size(bufferSize){
    	//: startKey(startKey), startPos(startPos), endKey(endKey), slope(slope),buffer_size(0), max_buffer_size(bufferSize){
            //buffer.reserve(bufferSize);
            keys.reserve(seg_keys.size());
            for (auto i = 0; i < seg_keys.size(); i++){
            	keys.emplace_back(seg_keys[i].first, seg_keys[i].second);
            }
    }

    KeyType get_start_key() const
        {
            return startKey;
        }

    std::vector<std::pair<KeyType,PosType>> get_segment_keys() const{
    	return keys;
    }

    std::pair<long double, long double> get_slope_intercept() const
        {
            return {slope, startPos};
        }
    bool insertToBuffer(const KeyType &key, const PosType &pos){
    	if (buffer_size >= max_buffer_size)
    		return false;
    	this->buffer.insert({key, pos});
    	this->buffer_size += 1;
    	return true;
    }
    //function to merge keys and buffer of the segment together and insert the new key, in the case of buffer overflow
    std::vector<std::pair<KeyType, PosType>> mergeSegmentKeysAndBuffer(const KeyType &new_key, const PosType &new_pos) {
        // Result vector to store merged keys and buffer
    	buffer.insert({new_key,new_pos});
        std::vector<std::pair<KeyType, PosType>> mergedKeys;
        mergedKeys.reserve(keys.size() + buffer.size());

        // Iterators for keys vector and buffer map
        auto keysIt = keys.begin();
        auto bufferIt = buffer.begin();

        // Merge the two sorted containers
        while (keysIt != keys.end() && bufferIt != buffer.end()) {
            if (keysIt->first <= bufferIt->first) {
                mergedKeys.emplace_back(*keysIt);
                ++keysIt;
            } else {
                mergedKeys.emplace_back(bufferIt->first, bufferIt->second);
                ++bufferIt;
            }
        }

        // Add remaining elements from keys, if any
        while (keysIt != keys.end()) {
            mergedKeys.emplace_back(*keysIt);
            ++keysIt;
        }

        // Add remaining elements from buffer, if any
        while (bufferIt != buffer.end()) {
            mergedKeys.emplace_back(bufferIt->first, bufferIt->second);
            ++bufferIt;
        }

        // Clear the buffer
        buffer.clear();
        buffer_size = 0;

        // Return the merged and sorted vector
        return mergedKeys;
    }

//    // Add a key and position to the buffer
//    void addToBuffer(KeyType key, PosType pos) {
//        buffer.emplace_back(key, pos);
//    }
//
//    // Check if the buffer is full
//    bool isBufferFull() const {
//        return buffer.size() == buffer.capacity();
//    }
//
//    // Merge buffer with the segment's initial data
//    void mergeBuffer(std::vector<std::pair<KeyType, PosType>>& segmentData) {
//        segmentData.insert(segmentData.end(), buffer.begin(), buffer.end());
//        buffer.clear();
//    }
    //check keys and buffer to see if key is present or not
    bool containsKey(const KeyType& key) const{
    	// Check the vector first
    	    for (const auto& pair : keys) {
    	        if (pair.first == key) {
    	            return true;
    	        }
    	    }

    	    // Check the map
    	    return buffer.count(key) > 0;
    }
    std::map<KeyType,PosType> get_buffer() const{
    	//std::cout<<"The size of the buffer is: "<<buffer.size()<<std::endl;
    	return buffer;
    }
};





#endif /* INDEX_SEGMENT_HPP_ */
