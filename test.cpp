/*
 * test.cpp
 *
 *  Created on: 30-Nov-2024
 *      Author: akhilesh
 */


#include <iostream>
#include <vector>

//#include "Index/SpatialIndex.hpp"

//int main()
//{
//	// Sample data points (x, y)
//	std::vector<std::pair<double, double>> data = { { 1, 2.0 }, { 2, 3.0 }, { 3,
//			4.0 }, { 4, 7.0 }, { 5, 10.0 } };
//	SpatialIndex<uint64_t, size_t> si(32.0, 16);
//	si.buildIndex(data);
//	std::pair<double, double> q;// = std::make_pair(3, 4.0);
//	int loc;// = si.lookup(q);
////	if (loc == -1)
////		std::cout << "key not found" << std::endl;
////	else
////		std::cout << "key found at position:" << loc << std::endl;
////	for(auto q:data){
////		loc=si.lookup(q);
////		if (loc == -1)
////				std::cout << "key not found" << std::endl;
////			else
////				std::cout << "key found at position:" << loc << std::endl;
////	}
//
////	q = std::make_pair(2.5, 4.5);
////	if (si.insertData(q)) {
////		std::cout << "Insert Successful" << std::endl;
////	} else {
////		std::cout << "Insert Unsuccessful" << std::endl;
////	}
////
////	loc = si.lookup(q);
////	if (loc == -1)
////		std::cout << "key not found" << std::endl;
////	else
////		std::cout << "key found at position:" << loc << std::endl;
//
//	for(auto q:data){
//			loc=si.lookup(q);
//			if (loc == -1)
//					std::cout << "key not found" << std::endl;
//				else
//					std::cout << "key found at position:" << loc << std::endl;
//		}
//		q = std::make_pair(2.5, 4.5);
//		if (si.insertData(q)) {
//			std::cout << "Insert Successful" << std::endl;
//		} else {
//			std::cout << "Insert Unsuccessful" << std::endl;
//		}
//		q = std::make_pair(1, 2.0);
//		si.deleteKey(q);
//		si.updateKey(q,std::make_pair(1,2.0));
//		si.updateKey(std::make_pair(2,3.0),std::make_pair(1,2.0));
//	for(auto q:data){
//			loc=si.lookup(q);
//			if (loc == -1)
//					std::cout << "key not found" << std::endl;
//				else
//					std::cout << "key found at position:" << loc << std::endl;
//		}
//
////	auto rq=si.rangeQuery(std::make_pair(1, 2.0),std::make_pair(5, 10.0));
////
////	std::cout<<"Range query executed successfully, search result count"<<rq.size()<<std::endl;
//	//q=std::make_pair(1, 2.0);
//
////	rq=si.kNNQuery(q,2,1.0,1.0);
////
////	std::cout<<"KNN query executed successfully, search result count"<<rq.size()<<std::endl;
//
//	///index built and lookup performed.. further function would be insert, delete and update
//	return 0;
//}
