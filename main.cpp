/*
 * main.cpp
 *
 *  Created on: 13-Dec-2024
 *      Author: Supriya Mishra
 */
#include<iostream>
#include<fstream>
#include <sstream>
#include <vector>
#include <random>

#include "Index/SpatialIndex.hpp"


using namespace std;

using Point = pair<double, double>;
using RangeQuery = pair<Point, Point>;

vector<pair<double, double>> readDataPoints(vector<pair<double, double>>& dataPoints, string filename, size_t N){
	    std::ifstream file(filename);
	    std::string line;
	    int rowCount = 1;
	    pair<double, double>p;

	    string lat, lon;
		if (!file.is_open()) {
			std::cerr << "Error opening file!" << std::endl;
			//return -1; // Return -1 if file can't be opened
		}

		// Read each line and count
		while (getline(file, line) && rowCount<=N) {
			stringstream ss(line);
			getline(ss, lat, ','); // Read first column
			getline(ss, lon, ','); // Read second column
			try {
				double value1 = stod(lat); // Convert first column to double
				double value2 = stod(lon); // Convert second column to double

				// Print the values or store them as needed
				//std::cout << "Column 1: " << value1 << ", Column 2: " << value2
						//<< std::endl;
				dataPoints.push_back(make_pair(value1,value2));
			} catch (const std::invalid_argument &e) {
				std::cerr << "Invalid data in CSV file: " << e.what() << std::endl;
			} catch (const std::out_of_range &e) {
				std::cerr << "Value out of range in CSV file: " << e.what()
						<< std::endl;
			}
			rowCount++;
		}
		if(file.eof()){
			cout<<"dataset finished"<<endl;
		}
		cout<<"Records entered:"<<rowCount-1<<endl;
		file.close();
		return dataPoints;
}


// Function to generate sample range queries
vector<RangeQuery> sample_range_queries(const vector<Point>& points, int numQueries = 10) {
	const vector<double> selectivities = {0.001, 0.01, 0.05, 0.1};
	int c=0;
	if (points.empty()) {
        return {};
    }

//    random_device rd;
//    mt19937 gen(rd());

    double minX = points[0].first, maxX = points[0].first;
    double minY = points[0].second, maxY = points[0].second;

    for (const auto& p : points) {
        minX = min(minX, p.first);
        maxX = max(maxX, p.first);
        minY = min(minY, p.second);
        maxY = max(maxY, p.second);
    }

    vector<RangeQuery> rangeQueries;
	for (double selectivity : selectivities) {
		for (int i = 0; i < numQueries; ++i) {
			// Simplified box generation (uniform distribution approximation)
			double width = (maxX - minX) * sqrt(selectivity); // Adjust sqrt for 2D
			double height = (maxY - minY) * sqrt(selectivity);

			// Generate bottom-left corner within data bounds
//			uniform_real_distribution<> xDist(minX, maxX - width);
//			uniform_real_distribution<> yDist(minY, maxY - height);

			Point bottomLeft = points[c];
			//{xDist(gen), yDist(gen)};
			Point topRight = { bottomLeft.first + width, bottomLeft.second
					+ height };

			// Ensure top-right is within bounds (important!)
			topRight.first = (int)min(topRight.first, maxX);
			topRight.second = (int)min(topRight.second, maxY);

			rangeQueries.push_back( { bottomLeft, topRight });
			c++;
		}
	}

    return rangeQueries;
}

vector<pair<double, double>>sample_point_queries(vector<pair<double, double>>& points, size_t s=100){
	// seed the generator
//		std::random_device rd;
//	    std::mt19937 gen(rd());
//	    std::uniform_int_distribution<> uint_dist(0, points.size()-1);
		// generate random indices
		std::vector<pair<double, double>> samples;
		samples.reserve(s);
		for (size_t i = 0; i < s; ++i) {
			//auto idx = uint_dist(gen);
			//samples.emplace_back(points[idx]);
			samples.emplace_back(points[i]);
		}

		return samples;
}



int main(int argc, char **argv){

	string op=argv[1]; //operation code
	string dataset=argv[2]; //dataset file
	size_t N = std::stoi(argv[3]); //dataset size
	uint64_t error=stoi(argv[4]); //error parameter
	uint64_t buf=stoi(argv[5]); //buffer parameter
	string queryFile; //query dataset

	vector<pair<double, double>> dataPoints;
	vector<pair<double, double>> queryPoints;
	readDataPoints(dataPoints,dataset,N);
	//build the index
	SpatialIndex<uint64_t, size_t> si(error,buf);
	si.buildIndex(dataPoints);

	//manual query point
	vector<pair<double,double>> mq={{3574,16219},{19550,19870}};

//	si.rangeQuery(mq.front(),mq.back());
	auto point_queries = sample_point_queries(dataPoints,50);
	auto range_queries = sample_range_queries(dataPoints);
	//auto knn_queries = sample_knn_queries(dataPoints);

	if(op.compare("point")==0){
		//auto start = std::chrono::high_resolution_clock::now();
		si.lookupIndex(point_queries);
		//auto end = std::chrono::high_resolution_clock::now();
		//auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
			//	end - start);
		//cout << "time taken in 100 lookups in microseconds:"
				//<< duration.count()<<"Average latency of lookup query in microseconds"<<duration.count()/point_queries.size()<<endl;//<<"lookup throughput is"<<(100/duration.count()*1000000) << endl;
	}
	if(op.compare("range")==0){
		//si.performRangeQuery(range_queries);
		auto start = std::chrono::high_resolution_clock::now();
		for(auto q:range_queries){
			si.rangeQuery(q.first,q.second);
		}
		auto end = std::chrono::high_resolution_clock::now();
				auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
						end - start);
				//cout<<"total number of points in 50 range queries"<<qpoints<<endl;
				cout << "time taken in 50 range query in microseconds:"
						<< duration.count() <<"Average latency of range query in microseconds"<<duration.count()/range_queries.size()<< endl;

	}

	if(op.compare("insert")==0){
		queryFile=argv[6];
		readDataPoints(queryPoints,queryFile,500);
		auto insert_queries = sample_point_queries(queryPoints,50);
		auto start = std::chrono::high_resolution_clock::now();
		for (auto p : insert_queries) {
			si.insertData(p);
		}
		auto end = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
				end - start);
		cout << "time taken in 50 inserts in microseconds:"
				<< duration.count()<<"Insertion throughput " <<(insert_queries.size()*1000000)/duration.count()<< endl;
	}
	if(op.compare("delete")==0){
		//queryFile=argv[6];
		//readDataPoints(queryPoints,queryFile,500);
		auto delete_queries = sample_point_queries(dataPoints,50);
		auto start = std::chrono::high_resolution_clock::now();
		for (auto p : delete_queries) {
			si.deleteKey(p);
		}
		auto end = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
				end - start);
		cout << "time taken in 50 deletions in microseconds:"
				<< duration.count()<<"Deletion throughput"<<(delete_queries.size()*1000000)/duration.count()<< endl;

	}

	if(op.compare("update")==0){
			queryFile=argv[6];
			readDataPoints(queryPoints,queryFile,500);
			auto old_points = sample_point_queries(dataPoints,50);
			auto new_points =sample_point_queries(queryPoints,50);
			vector<pair<pair<double, double>,pair<double, double>>>updateQueries;

			for (int i=0;i<50;i++) {
				updateQueries.emplace_back(old_points[i],new_points[i]);
			}
			auto start = std::chrono::high_resolution_clock::now();
			si.performUpdateQuery(updateQueries);
			auto end = std::chrono::high_resolution_clock::now();
					auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
							end - start);
			cout << "time taken in 50 updates in microseconds:"
							<< duration.count()<<"Update throughput"<<(updateQueries.size()*1000000)/duration.count()<< endl;
		}

	if(op.compare("knn")==0){
			auto start = std::chrono::high_resolution_clock::now();
			for(auto q:point_queries){
				si.kNNQuery(q,5,1.0,1.0);
			}
			auto end = std::chrono::high_resolution_clock::now();
			auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
							end - start);
			cout << "time taken in knn queries in microseconds:"
									<< duration.count() <<"Average latency of KNN query in microseconds"<<duration.count()/point_queries.size()<< endl;
		}

//	string filename1 = "/home/akhilesh/Downloads/current_projects/uniform10k.csv";
//	string filename2="/home/akhilesh/Downloads/current_projects/uniform1k.csv";
//	uint64_t error=128;
//	uint64_t buf=64;
//	int n=10;
//	random_device rd;
//	mt19937 gen(rd());
//	uniform_real_distribution<> dis(0.0, 1.0);
	//vector<pair<double, double>> queryPoints;
	//pair<double, double>iQuery;
	//cout<<"Error and buffer size respectively"<<error<<","<<buf<<endl;
	//SpatialIndex<uint64_t, size_t> si(error,buf);//si(32.0, 10);
//	si.buildIndex(readDatasetForBulkLoad(filename1));
//	queryPoints=readDatasetForQuery(filename2,1000);
//	auto start = std::chrono::high_resolution_clock::now();
//	for(auto p:queryPoints){
//		si.insertData(p);
//	}
//	auto end = std::chrono::high_resolution_clock::now();
//	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
//	cout<<"time taken in 1000 inserts in milliseconds:"<< duration.count()<<endl;
//	//cout<<"index built successfully"<<endl;
//	start = std::chrono::high_resolution_clock::now();
//	si.lookupIndex(queryPoints);
//	end = std::chrono::high_resolution_clock::now();
//	duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
//	cout<<"time taken in 1000 lookups in milliseconds:"<< duration.count()<<endl;
////	while(n<=1000){
//		si.lookupIndex(readDatasetForLookup(filename,n));
//		n=n*10;
//	}
//	auto start = std::chrono::high_resolution_clock::now();
//	for(int i=0;i<1000;i++){
//		iquery={dis(gen), dis(gen)};
//		si.insertData(iquery);
//
//	}
//	auto end = std::chrono::high_resolution_clock::now();
//
//				//std::chrono::duration<double> duration = end - start;
//				auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
//						end - start);
//				cout<<"time taken in 10 inserts in milliseconds:"<< duration.count()<<endl;


	return 0;
}


