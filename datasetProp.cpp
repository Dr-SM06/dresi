/*
 * datasetProp.cpp
 *
 *  Created on: 13-Dec-2024
 *      Author: akhilesh
 */
#include <iostream>
#include <fstream>
#include <string>

//int countRowsInCSV(const std::string& filename) {
//    std::ifstream file(filename);
//    std::string line;
//    int rowCount = 0;
//
//    // Check if file is open
//    if (!file.is_open()) {
//        std::cerr << "Error opening file!" << std::endl;
//        return -1; // Return -1 if file can't be opened
//    }
//
//    // Read each line and count
//    while (std::getline(file, line)) {
//        rowCount++;
//    }
//
//    file.close();
//    return rowCount;
//}
//
//int main() {
//    std::string filename = "/home/akhilesh/Downloads/current_projects/normal_dataset.csv"; // Replace with your CSV file path
//    int rows = countRowsInCSV(filename); Number of rows in the CSV file: 553471
//
//    if (rows != -1) {
//        std::cout << "Number of rows in the CSV file: " << rows << std::endl;
//    }
//
//    return 0;
//}
//vector<pair<double, double>> readDatasetForQuery(string filename,int n){
//	std::ifstream file(filename);
//		    std::string line;
//		    int rowCount = 0;
//		    vector<pair<double, double>> queryPoints;
//		    pair<double, double>p;
//
//		    string lat, lon;
//			if (!file.is_open()) {
//				std::cerr << "Error opening file!" << std::endl;
//				//return -1; // Return -1 if file can't be opened
//			}
//
//			// Read each line and count
//			while (getline(file, line) && rowCount<=n) {
//				stringstream ss(line);
//				getline(ss, lat, ','); // Read first column
//				getline(ss, lon, ','); // Read second column
//				try {
//					double value1 = stod(lat); // Convert first column to double
//					double value2 = stod(lon); // Convert second column to double
//
//					// Print the values or store them as needed
//					//std::cout << "Column 1: " << value1 << ", Column 2: " << value2
//							//<< std::endl;
//					queryPoints.push_back(make_pair(value1,value2));
//				} catch (const std::invalid_argument &e) {
//					std::cerr << "Invalid data in CSV file: " << e.what() << std::endl;
//				} catch (const std::out_of_range &e) {
//					std::cerr << "Value out of range in CSV file: " << e.what()
//							<< std::endl;
//				}
//				rowCount++;
//			}
//			if(file.eof()){
//				cout<<"dataset finished"<<endl;
//			}
//			cout<<"Records searched:"<<rowCount-1<<endl;
//			file.close();
//			return queryPoints;
//}

//vector<pair<double, double>>sample_range_queries(vector<pair<double, double>>& points, size_t s=100){
//
//}
//
//vector<pair<double, double>>sample_knn_queries(vector<pair<double, double>>& points, size_t s=100){
//
//}


