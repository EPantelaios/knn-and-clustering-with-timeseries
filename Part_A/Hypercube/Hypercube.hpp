#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include <cstring>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <fstream>
#include <chrono>
#include "../Common/H_hash.hpp"
#include "../Common/utils.hpp"

using namespace std;
using namespace std::chrono;


typedef struct knn_info_hyper{//helper struct for information about knn

    string item_id;
    double distance;
}knn_info_hyper;


class Hypercube
{
    private:
        vector<H_hash> h_hash;
        vector<Point_hc *> points_set; //store points for verification for the actual distance
        unordered_multimap<int, Point_hc *>  hyperc; //the structure from int(key) to the corresponding points
        vector<map<int, int>> f_function;

        int k; // d'
        int N; //#nn def 1
        int M; // max number of points def 10
        int probes; // number of neighbouring buckets to visit def 2
        int R; //range def 10000
        int original_dimensions; 
        int w; //w for the hast function (2-6 by def)

    public:
        Hypercube(int k, int M, int probes, int N, int R, int w, int d,vector <Point_hc *> points_set);

        void hypercube_insert();

        int range(Point_hc query,int R,vector<string> &range_search_names); //typical range search

        int knn(Point_hc query,int N,vector <knn_info_hyper> &point_info); //k nearest neighbour if N=1 then it is basically nn

        int find_bit(int num,int pos,int k);

        int reverse_to_int(string mystring);//string to int

        void hamming_distance_strings(string mystring,int cur_height,int hamming_distance, vector<string> &all_comb,int k,int start_point);

        vector<int> hamming(int l,int k,vector <string> &everything,int distance); //returns a vector with the
        //neighbouring buckets of the hypercube

        int init_points(vector<vector<int>> keys,int centroid_num); 

        void print_query_results(Point_hc query, string output_file, string algorithm, print_result_info &info);

        ~Hypercube();
        //sorting
        void swap(knn_info_hyper &x, knn_info_hyper &y);

        void bubbleSort(vector<knn_info_hyper> &items, int n);
        //bruteforce
        int kNN_bruteforce(Point_hc query, int N, vector<knn_info_hyper> &min_distances); 

        int NN_bruteforce(Point_hc query, knn_info_hyper &info);

};