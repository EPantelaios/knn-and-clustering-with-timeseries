#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include <cstring>
#include <unordered_map>
#include <algorithm>
#include <fstream>
#include <chrono>
#include "../Common/G_hash.hpp"
#include "../Common/utils.hpp"

using namespace std;
using namespace std::chrono;

class LSH{

    private:
        vector<unordered_multimap<int, Point_lsh *>> hashtables;
        vector<G_hash> g_hash;
        vector<Point_lsh *> points_set; //store points for verification for the actual distance
        int L;
        int k;
        int w;
        int dimension;
        double R;
        int table_size;

    public:

        LSH(vector<Point_lsh *> dataset, int L, int k, int w, int d, double R);

        ~LSH();
        //insert point to L hashtables and pointset
        int insert_lsh();

        int find_nearest_neighbor(Point_lsh query, nn_info &info, vector<unsigned long long int> IDs);

        int find_N_nearest_neighbor(Point_lsh query, vector<nn_info> &min_vector, int N);

        int NN_bruteforce(Point_lsh query, nn_info &info);

        int kNN_bruteforce(Point_lsh query, int N, vector<nn_info> &min_distances);

        int range_search(Point_lsh query, double R, vector<string> &names);

        void swap(nn_info &x, nn_info &y);

        void bubbleSort(vector<nn_info> &items, int n);

        bool duplicate_point(vector<nn_info> points, string item, int N);

        int N_min_bucket(vector<nn_info> &min_vector, float distance, string item_id, int N);

        void print_query_results(Point_lsh query, string output_file, string algorithm, print_result_info &info);
};