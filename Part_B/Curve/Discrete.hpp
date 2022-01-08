#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include <cstring>
#include <unordered_map>
#include <algorithm>
#include <fstream>
#include <assert.h>
#include <chrono>
#include <iterator>
//#include "../Common/G_hash.hpp"
//#include "../Common/utils.hpp"
#include "../Vector/Clustering.hpp"

using namespace std;
using namespace std::chrono;




typedef struct nn_info{

    string item_id;
    float distance;
}nn_info;


class Discrete{

    private:
        vector<unordered_multimap<int, Point_frechet *>> hashtables;
        vector<Point_frechet *> origin_curves; //store curves for verification for the minimum discrete Frechet distance
        vector<vector<Point_frechet>> modified_curves; //snap curves L times
        vector<vector<Point_frechet *>> points_set; //store the converted curves to vectors

        vector<Point_frechet> modified_query_curves;
        vector<Point_frechet> query_vectors;

        vector<G_hash> g_hash;
        int L;
        int k;
        int w;
        int dimension;
        int table_size;
        double delta;
        vector<double> t;
        double **C;

    public:

        Discrete(vector<Point_frechet *> dataset, int L, int k, int w, int d, double delta);

        ~Discrete();
        //insert point to L hashtables and pointset
        int insert_lsh();

        int find_nearest_neighbor(Point_frechet query, nn_info &min, vector<unsigned long long int> IDs);

        int NN_bruteforce(Point_frechet query, nn_info &min);

        int hash_all_input();

        int snap_curves();

        int remove_duplicates();

        int convert_input_curves_to_vectors();

        Point_frechet *convert_curve_to_vector(Point_frechet curve);

        int padding();

        int convert_query_to_vector(Point_frechet &curve);

        double frechet_distance(vector<curve_struct> p, vector<curve_struct> q);

        int print_query_results(Point_frechet query, string output_file, string algorithm, print_result_info &info);

        int testing_dimensions_pointset();

        int clear_query_ds();

        int range_search_discrete(vector<centroid_info> &centroids, int cur_index, double R, vector<vector<int>> &keys);

        vector<Point_frechet *> * get_points_set_address();

        vector<vector<Point_frechet *>> get_points_set();

        int insert_K_points(vector<centroid_info> centroids, int total_clusters, int index);

        int init_points(vector<vector<int>> keys, int total_clusters);

        int get_dimension();
        
};