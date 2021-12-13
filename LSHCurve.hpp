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
#include "Curve.hpp"
#include "G_hash.hpp"
#include "utils.hpp"

using namespace std;
using namespace std::chrono;


typedef struct nn_info{

    string item_id;
    float distance;
}nn_info;

class LSHCurve{

    private:
        vector<unordered_multimap<int, Point>> hashtables;
        vector<Curve *> origin_curves; //store curves for verification for the minimum discrete Frechet distance
        vector< vector<Curve> > modified_curves; //snap curves L times
        vector< vector<Point> > points_set; //store the converted curves to vectors
        vector<Curve> modified_query_curves;
        vector<Point> query_vectors;
        vector<G_hash> g_hash;
        int L;
        int k;
        int w;
        int dimension;
        int table_size;
        double delta;
        vector<double> t;

    public:

        LSHCurve(vector<Curve *> dataset, int L, int k, int w, int d, double delta);

        ~LSHCurve();
        //insert point to L hashtables and pointset
        int insert_lsh();

        int find_nearest_neighbor(Curve query, nn_info &min, int cnt);

        int NN_bruteforce(Curve query, nn_info &min, int cnt);

        int hash_all_input();

        int snap_curves();

        int remove_duplicates();

        int convert_input_curves_to_vectors();

        Point convert_curve_to_vector(Curve curve);

        int padding();

        int convert_query_to_vector(Curve curve);

        double frechet_distance(vector<curve> p, vector<curve> q);

        int print_query_results(Curve query, string output_file, string algorithm, print_result_info &info);

        int testing_dimensions_pointset();

        int clear_query_ds();
};