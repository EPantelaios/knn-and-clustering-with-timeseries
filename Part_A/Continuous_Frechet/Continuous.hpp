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
#include "Curve_cont.hpp"
#include "../Common/G_hash.hpp"
#include "../Common/utils.hpp"
#include "Fred/include/curve.hpp"
#include "Fred/include/frechet.hpp"


#define EPSILON 1  //value for filtering preprocessing of initial curves


using namespace std;
using namespace std::chrono;


class Continuous{

    private:
        unordered_multimap<int, Point_cont> hashtable;
        vector<Curve_cont *> origin_curves; //store curves for verification for the minimum discrete Frechet distance
        vector<Curve_cont> modified_curves; //snapped curves
        vector<Point_cont> points_set; //store the converted curves to vectors
        Curve_cont origin_query_curve;
        Point_cont query_vector;
        G_hash g_hash;
        int k;
        int w;
        int origin_dimension;
        int table_size;
        double delta;
        double t;

    public:

        Continuous(vector<Curve_cont *> dataset, int k, int w, int d, double delta);

        ~Continuous();
        //insert point to L hashtables and pointset
        int insert_lsh();

        int find_nearest_neighbor(Curve query, nn_info &min, unsigned long long int ID);

        int NN_bruteforce(Curve query, nn_info &min);

        int hash_all_input();

        int filtering();

        int snap_curves();

        int minima_maxima();

        int convert_input_curves_to_vectors();

        Point_cont convert_curve_to_vector(Curve_cont curve);

        int padding();

        int convert_query_to_vector(Curve_cont &curve);

        double frechet_distance(vector<double> p, vector<double> q);

        int print_query_results(Curve_cont query, string output_file, string algorithm, print_result_info &info);

        int testing_dimensions_pointset();

};