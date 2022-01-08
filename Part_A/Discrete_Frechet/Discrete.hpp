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
#include "Curve_discrete.hpp"
#include "../Common/G_hash.hpp"
#include "../Common/utils.hpp"

using namespace std;
using namespace std::chrono;


class Discrete{

    private:
        vector<unordered_multimap<int, Point_discrete>> hashtables;
        vector<Curve_discrete *> origin_curves; //store curves for verification for the minimum discrete Frechet distance
        vector< vector<Curve_discrete> > modified_curves; //snap curves L times
        vector< vector<Point_discrete> > points_set; //store the converted curves to vectors
        vector<Curve_discrete> modified_query_curves;
        vector<Point_discrete> query_vectors;
        vector<G_hash> g_hash;
        double **C; //two dimension array for computing discrete frechet distance
        int L;
        int k;
        int w;
        int dimension;
        int table_size;
        double delta;
        vector<double> t;

    public:

        Discrete(vector<Curve_discrete *> dataset, int L, int k, int w, int d, double delta);

        ~Discrete();
        //insert point to L hashtables and pointset
        vector<Point_discrete> get_points_set(int index);
        
        int insert_lsh();

        int find_nearest_neighbor(Curve_discrete query, nn_info &min, vector<unsigned long long int> IDs);

        int NN_bruteforce(Curve_discrete query, nn_info &min);

        int hash_all_input();

        int snap_curves();

        int remove_duplicates();

        int convert_input_curves_to_vectors();

        Point_discrete convert_curve_to_vector(Curve_discrete curve);

        int padding();

        int convert_query_to_vector(Curve_discrete &curve);

        double frechet_distance(vector<curve_struct> p, vector<curve_struct> q);

        int print_query_results(Curve_discrete query, string output_file, string algorithm, print_result_info &info);

        int testing_dimensions_pointset();

        int clear_query_ds();
};