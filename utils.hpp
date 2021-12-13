#pragma once
#include <iostream>
#include <cstring>
#include <cassert>
#include <vector>
#include <cmath>
#include "Curve.hpp"

using namespace std;


typedef struct print_result_info{

    double total_time;
    double total_time_true;
    double maf;

}print_result_info;

int check_arguments(int argc, char **argv, string &input_file, string &query_file, string &output_file, int &k,
                    int &L, int &M, int &probes, string &algorithm, string &metric, double &delta);

bool file_exist(const string &filename);

unsigned long long int euclidean_mod(long long int a, unsigned long long int b);

double euclidean_distance(vector<double> x, vector<double> y);