#pragma once
#include <iostream>
#include <cstring>
#include <cassert>
#include <vector>

using namespace std;

typedef struct print_result_info{

    double total_time;
    double total_time_true;
    double maf;
    double total_af;

}print_result_info;


int check_arguments(int argc, char **argv, string &input_file, string &conf_file, string &output_file,
                    string &update_method, string &assignment_method, int &complete_param, int &silhouette_param);

bool file_exist(const string &filename);

unsigned long long int euclidean_mod(long long int a, unsigned long long int b);

double euclidean_distance(vector<double> x, vector<double> y);