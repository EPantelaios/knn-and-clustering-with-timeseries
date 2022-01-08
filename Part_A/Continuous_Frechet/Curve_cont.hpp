#pragma once
#include <iostream>
#include <vector>
#include <list>
#include <cmath>

using namespace std;

class Curve_cont{

    private:
        string item_id;
        vector<double> polygonal_curve;
        vector<double> filter_curve;
        int dimension;

    public:

        Curve_cont();
        Curve_cont(string id);
        ~Curve_cont();
        string get_item_id();
        int add_curve_coordinate(double x);
        int set_curve_coordinate(int index, double c);
        int erase_coordinate(int index);
        vector<double> get_curve();
        vector<double> *get_curve_address();
        void print_curves();
        void print_filtered_curves();
        void set_filter_curve(vector<double> curve);
        vector<double> get_filtered_curve();
        void set_dimension(int d);
        int get_dimension();
};