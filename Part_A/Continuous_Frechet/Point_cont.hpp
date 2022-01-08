#pragma once
#include <iostream>
#include <vector>
#include <list>
#include <cmath>
#include "Curve_cont.hpp"

using namespace std;


class Point_cont{

    private:
        string item_id;
        vector<double> coordinates;
        int dimension;
        unsigned long long int ID;
        vector<double> filtered_curve;       //coordinated of the initial curve

    public:
        Point_cont();
        Point_cont(string id);
        ~Point_cont();
        int add_coordinate(double x);
        string get_item_id();
        vector<double> get_coordinates();
        void print_coordinates();
        void print_filtered_curve();
        void set_dimension(int d);
        int get_dimension();
        void set_ID(unsigned long long int id);
        unsigned long long int get_ID();
        vector<double> get_filtered_curve();
        void set_filtered_curve(vector<double> curve);

};