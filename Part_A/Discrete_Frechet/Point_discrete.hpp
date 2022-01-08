#pragma once
#include <iostream>
#include <vector>
#include <list>
#include <cmath>
#include "Curve_discrete.hpp"

using namespace std;


class Point_discrete{

    private:
        string item_id;
        vector<double> coordinates;
        int dimension;
        unsigned long long int ID;
        vector<curve_struct> origin_curve;       //coordinated of the initial curve

    public:

        Point_discrete(string id);
        ~Point_discrete();
        int add_coordinate(double x);
        string get_item_id();
        vector<double> get_coordinates();
        void print_coordinates();
        void print_origin_curve();
        void set_dimension(int d);
        int get_dimension();
        void set_ID(unsigned long long int id);
        unsigned long long int get_ID();
        vector<curve_struct> get_origin_curve();
        void set_origin_curve(vector<curve_struct> curve);

};