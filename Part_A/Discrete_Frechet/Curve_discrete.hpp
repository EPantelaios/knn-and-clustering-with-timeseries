#pragma once
#include <iostream>
#include <vector>
#include <list>
#include <cmath>

using namespace std;


typedef struct curve_struct{

    double x;
    double y;

}curve_struct;

class Curve_discrete{

    private:
        string item_id;
        vector<curve_struct> polygonal_curve;
        vector<curve_struct> init_curve;
        int dimension;

    public:

        Curve_discrete(string id);
        ~Curve_discrete();
        string get_item_id();
        int add_curve_coordinate(curve_struct x);
        int erase_coordinate(int index);
        vector<curve_struct> get_curve();
        void print_curves();
        void set_init_curve(vector<curve_struct> curve);
        vector<curve_struct> get_init_curve();
        void set_dimension(int d);
        int get_dimension();
};