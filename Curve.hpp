#pragma once
#include <iostream>
#include <vector>
#include <list>
#include <cmath>

using namespace std;


typedef struct curve{

    double x;
    double y;

}curve;

class Curve{

    private:
        string item_id;
        vector<curve> polygonal_curve;
        vector<curve> init_curve;
        int dimension;

    public:

        Curve(string id);
        ~Curve();
        string get_item_id();
        int add_curve_coordinate(curve x);
        int erase_coordinate(int index);
        vector<curve> get_curve();
        void print_curves();
        void set_init_curve(vector<curve> curve);
        vector<curve> get_init_curve();
        void set_dimension(int d);
        int get_dimension();
};