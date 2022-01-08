#pragma once
#include <iostream>
#include <vector>
#include <map>
#include <cmath>

using namespace std;


class Point_hc{

    private:
        string item_id;
        vector<double> coordinates;
        int dimension;

    public:

        Point_hc(string id);
        ~Point_hc();
        int add_coordinate(double x);
        string get_item_id();
        vector<double> get_coordinates();
        double get_a_coordinate(int i);
        void print_coordinates();
        void set_dimension(int d);
        int get_dimension();
};