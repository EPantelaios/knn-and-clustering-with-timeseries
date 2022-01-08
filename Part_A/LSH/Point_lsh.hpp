#pragma once
#include <iostream>
#include <vector>
#include <list>
#include <cmath>

using namespace std;

class Point_lsh{

    private:
        string item_id;
        vector<double> coordinates;
        int dimension;
        vector<unsigned long long int> ID;

    public:

        Point_lsh(string id);
        ~Point_lsh();
        int add_coordinate(double x);
        string get_item_id();
        vector<double> get_coordinates();
        void print_coordinates(Point_lsh x);
        void set_dimension(int d);
        int get_dimension();
        void add_ID(unsigned long long int id);
        unsigned long long int get_ID(int i);

};