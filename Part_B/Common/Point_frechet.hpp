#pragma once
#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <cmath>

using namespace std;


typedef struct curve_struct{

    double x;
    double y;

}curve_struct;

class Point_frechet{

    private:
        string item_id;
        vector<curve_struct> origin_curve;    //coordinated of the initial curve
        vector<double> coordinates;
        int dimension;

        //For Lloyd
        float min_distance_from_centroids;
        float square_dist;
        float sum_square_dist;

        //For Reverse Range Search LSH and Hypercube
        int marked;  //For LSH --> -1 = unmarked, 0 = belongs to centroid 0, 1 = belongs to centroid 1, ...
        float cur_distance; //Distance from current centroid

    public:

        Point_frechet(string id);
        int add_coordinate(double x);
        string get_item_id();
        vector<double> get_coordinates();
        double get_a_coordinate(int i);
        void print_coordinates();
        void print_coordinates_output(ofstream &out);
        void print_curves();
        void print_curves_output(ofstream &out);
        void set_dimension(int d);
        int get_dimension();
        void set_min_distance(float i);
        float get_min_distance();
        void set_square_dist(float i);
        float get_square_dist();
        void set_sum_square_dist(float i);
        float get_sum_square_dist();
        float init_sum_square_distance();
        void set_cur_distance(double i);
        double get_cur_distance();
        void set_marked(int i);
        int get_marked();
        curve_struct get_curve_coordinate(int i);
        vector<curve_struct> get_curve();
        void set_curve(vector<curve_struct> curve);
        int add_curve_coordinate(curve_struct c);
        int erase_curve_coordinate(int index);
        vector<curve_struct> * get_curve_address();
};