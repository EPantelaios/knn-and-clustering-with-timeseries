#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <chrono>
#include "Discrete.hpp"

using namespace std;
using namespace std::chrono;


class Clustering_Frechet{

    private:
        vector<centroid_info> centroids; //K Centroids
        vector<Point_frechet *> origin_curves;
        vector<Point_frechet *> *points_set_ptr;
        int total_clusters;
        string update_method;
        string assignment_method; //Classic (Lloyds), LSH Vector, Hypercube, LSH Frechet
        int L;
        int k;
        int dimensions;
        int w;
        int max_number_hypercube;
        int hypercube_dimensions;
        int probes;
        double delta;
        double **C;
    public:

        Clustering_Frechet(vector<Point_frechet *> origin_curves, int total_clusters, string update_method, string assignment_method,
                            int L, int k, int w, int max_number_hypercube, int hypercube_dimensions, int probes, double delta);

        ~Clustering_Frechet();

        //Implement initialization++ (k-Means++)
        int initialization();

        int Lloyd_assignment();
        
        Discrete* LSH_Frechet(int L, int k, int w);

        int update_means();

        //For each centroid if all coordinates have distance less than epsilon from the previous iteration
        //return false and stop the process
        bool centroids_not_same_position(vector<Point_frechet *> means, float epsilon);

        double Cluster_Silhouette(vector<double> &s_i, double &stotal);

        double Silhouette(double a,double b);

        int lower_bound_distance(vector<Point_frechet *> points_set, int total_size, float value);

        double calculate_min_radius(vector<centroid_info> centroids);

        int clear_points();

        void print_results(int complete_param, int silhouette_param, string output_file);

        double frechet_distance_cluster(vector<curve_struct> *p, vector<curve_struct> *q);

};

struct back_track_path
{

    int t;
    int *path;
    double value;
};
back_track_path *  frechet_distance_and_path(vector <curve_struct> *p,vector <curve_struct> * q ); //optimal backtracking_path
vector <curve_struct> find_mean_curve(back_track_path *the_path,vector <curve_struct> p,vector <curve_struct> q); //returns a mean curve between two curves 

struct comp_bin_node
{ //tree node
    comp_bin_node * left_child;
    comp_bin_node * right_child;
    comp_bin_node * parent;
    Point_frechet *value; //this is going to be curve
    int cur_height;

};
comp_bin_node * comp_bin_node_init(Point_frechet * value,int height);
int *find_path(int tree_counter,int &counter);
void comp_bin_add_node(int tree_counter,Point_frechet* value,comp_bin_node * head,int height);


class comp_bin_tree
{//tree for mean curve

    public:
        comp_bin_node *head;
        int tree_counter;
        int height;
        comp_bin_tree(vector <Point_frechet *>  *values);
        vector <curve_struct>  comp_bin_tree_compute(comp_bin_node *the_node,int dimensions);
        void delete_tree(comp_bin_node *the_node);
};


double mean(double a, double b);
int filtering(vector <curve_struct> *mean_curve,int dimensions);
