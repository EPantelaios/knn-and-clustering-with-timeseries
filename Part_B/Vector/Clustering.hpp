#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <chrono>
#include "Hypercube.hpp"
#include "../Common/utils.hpp"

using namespace std;
using namespace std::chrono;


class Clustering{

    private:
        vector<centroid_info> centroids; //K Centroids
        vector<Point_frechet*> points_set;
        int total_clusters;
        string method;     //Classic (Lloyds), LSH, Hypercube
        int L;
        int k;
        int w;
        int max_number_hypercube;
        int hypercube_dimensions;
        int probes;


    public:
        Clustering(vector<Point_frechet*> dataset, int total_clusters, string method, int L, int k, int w, int max_number_hypercube,
                   int hypercube_dimensions, int probes);

        //Implement initialization++ (k-Means++)
        int initialization();

        int Lloyd_assignment();

        LSH* LSH_reverse_assignment(int L, int k, int w);

        Hypercube* Hypercube_reverse_assignment(int k,int M,int probes,int w);

        int update_means();

        //For each centroid if all coordinates have distance less than epsilon from the previous iteration
        //return false and stop the process
        bool centroids_not_same_position(vector<Point_frechet*> means, float epsilon);

        double Cluster_Silhouette(vector<double> &s_i, double &stotal);

        double Silhouette(double a,double b);

        int lower_bound_distance(vector<Point_frechet *> points_set, int total_size, float value);

        double calculate_min_radius(vector<centroid_info> centroids);

        int clear_points();

        void print_results(int complete_param, string output_file);
};