#ifndef _HYPERCUBE_H_
#define _HYPERCUBE_H_
#include <iostream>
#include <vector>
#include <cmath>
#include <cstring>
#include <unordered_map>
#include <algorithm>
#include <map>
#include "../Common/utils.hpp"
#include "LSH.hpp"


typedef struct knn_info_hyper{

    string item_id;
    int distance;
}knn_info_hyper;


class Hypercube
{
    private:
            vector<H_hash> h_hash;
        vector<Point_frechet *> points_set; //store points for verification for the actual distance
        unordered_multimap<int, Point_frechet *>  hyperc;
        //f_function must be 2d?
        //can be calculated at the beginning
        vector<map<int,int>> f_function; 

        int k; // new dimensions 14
        int N; //no nn 1
        int M; // max number of points 14
        int probes; // number of neighbouring buckets to visit 2
        int R; //range 10000
        int original_dimensions; 
        int w; 

    public:
        Hypercube(int k, int M, int probes, int N, int R, int w, int d,vector <Point_frechet *> points_set);
        void hypercube_insert();
        int range(Point_frechet query,int R,vector <knn_info_hyper> &point_info);
        int knn(Point_frechet query,int N,vector <knn_info_hyper> &point_info);
        int find_bit(int num,int pos,int k);
        int reverse_to_int(string mystring);
        void hamming_distance_strings(string mystring,int cur_height,int hamming_distance,vector<string> &all_comb,int k,int start_point);
        vector<int> hamming(int l,int k,vector <string> &everything,int distance);
        int range_search(vector<centroid_info> &centroids,int i,double R,vector<vector<int>> &keys);
        int init_points(vector<vector<int>> keys,int centroid_num);
        int hypercube_insert_point(vector<centroid_info> centroids, int total_clusters,int index);

        ~Hypercube();

};

#endif
