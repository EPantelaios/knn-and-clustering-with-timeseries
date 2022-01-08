#include <iostream>
#include <vector>
#include <cmath>
#include <cstring>
#include <unordered_map>
#include <algorithm>
#include "../Common/G_hash.hpp"
#include "../Common/utils.hpp"

using namespace std;


typedef struct centroid_info{

    Point_frechet *centroid; //K Centroids;
    vector<Point_frechet *> points; //points which belong to the centroid
    int points_size;

}centroid_info;


typedef struct knn_info{

    string item_id;
    float distance;
}knn_info;



class LSH{

    private:
        vector< unordered_multimap<int, Point_frechet *> > hashtables;
        vector<G_hash> g_hash;
        vector<Point_frechet *> points_set; //store points for verification for the actual distance
        int L;
        int k;
        int w;
        int dimension;
        int table_size;

    public:

        LSH(vector<Point_frechet *> dataset, int L, int k, int w, int d);

        ~LSH();

        //insert point to L hashtables and pointset
        int insert_lsh();

        int insert_K_points(vector<centroid_info> centroids, int total_clusters, int index);

        int init_points(vector<vector<int>> keys, int total_clusters);

        int range_search(vector<centroid_info> &centroids, int cur_index, double R, vector<vector<int>> &keys);

};