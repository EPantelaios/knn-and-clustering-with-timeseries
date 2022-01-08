#include <iostream>
#include <vector>
#include <cmath>
#include <cstring>
#include <unordered_map>
#include <algorithm>
#include "../Common/utils.hpp"
#include "LSH.hpp"


LSH::LSH(vector<Point_frechet *> dataset, int L, int k, int w, int d){

    this->points_set = dataset;
    this->L=L;
    this->k=k;
    this->w=w;
    this->dimension=d;
    this->table_size= points_set.size()/8;
    
    if(table_size==0)
        table_size=1;

    for(int i=0;i<L;i++){
        
        G_hash g(k, w, dimension, this->table_size);
        g_hash.push_back(g);
        unordered_multimap<int, Point_frechet *> hashtable;
        hashtables.push_back(hashtable);
    }
}

LSH::~LSH(){

    cout << "Successfully Destroyed LSH" << endl;
    for(int i=0; i<L; i++){
        
        hashtables[i].clear();
    }

    for(vector<Point_frechet *>::iterator it=points_set.begin(); it != points_set.end(); it++){
    
        delete *it;
    }
}

//insert point to L hashtables and pointset
int LSH::insert_lsh(){

    int key;
    int size = points_set.size();
    for(int i=0; i<L; i++){
        
        for(int j=0;j<size;j++){
            
            key = g_hash[i].G_compute(*points_set[j]);
            hashtables[i].insert(make_pair(key, points_set[j]));
        }
    }

    return 0;
}


int LSH::insert_K_points(vector<centroid_info> centroids, int total_clusters, int index){

    int key;
    for(int i=0; i<L; i++){

        key = g_hash[i].G_compute(*centroids[index].centroid);
        hashtables[i].insert(make_pair(key, centroids[index].centroid));

    }

    return 0;
}



int LSH::init_points(vector<vector<int>> keys, int total_clusters){

    for(int i=0;i<total_clusters;i++){

        for(int j=0;j<keys[i].size();j++){

            auto its = hashtables[j].equal_range(keys[i][j]);
            for(auto it = its.first; it != its.second; ++it){

                it->second->set_marked(-1);
                it->second->set_cur_distance(INT32_MAX);
            }
        }
    }
    return 0;
}



int LSH::range_search(vector<centroid_info> &centroids, int cur_index, double R, vector<vector<int>> &keys){

    int key=0, cnt_new_points=0;
    float tmp_distance=0.0;

    for(int i=0; i<L; i++){

        key = g_hash[i].G_compute(*centroids[cur_index].centroid);
        
        if(keys[cur_index].size()<L){

            keys[cur_index].push_back(key);
        }

        auto its = hashtables[i].equal_range(key);

        for(auto it = its.first; it != its.second; ++it){

            //Skip this iteration. Point has already been assigned to some cluster
            if(it->second->get_marked() == cur_index){

                continue;
            }

            tmp_distance = euclidean_distance(centroids[cur_index].centroid->get_coordinates(), it->second->get_coordinates());

            if(tmp_distance < R){
                
                if(it->second->get_marked() != -1){
                    
                    //if two centroids have the same point, choose the one with smaller distance
                    if(tmp_distance < it->second->get_cur_distance()){

                        it->second->set_cur_distance(tmp_distance);
                        int tmp_index = it->second->get_marked();
                        it->second->set_marked(cur_index);
                        centroids[cur_index].points.push_back(it->second);
                        cnt_new_points++;

                        for(int i=0;i<centroids[tmp_index].points.size();i++){

                            if(it->second->get_item_id() == centroids[tmp_index].points[i]->get_item_id()){
                                
                                centroids[tmp_index].points.erase(centroids[tmp_index].points.begin() + i);
                                break;
                            }
                        }
                    }

                }
                else{
                    
                    it->second->set_cur_distance(tmp_distance);
                    it->second->set_marked(cur_index);
                    centroids[cur_index].points.push_back(it->second);
                    cnt_new_points++;
                }

            }
        }
    }

    return cnt_new_points;
}
