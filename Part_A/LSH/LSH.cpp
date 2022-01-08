#include "LSH.hpp"


LSH::LSH(vector<Point_lsh *> dataset, int L, int k, int w, int d, double R){

    this->points_set = dataset;
    this->L=L;
    this->k=k;
    this->w=w;
    this->dimension=d;
    this->R=R;
    this->table_size= points_set.size()/16;

    if(table_size==0)
        table_size=1;

    for(int i=0;i<L;i++){
        
        G_hash g(k, w, dimension, this->table_size);
        g_hash.push_back(g);
        unordered_multimap<int, Point_lsh *> hashtable;
        hashtables.push_back(hashtable);
    }
}

LSH::~LSH(){

    for(vector<Point_lsh *>::iterator it=points_set.begin(); it != points_set.end(); it++){
    
        delete *it;
    }
}


//insert point to L hashtables and pointset
int LSH::insert_lsh(){

    int key;
    int size = points_set.size();
    unsigned long long int tmp;
    for(int i=0; i<L; i++){
        
        for(int j=0;j<size;j++){
            
            tmp = g_hash[i].ID_compute(*points_set[j]);
            points_set[j]->add_ID(tmp);
            
            key = tmp % this->table_size;
            hashtables[i].insert(make_pair(key, points_set[j]));

        }
    }

    return 0;
}


int LSH::find_nearest_neighbor(Point_lsh query, nn_info &info, vector<unsigned long long int> IDs){

    info.distance = INT32_MAX;
    int key=0;
    unsigned long long int id;
    double tmp_distance=0.0;
    int flag_ID=1;

    for(int i=0; i<L; i++){

        key = IDs[i] % table_size;
        auto its = hashtables[i].equal_range(key);
        flag_ID=1;

        for (auto it = its.first; it != its.second; ++it){

            if(IDs[i]== it->second->get_ID(i)){
                
                tmp_distance = euclidean_distance(query.get_coordinates(), it->second->get_coordinates());

                if(info.distance > tmp_distance){

                    info.distance  = tmp_distance;
                    info.item_id = it->second->get_item_id();
                    flag_ID=0;
                }
            }
        }

        //only if not find ID(p) == ID(q)
        if(flag_ID == 1){

            for(auto it = its.first; it != its.second; ++it){

                tmp_distance  = euclidean_distance(query.get_coordinates(), it->second->get_coordinates());
                if(info.distance > tmp_distance){

                    info.distance = tmp_distance;
                    info.item_id = it->second->get_item_id();
                }
            }
        }
    }

    return 0;
}



int LSH::find_N_nearest_neighbor(Point_lsh query, vector<nn_info> &min_vector, int N){

    nn_info aux;
    int key=0, cnt=0;
    double tmp_distance=0.0;
    unsigned long long int id;

    for(int i=0; i<L; i++){

        id = g_hash[i].ID_compute(query);
        key = id % table_size;
        auto its = hashtables[i].equal_range(key);

        for(auto it = its.first; it != its.second; ++it){

            tmp_distance = euclidean_distance(query.get_coordinates(), it->second->get_coordinates());
            //only for the first N points
            if(cnt < N){

                aux.distance = tmp_distance;
                aux.item_id = it->second->get_item_id();
                min_vector.push_back(aux);
                cnt++;
            }
            else{

                N_min_bucket(min_vector, tmp_distance, it->second->get_item_id(), min_vector.size());
            }

        }
    }

    bubbleSort(min_vector, min_vector.size());

    return 0;
}



int LSH::NN_bruteforce(Point_lsh query, nn_info &info){
    
    info.distance=INT32_MAX;
    double tmp_distance=0;

    for(int i=0;i<points_set.size();i++){

        tmp_distance = euclidean_distance(points_set[i]->get_coordinates(), query.get_coordinates());
        if(info.distance>tmp_distance){

            info.distance=tmp_distance;
            info.item_id = points_set[i]->get_item_id();
        }
    }

    return 0;
}



int LSH::kNN_bruteforce(Point_lsh query, int N, vector<nn_info> &min_distances){

    nn_info tmp;

    for(int i=0;i<points_set.size();i++){
        
        tmp.distance = euclidean_distance(points_set[i]->get_coordinates(), query.get_coordinates());
        tmp.item_id = points_set[i]->get_item_id();
        min_distances.push_back(tmp);
    }

    bubbleSort(min_distances, min_distances.size());

    return 0;
}


int LSH::range_search(Point_lsh query, double R, vector<string> &names){

    int key=0, flag=1;
    float tmp_distance=0.0;
    unsigned long long int id;
    string tmp_name;

    for(int i=0; i<L; i++){

        id = g_hash[i].ID_compute(query);
        key = id % table_size;
        auto its = hashtables[i].equal_range(key);

        for(auto it = its.first; it != its.second; ++it){

            flag=1;

            tmp_distance = euclidean_distance(query.get_coordinates(), it->second->get_coordinates());

            if(tmp_distance < R){
                
                tmp_name = it->second->get_item_id();
                for(int i=0;i<names.size();i++){

                    if(tmp_name == names[i]){

                        flag=0;
                        break;
                    }
                }

                if(flag==1){

                    names.push_back(it->second->get_item_id());
                }
            }
        }
    }
    return 0;
}




int LSH::N_min_bucket(vector<nn_info> &min_vector, float distance, string item_id, int N){

    bubbleSort(min_vector, N);
    for(int i=N-1;i>=0;i--){

        if(min_vector[i].distance > distance){
            
            if(duplicate_point(min_vector, item_id, N) == false){
                
                min_vector[i].distance = distance;
                min_vector[i].item_id = item_id;
                
            }
            return 0;
        }
    }
    return 1;
}



bool LSH::duplicate_point(vector<nn_info> points, string item, int N){

    for(int i=0;i<N;i++){

        if(points[i].item_id == item)
            return true;
    }

    return false;
}



void LSH::swap(nn_info &x, nn_info &y){ 

    nn_info temp = x; 
    x = y; 
    y = temp; 
}


void LSH::bubbleSort(vector<nn_info> &items, int n){ 

    int i, j, swapped; 

    for(i = 0; i < n-1; i++){

            swapped = 0; 
            for (j = 0; j < n-i-1; j++){ 

                if (items[j].distance > items[j+1].distance){ 
                    
                    swap(items[j], items[j+1]); 
                    swapped = 1; 
                } 
            } 
            if(swapped == 0) 
                break; 
    } 
}



void LSH::print_query_results(Point_lsh query, string output_file, string algorithm, print_result_info &info){
    
    ofstream out;
    double current_maf;
    nn_info nn_distance;
    nn_info bruteforce_distance;
    double time, time_true;
    vector<unsigned long long int> IDs;

    for(int i=0;i<L;i++){

        IDs.push_back(g_hash[i].ID_compute(query));
    }

    auto time_1 = high_resolution_clock::now();
    find_nearest_neighbor(query, nn_distance, IDs);
    auto time_2 = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(time_2 - time_1);
    time = duration.count();
    info.total_time += time;   

    auto time_true_1 = high_resolution_clock::now();
    NN_bruteforce(query, bruteforce_distance);
    auto time_true_2 = high_resolution_clock::now();
    auto duration_true = duration_cast<microseconds>(time_true_2 - time_true_1);
    time_true = duration_true.count();
    info.total_time_true += time_true;


    if(bruteforce_distance.distance != 0){

        current_maf = nn_distance.distance / bruteforce_distance.distance;
        if(current_maf > info.maf){

            info.maf = current_maf;
        }
        info.total_af += current_maf;
    }

    out.open(output_file, ios::app);
    if (!out){

        cout << "File not created!";
    }
    else{

        out << "Query: " << query.get_item_id() << endl;
        out << "Algorithm: " << algorithm << endl;
        out << "Approximate Nearest neighbor: " << nn_distance.item_id << endl;
        out << "True Nearest neighbor: " << bruteforce_distance.item_id << endl;
        out << "distanceApproximate: " << nn_distance.distance << endl;
        out << "distanceTrue: " << bruteforce_distance.distance << endl;
        out << endl;
        out.close(); 
    }
}