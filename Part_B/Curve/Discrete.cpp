#include "Discrete.hpp"

Discrete::Discrete(vector<Point_frechet *> dataset, int L, int k, int w, int d, double delta){

    this->origin_curves = dataset;
    this->L=L;
    this->k=k;
    this->w=w;
    this->dimension=d;
    this->delta=delta;
    this->table_size= origin_curves.size()/16;

    C = new double*[2*d];
    for(int i=0;i<2*d;i++){

    C[i] = new double[2*d];
    }

    //this->table_size=5;

    if(table_size==0)
        table_size=1;

    for(int i=0;i<L;i++){
        
        G_hash g(k, w, dimension*2, this->table_size);
        g_hash.push_back(g);

        unordered_multimap<int, Point_frechet *> hashtable;
        hashtables.push_back(hashtable);

        vector<Point_frechet> point;
        modified_curves.push_back(point);

        vector<Point_frechet *> point_ptr;
        points_set.push_back(point_ptr);

        Point_frechet init_curve("");
        modified_query_curves.push_back(init_curve);
    }

    random_device rand_dev;
    mt19937 gen(rand_dev());
    uniform_real_distribution<float> x(0, delta);
    
    int max_iter = L*2;
    for(int i=0;i<max_iter;i++){

        t.push_back(x(gen));
    }
}

Discrete::~Discrete(){

    for(int i=0;i<2*dimension;i++){

        delete [] C[i];
    }
    delete [] C;
    for(vector<Point_frechet *>::iterator it=origin_curves.begin(); it != origin_curves.end(); it++){
    
        delete *it;
    }
    for(int i=0;i<L;i++){

        int cnt=0;
        for(vector<Point_frechet*>::iterator it = points_set[i].begin(); it != points_set[i].end(); it++){

            delete *it;
            cnt++;
        }
    }

}


vector<Point_frechet *> * Discrete::get_points_set_address(){

    return &points_set[0];
}

vector<vector<Point_frechet *>> Discrete::get_points_set(){

    return points_set;
}




//insert point to L hashtables and pointset
int Discrete::insert_lsh(){

    unsigned long long int key, size;

    for(int i=0; i<L; i++){

        size = points_set[i].size();
        for(int j=0;j<size;j++){

            key = g_hash[i].G_compute(*points_set[i][j]);
            hashtables[i].insert(make_pair(key, points_set[0][j]));
        }
    }

    return 0;
}





int Discrete::range_search_discrete(vector<centroid_info> &centroids, int cur_index, double R, vector<vector<int>> &keys){

    int key=0, cnt_new_points=0;
    double tmp_distance=0.0;


    for(int i=0; i<L; i++){
        
        
        if(keys[cur_index].size()<L){
            key = g_hash[i].G_compute(query_vectors[i]);

            keys[cur_index].push_back(key);
        }

        auto its = hashtables[i].equal_range(keys[cur_index][i]);

        for(auto it = its.first; it != its.second; ++it){

            //Skip this iteration. Point has already been assigned to some cluster
            //cout << "marked: " << it->second.get_marked() << endl;
            if(it->second->get_marked() == cur_index){

                continue;
            }

            tmp_distance = frechet_distance(centroids[cur_index].centroid->get_curve(), it->second->get_curve());

            if(tmp_distance < R){
                
                if(it->second->get_marked() != -1){
                    
                    //if two centroids have the same point, choose the one with smaller distance
                    if(tmp_distance < it->second->get_cur_distance()){

                        it->second->set_cur_distance(tmp_distance);
                        int tmp_index = it->second->get_marked();
                        it->second->set_marked(cur_index);/////check check

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





int Discrete::insert_K_points(vector<centroid_info> centroids, int total_clusters, int index){

    unsigned long long int key;

    for(int i=0; i<L; i++){


        key = g_hash[i].G_compute(query_vectors[i]);
        hashtables[i].insert(make_pair(key, centroids[index].centroid));
    }

    return 0;
}



int Discrete::init_points(vector<vector<int>> keys, int total_clusters){

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




int Discrete::hash_all_input(){

    assert(!snap_curves());

    assert(!remove_duplicates());

    assert(!convert_input_curves_to_vectors());
    assert(!padding());
    testing_dimensions_pointset();
    assert(!insert_lsh());
    
    return 0;
}



int Discrete::snap_curves(){

    int cnt_curve=0;
    double x, y;
    curve_struct new_curve;


    for(int i=0;i<L;i++){

        cnt_curve=0;

        for(vector<Point_frechet *>::iterator it = origin_curves.begin(); it != origin_curves.end(); it++){

            Point_frechet Y((*it)->get_item_id());
            modified_curves[i].push_back(Y);

            for(int j=0;j<dimension;j++){

                x = (*it)->get_curve()[j].x;
                y = (*it)->get_curve()[j].y;

                new_curve.x = floor(abs(x - t[i]) / delta + 1.0/2.0) * delta + t[i];
                new_curve.y = floor(abs(y - t[i+1]) / delta + 1.0/2.0) * delta + t[i+1];

                modified_curves[i][cnt_curve].add_curve_coordinate(new_curve);
            }
            
            //store the the initial curve before proceeding with the conversions
            vector<curve_struct> tmp_curve = (*it)->get_curve();
            modified_curves[i][cnt_curve].set_curve(tmp_curve);

            cnt_curve++;
        }

    }

    return 0;
}



//Remove consecutive duplicates (x,y)
int Discrete::remove_duplicates(){

    int cnt_curve=0, index;
    double x, y, next_x, next_y;
    int cur_index=0, max_index=0;

    for(int i=0;i<L;i++){
        
        int cnt=0;
        for(vector<Point_frechet>::iterator it = modified_curves[i].begin(); it != modified_curves[i].end(); it++){
            
            int tmp_cnt=0;

            cur_index=0;
            max_index=it->get_dimension();

            while(cur_index<max_index-1){
                
                x = it->get_curve()[cur_index].x;
                y = it->get_curve()[cur_index].y;
                next_x = it->get_curve()[cur_index+1].x;
                next_y = it->get_curve()[cur_index+1].y;
                
                if(x == next_x && y == next_y){
                    
                    it->erase_curve_coordinate(cur_index);

                    max_index--;
                    tmp_cnt++;
                }
                else{

                    cur_index++;
                }
            }

            //cout << "L = " << i+1 << "\tAfairethikan " << tmp_cnt << " stoixeia\n";
        }

    }

    return 0;
}




//Concatenate all input curves
int Discrete::convert_input_curves_to_vectors(){

    for(int i=0;i<L;i++){

        int cnt=0;
        for(vector<Point_frechet>::iterator it = modified_curves[i].begin(); it != modified_curves[i].end(); it++){

            Point_frechet *new_point = convert_curve_to_vector(*it);
            points_set[i].push_back(new_point);
            cnt++;
        }
    }

    return 0;
}



//Concatenate curve
Point_frechet *Discrete::convert_curve_to_vector(Point_frechet curve){

    Point_frechet *new_point = new Point_frechet(curve.get_item_id());

    for(int i=0;i<curve.get_curve().size();i++){

        new_point->add_coordinate(curve.get_curve()[i].x);
        new_point->add_coordinate(curve.get_curve()[i].y);
    }

    new_point->set_curve(curve.get_curve());

    return new_point;
}






int Discrete::get_dimension()
{
    return dimension;
}




int Discrete::padding(){

    int cur_size=0, size = dimension * 2;

    for(int i=0;i<L;i++){

        for(vector<Point_frechet *>::iterator it = points_set[i].begin(); it != points_set[i].end(); it++){

            cur_size = (*it)->get_dimension();
            while(cur_size < size){

                (*it)->add_coordinate(100000);
                (*it)->add_coordinate(100000);
                cur_size = cur_size + 2;
            }
        }
    }

    return 0;
}




int Discrete::testing_dimensions_pointset(){

    int size = dimension * 2;

    for(int i=0;i<L;i++){

        //if(i==0) cout << "\nSize: " << size << endl;

        for(vector<Point_frechet *>::iterator it = points_set[i].begin(); it != points_set[i].end(); it++){

            if((*it)->get_dimension() != size){

                assert(1==0);
            }
        }

    }

    return 0;
}







int Discrete::convert_query_to_vector(Point_frechet &query){

    double x, y, next_x, next_y;
    curve_struct new_curve;
    int cnt_curve=0, index;
    int cur_index=0, max_index=0;
    //cout << "This is the query size " << query.get_curve().size() << endl;
    //cout << " and this is the dimension " << dimension <<endl;
    clear_query_ds();

    for(int i=0;i<L;i++){
        
        //snapping
        for(int j=0;j<dimension;j++){
            x = query.get_curve()[j].x;
            y = query.get_curve()[j].y;

            new_curve.x = floor(abs(x - t[i]) / delta + 1.0/2.0) * delta + t[i];
            new_curve.y = floor(abs(y - t[i+1]) / delta + 1.0/2.0) * delta + t[i+1];

            modified_query_curves[i].add_curve_coordinate(new_curve);
        }

        //store the the initial curve before proceeding with the conversions
        vector<curve_struct> tmp_curve = query.get_curve();
        modified_query_curves[i].set_curve(tmp_curve);




        cur_index=0;
        max_index=modified_query_curves[i].get_dimension();
        //Remove consecutive duplicates (x,y)
        while(cur_index<max_index-1){
            
            x = modified_query_curves[i].get_curve()[cur_index].x;
            y = modified_query_curves[i].get_curve()[cur_index].y;
            next_x = modified_query_curves[i].get_curve()[cur_index+1].x;
            next_y = modified_query_curves[i].get_curve()[cur_index+1].y;
            
            if(x == next_x && y == next_y){

                modified_query_curves[i].erase_curve_coordinate(cur_index);
                max_index--;
            }
            else{

                cur_index++;
            }

        }


        //Convert curve to vector
        Point_frechet *new_point = convert_curve_to_vector(modified_query_curves[i]);
        query_vectors.push_back(*new_point);

        //this is probably important
        //delete new_point;

        int cur_size=0, size=dimension * 2;

        cur_size = query_vectors[i].get_dimension();
        //cout << "cur_size: " << cur_size << endl;
        while(cur_size < size){

            query_vectors[i].add_coordinate(100000);
            query_vectors[i].add_coordinate(100000);
            cur_size = cur_size + 2;
        }
    }

    //cout << "QUERY VECTOR:\n";

    //query_vectors[0].print_coordinates();

    //cout << endl << endl<<endl;



    return 0;
}



//clear Data Structure for query
int Discrete::clear_query_ds(){

    query_vectors.clear();
    modified_query_curves.clear();
    
    for(int i=0;i<L;i++){

        Point_frechet init_curve("");
        modified_query_curves.push_back(init_curve);
    }
    return 0;
}






inline double euclidean_distance_curve(curve_struct p, curve_struct q){

    //sqrt( (x1-x2)^2 + (y1-y2)^2 )
    double x = p.x - q.x;
    double y = p.y - q.y;
    return sqrt(x*x + y*y);
}




double Discrete::frechet_distance(vector<curve_struct> p, vector<curve_struct> q){

    int m1 = p.size(), m2 = q.size();
    double euclidean_distane, current_index, left, up, diagonal;
    double final_value;

    for(int i=0;i<m1;i++){

        for(int j=0;j<m2;j++){

            if(i==0 && j==0){
                
                C[i][j] = euclidean_distance_curve(p[i], q[j]);
                continue;
            }
            else if(i==0 && j>0){

                current_index = C[i][j-1];
            }
            else if(i>0 && j==0){

                current_index = C[i-1][j];
            }
            else if(i>0 && j>0){
                
                left = C[i-1][j];
                up = C[i][j-1];
                diagonal = C[i-1][j-1];
                
                if(left <= up && left <= diagonal) current_index=left;
                else if(up <=left && up <= diagonal) current_index=up;
                else if(diagonal <= left && diagonal <= up) current_index=diagonal;
            }
            else{

                assert(1==0);
            }

            euclidean_distane = euclidean_distance_curve(p[i], q[j]);

            if(euclidean_distane > current_index){

                C[i][j]=euclidean_distane;
            }
            else{

                C[i][j]=current_index;
            }
        }
    }
    final_value=C[m1-1][m2-1];

    return final_value;
}









/*int Discrete::print_query_results(Curve query, string output_file, string algorithm, print_result_info &info){

    ofstream out;
    double current_maf;
    nn_info nn_distance;
    nn_info bruteforce_distance;
    double time, time_true;
    vector<unsigned long long int> IDs;

    convert_query_to_vector(query); 

    for(int i=0;i<L;i++){

        IDs.push_back(g_hash[i].ID_compute(query_vectors[i]));
    }

    auto time_1 = high_resolution_clock::now();
    find_nearest_neighbor(query, nn_distance, IDs);
    auto time_2 = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(time_2 - time_1);
    time = duration.count();
    info.total_time += time;
   

    auto time_true_1 = high_resolution_clock::now();
    NN_bruteforce(query, bruteforce_distance);
    auto time_true_2 = high_resolution_clock::now();
    auto duration_true = duration_cast<milliseconds>(time_true_2 - time_true_1);
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

    return 0;
} */
