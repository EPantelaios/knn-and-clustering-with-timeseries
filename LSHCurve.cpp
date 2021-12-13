#include "LSHCurve.hpp"

LSHCurve::LSHCurve(vector<Curve *> dataset, int L, int k, int w, int d, double delta){

    this->origin_curves = dataset;
    this->L=L;
    this->k=k;
    this->w=w;
    this->dimension=d;
    this->delta=delta;
    this->table_size= origin_curves.size()/16;

    if(table_size==0)
        table_size=1;

    for(int i=0;i<L;i++){
        
        G_hash g(k, w, dimension*2, this->table_size);
        g_hash.push_back(g);

        unordered_multimap<int, Point> hashtable;
        hashtables.push_back(hashtable);

        vector<Curve> curve;
        modified_curves.push_back(curve);
        
        vector<Point> point;
        points_set.push_back(point);

        Curve init_curve("");
        modified_query_curves.push_back(init_curve);
    }

    random_device rand_dev;
    mt19937 gen(rand_dev());
    uniform_real_distribution<float> x(0, delta);
    
    int max_iter = L*2;
    for(int i=0;i<max_iter;i++){

        t.push_back(x(gen));
    }

    int cnt=1;
    for(vector<double>::iterator it=t.begin(); it != t.end(); it++){
    
        cout << "t_"  << cnt << ":    " <<  *it << endl;
        
        cnt++;
    }

}

LSHCurve::~LSHCurve(){

    for(vector<Curve *>::iterator it=origin_curves.begin(); it != origin_curves.end(); it++){
    
        delete *it;
    }

}




//insert point to L hashtables and pointset
int LSHCurve::insert_lsh(){

    int key, size;
    unsigned long long int tmp;


    for(int i=0; i<L; i++){

        size = points_set[i].size();

        for(int j=0;j<size;j++){

            tmp = g_hash[i].ID_compute(points_set[i][j]);
            points_set[i][j].set_ID(tmp);

            key = tmp % this->table_size;
            hashtables[i].insert(make_pair(key, points_set[i][j]));
        }
    }

    return 0;
}



int LSHCurve::find_nearest_neighbor(Curve query, nn_info &min, int cnt){

    min.distance = INT32_MAX;
    int key=0;
    unsigned long long int id;
    double tmp_distance=0.0;
    int flag_ID=1;

    convert_query_to_vector(query); 
    

    for(int i=0; i<L; i++){

        id = g_hash[i].ID_compute(query_vectors[i]);
        key = id % table_size;
        auto its = hashtables[i].equal_range(key);
        flag_ID=1;

        for (auto it = its.first; it != its.second; ++it){

            if(id == it->second.get_ID()){
                
                printf("MPHKE STO ID\n");

                tmp_distance = frechet_distance(query.get_curve(), it->second.get_origin_curve());

                if(min.distance > tmp_distance){
                    min.distance = tmp_distance;
                    min.item_id = it->second.get_item_id();
                    flag_ID = 0;
                }
            }
        }


        //only if not find ID(p) == ID(q)
        if(flag_ID == 1){

            for(auto it = its.first; it != its.second; ++it){

                tmp_distance = frechet_distance(query.get_curve(), it->second.get_origin_curve());

                if(min.distance>tmp_distance){

                    min.distance = tmp_distance;
                    min.item_id=it->second.get_item_id();
                }
            }
        }
    }

    cout << "Approximate Nearest neighbor:\n";
    cout << "Item_id:  " << min.item_id << endl;
    cout << "Distance: " << min.distance << endl;

    return 0;
}




int LSHCurve::NN_bruteforce(Curve query, nn_info &min, int cnt){
    
    min.distance = INT32_MAX;
    double tmp_distance=0.0;
    vector<double> distance777;

    //cout << "ORIGIN CURVES\n\n\n";
    //origin_curves[0]->print_curves();
    
    for(int i=0;i<origin_curves.size();i++){

        tmp_distance = frechet_distance(origin_curves[i]->get_curve(), query.get_curve());

        if(min.distance>tmp_distance){

            min.distance=tmp_distance;
            min.item_id = origin_curves[i]->get_item_id();
        }
    }

   
    cout << "Real Nearest neighbor:\n";
    cout << "Item_id:  " << min.item_id << endl;
    cout << "Distance: " << min.distance << endl;

    return 0;
}






int LSHCurve::hash_all_input(){

    assert(!snap_curves());
    assert(!remove_duplicates());
    assert(!convert_input_curves_to_vectors());
    assert(!padding());
    testing_dimensions_pointset();
    //exit(1);
    assert(!insert_lsh());

    return 0;
}



int LSHCurve::snap_curves(){

    int cnt_curve=0;
    double x, y;
    curve new_curve;


    for(int i=0;i<L;i++){

        cnt_curve=0;

        for(vector<Curve *>::iterator it = origin_curves.begin(); it != origin_curves.end(); it++){

            Curve Y((*it)->get_item_id());
            modified_curves[i].push_back(Y);

            for(int j=0;j<dimension;j++){

                x = (*it)->get_curve()[j].x;
                y = (*it)->get_curve()[j].y;

                new_curve.x = floor(abs(x - t[i]) / delta + 1.0/2.0) * delta + t[i];
                new_curve.y = floor(abs(y - t[i+1]) / delta + 1.0/2.0) * delta + t[i+1];

                modified_curves[i][cnt_curve].add_curve_coordinate(new_curve);
            }
            
            //store the the initial curve before proceeding with the conversions
            vector<curve> tmp_curve = (*it)->get_curve();
            modified_curves[i][cnt_curve].set_init_curve(tmp_curve);

            cnt_curve++;
        }

    }

    return 0;
}



//Remove consecutive duplicates (x,y)
int LSHCurve::remove_duplicates(){

    int cnt_curve=0, index;
    double x, y, next_x, next_y;
    int cur_index=0, max_index=0;

    for(int i=0;i<L;i++){
        
        int cnt=0;
        for(vector<Curve>::iterator it = modified_curves[i].begin(); it != modified_curves[i].end(); it++){
            
            int tmp_cnt=0;

            cur_index=0;
            max_index=it->get_dimension();

            while(cur_index<max_index-1){
                
                x = it->get_curve()[cur_index].x;
                y = it->get_curve()[cur_index].y;
                next_x = it->get_curve()[cur_index+1].x;
                next_y = it->get_curve()[cur_index+1].y;
                
                if(x == next_x && y == next_y){
                    
                    it->erase_coordinate(cur_index);

                    max_index--;
                    tmp_cnt++;
                }
                else{

                    cur_index++;
                }
            }

            cout << "L = " << i+1 << "\tAfairethikan " << tmp_cnt << " stoixeia\n";
        }

    }

    return 0;
}




//Concatenate all input curves
int LSHCurve::convert_input_curves_to_vectors(){

    for(int i=0;i<L;i++){

        for(vector<Curve>::iterator it = modified_curves[i].begin(); it != modified_curves[i].end(); it++){

            Point new_point = convert_curve_to_vector(*it);
            points_set[i].push_back(new_point);
        }
    }

    return 0;
}



//Concatenate curve
Point LSHCurve::convert_curve_to_vector(Curve curve){

    Point new_point(curve.get_item_id());

    for(int i=0;i<curve.get_curve().size();i++){

        new_point.add_coordinate(curve.get_curve()[i].x);
        new_point.add_coordinate(curve.get_curve()[i].y);
    }

    new_point.set_origin_curve(curve.get_init_curve());

    return new_point;
}











int LSHCurve::padding(){

    int cur_size=0, size = dimension * 2;

    for(int i=0;i<L;i++){

        for(vector<Point>::iterator it = points_set[i].begin(); it != points_set[i].end(); it++){

            cur_size = it->get_dimension();
            while(cur_size < size){

                it->add_coordinate(100000);
                it->add_coordinate(100000);
                cur_size = cur_size + 2;
            }
        }
    }

    return 0;
}




int LSHCurve::testing_dimensions_pointset(){

    int size = dimension * 2;

    for(int i=0;i<L;i++){

        //if(i==0) cout << "\nSize: " << size << endl;

        for(vector<Point>::iterator it = points_set[i].begin(); it != points_set[i].end(); it++){

            if(it->get_dimension() != size){

                assert(1==0);
            }
        }

    }

    return 0;
}







int LSHCurve::convert_query_to_vector(Curve query){

    double x, y, next_x, next_y;
    curve new_curve;
    int cnt_curve=0, index;
    int cur_index=0, max_index=0;

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
        vector<curve> tmp_curve = query.get_curve();
        modified_query_curves[i].set_init_curve(tmp_curve);




        cur_index=0;
        max_index=modified_query_curves[i].get_dimension();
        //Remove consecutive duplicates (x,y)
        while(cur_index<max_index-1){
            
            x = modified_query_curves[i].get_curve()[cur_index].x;
            y = modified_query_curves[i].get_curve()[cur_index].y;
            next_x = modified_query_curves[i].get_curve()[cur_index+1].x;
            next_y = modified_query_curves[i].get_curve()[cur_index+1].y;
            
            if(x == next_x && y == next_y){

                modified_query_curves[i].erase_coordinate(cur_index);
                max_index--;
            }
            else{

                cur_index++;
            }

        }




        //Convert curve to vector
        Point new_point = convert_curve_to_vector(modified_query_curves[i]);
        query_vectors.push_back(new_point);



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
int LSHCurve::clear_query_ds(){

    query_vectors.clear();
    modified_query_curves.clear();
    
    for(int i=0;i<L;i++){

        Curve init_curve("");
        modified_query_curves.push_back(init_curve);
    }
    return 0;
}






inline double euclidean_distance_curve(double p, double q){

    return sqrt((p-q) * (p-q));
}




double LSHCurve::frechet_distance(vector<curve> p, vector<curve> q){

    int m1 = p.size(), m2 = q.size();
    double euclidean_distane, current_index, left, up, diagonal;

    double C[m1][m2];

    for(int i=0;i<m1;i++){

        for(int j=0;j<m2;j++){

            if(i==0 && j==0){
                
                C[i][j] = euclidean_distance_curve(p[i].y, q[j].y);
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

            euclidean_distane = euclidean_distance_curve(p[i].y, q[j].y);

            if(euclidean_distane > current_index){

                C[i][j]=euclidean_distane;
            }
            else{

                C[i][j]=current_index;
            }
        }
    }

    return C[m1-1][m2-1];
}









int LSHCurve::print_query_results(Curve query, string output_file, string algorithm, print_result_info &info){

    ofstream out;
    double current_maf;
    nn_info nn_distance;
    nn_info bruteforce_distance;
    double time, time_true;

    auto time_1 = high_resolution_clock::now();
    find_nearest_neighbor(query, nn_distance, 0);
    auto time_2 = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(time_2 - time_1);
    time = duration.count();
    info.total_time += time;



    auto time_true_1 = high_resolution_clock::now();
    NN_bruteforce(query, bruteforce_distance, 0);
    auto time_true_2 = high_resolution_clock::now();
    auto duration_true = duration_cast<milliseconds>(time_true_2 - time_true_1);
    time_true = duration_true.count();
    info.total_time_true += time_true;


    if(bruteforce_distance.distance != 0){

        current_maf = nn_distance.distance / bruteforce_distance.distance;
        if(current_maf > info.maf)
            info.maf = current_maf;
    }


    out.open(output_file, ios::app);
    if (!out) {
        cout << "File not created!";
    }
    else {

        out << "Query: " << query.get_item_id() << endl;
        out << "Algorithm: " << algorithm << endl;
        out << "Approximate Nearest neighbor: " << nn_distance.item_id << endl;;
        out << "distanceApproximate: " << nn_distance.distance << endl;
        out << "distanceTrue: " << bruteforce_distance.distance << endl;
        out << endl;
        out.close(); 
    }          

    return 0;
}