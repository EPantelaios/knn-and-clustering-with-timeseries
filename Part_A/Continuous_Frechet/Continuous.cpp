#include "Continuous.hpp"

Continuous::Continuous(vector<Curve_cont *> dataset, int k, int w, int d, double delta){

    this->origin_curves = dataset;
    this->k=k;
    this->w=w;
    this->origin_dimension=d;
    this->delta=delta;
    this->table_size=origin_curves.size()/32;

    if(table_size==0)
        table_size=1;

    G_hash g(k, w, origin_dimension, this->table_size);
    g_hash = g;

    random_device rand_dev;
    mt19937 gen(rand_dev());
    uniform_real_distribution<float> x(0, delta);
    
    t = x(gen);
}



Continuous::~Continuous(){

    for(vector<Curve_cont *>::iterator it=origin_curves.begin(); it != origin_curves.end(); it++){
    
        delete *it;
    }
}



//insert point to LSH hashtable
int Continuous::insert_lsh(){

    int key, size;
    unsigned long long int tmp;

    size = points_set.size();

    for(int i=0;i<size;i++){

        tmp = g_hash.ID_compute(points_set[i]);
        points_set[i].set_ID(tmp);

        key = tmp % this->table_size;
        hashtable.insert(make_pair(key, points_set[i]));
    }
    
    return 0;
}



int Continuous::find_nearest_neighbor(Curve query, nn_info &min, unsigned long long int ID){

    Frechet::Continuous::Distance result;
    min.distance = INT32_MAX;
    int key=0, size=0;

    key = ID % table_size;
    auto its = hashtable.equal_range(key);

    for(auto it = its.first; it != its.second; ++it){

        Curve curve_input(1);
        size = it->second.get_filtered_curve().size();
        for(int i=0;i<size;i++){

            Point p(1);
            p.set(0, it->second.get_filtered_curve()[i]);
            curve_input.push_back(p);
        }

        result = Frechet::Continuous::distance(curve_input, query);

        if(min.distance>result.value){

            min.distance = result.value;
            min.item_id=it->second.get_item_id();
        }
    }

    //In case bucket is empty, we have to search in other buckets
    if(min.distance==INT32_MAX){

        for(int i=0;i<table_size;i++){

            if(i==key)
                continue;

            auto its = hashtable.equal_range(i);
            
            for(auto it = its.first; it != its.second; ++it){

                Curve curve_input(1);
                size = it->second.get_filtered_curve().size();
                for(int i=0;i<size;i++){

                    Point p(1);
                    p.set(0, it->second.get_filtered_curve()[i]);
                    curve_input.push_back(p);
                }

                result = Frechet::Continuous::distance(curve_input, query);

                if(min.distance>result.value){

                    min.distance = result.value;
                    min.item_id=it->second.get_item_id();
                }
            }
        }
    }

    return 0;
}



int Continuous::NN_bruteforce(Curve query, nn_info &min){
    
    int size;
    min.distance = INT32_MAX;
    Frechet::Continuous::Distance result;

    for(int i=0;i<origin_curves.size();i++){

        Curve curve_input(1);
        size = origin_curves[i]->get_filtered_curve().size();

        for(int j=0;j<size;j++){

            Point p(1);
            p.set(0, origin_curves[i]->get_filtered_curve()[j]);
            curve_input.push_back(p);
        }
        result = Frechet::Continuous::distance(curve_input, query);
        if(min.distance>result.value){

            min.distance=result.value;
            min.item_id = origin_curves[i]->get_item_id();
        }
    }

    return 0;
}






int Continuous::hash_all_input(){

    assert(!filtering());
    assert(!snap_curves());
    assert(!minima_maxima());
    assert(!convert_input_curves_to_vectors());
    assert(!padding());
    assert(!testing_dimensions_pointset());
    assert(!insert_lsh());

    return 0;
}



int Continuous::filtering(){

    int prev_index=0, cur_index=1, next_index=2, total_size;
    vector<double> *curve;

    for(int i=0;i<origin_curves.size();i++){
        
        prev_index=0, cur_index=1, next_index=2;
        total_size=origin_curves[i]->get_dimension();
        curve = origin_curves[i]->get_curve_address();

        while(next_index < total_size){

            if(abs((*curve)[prev_index] - (*curve)[cur_index]) <= EPSILON &&
               abs((*curve)[cur_index] - (*curve)[next_index]) <= EPSILON){

                origin_curves[i]->erase_coordinate(cur_index);
                total_size--;
            }
            else{
                prev_index++;
                cur_index++;
                next_index++;
            }
        }

        origin_curves[i]->set_filter_curve(origin_curves[i]->get_curve());
    }
    
    return 0;
}





int Continuous::snap_curves(){

    int cnt_curve=0;
    double x, y;
    double new_curve;

    //Now the origin_curves have been filtered
    for(vector<Curve_cont *>::iterator it = origin_curves.begin(); it != origin_curves.end(); it++){

        Curve_cont Y((*it)->get_item_id());
        modified_curves.push_back(Y);

        for(int j=0;j<(*it)->get_dimension();j++){

            y = (*it)->get_curve()[j];
            new_curve = floor((y + t) / delta) * delta;  
            modified_curves[cnt_curve].add_curve_coordinate(new_curve);
        }
        
        //store the filtered curve before proceeding with the conversions
        vector<double> tmp_curve = (*it)->get_curve();
        modified_curves[cnt_curve].set_filter_curve(tmp_curve);
        cnt_curve++;
    }

    return 0;
}





//Remove coordinate between min and max per three coordinates
int Continuous::minima_maxima(){

    int prev_index=0, cur_index=1, next_index=2, total_size;
    vector<double> *curve;

    for(int i=0;i<modified_curves.size();i++){

        prev_index=0, cur_index=1, next_index=2;
        total_size=modified_curves[i].get_dimension();
        curve = modified_curves[i].get_curve_address();

        while(next_index < total_size){

            if(((*curve)[prev_index] <= (*curve)[cur_index] && (*curve)[cur_index] <= (*curve)[next_index]) || 
                ((*curve)[prev_index] >= (*curve)[cur_index] && (*curve)[cur_index] >= (*curve)[next_index])){

                modified_curves[i].erase_coordinate(cur_index);
                total_size--;
            }
            else{
                prev_index++;
                cur_index++;
                next_index++;
            }
        }
    }

    return 0;
}





//Concatenate all input curves
int Continuous::convert_input_curves_to_vectors(){

    for(vector<Curve_cont>::iterator it = modified_curves.begin(); it != modified_curves.end(); it++){

        Point_cont new_point = convert_curve_to_vector(*it);
        points_set.push_back(new_point);
    }

    return 0;
}




//Concatenate curve
Point_cont Continuous::convert_curve_to_vector(Curve_cont curve){

    Point_cont new_point(curve.get_item_id());

    for(int i=0;i<curve.get_curve().size();i++){

        new_point.add_coordinate(curve.get_curve()[i]);
    }

    new_point.set_filtered_curve(curve.get_filtered_curve());

    return new_point;
}





int Continuous::padding(){

    int cur_size=0, size = origin_dimension;

    for(vector<Point_cont>::iterator it = points_set.begin(); it != points_set.end(); it++){

        cur_size = it->get_dimension();
        while(cur_size < size){

            it->add_coordinate(100000);
            cur_size = cur_size + 1;
        }
    }

    return 0;
}




int Continuous::testing_dimensions_pointset(){

    int size = origin_dimension;

    for(vector<Point_cont>::iterator it = points_set.begin(); it != points_set.end(); it++){

        if(it->get_dimension() != size){

            assert(1==0);
        }
    }

    return 0;
}





int Continuous::convert_query_to_vector(Curve_cont &query){

    double new_curve, y;
    int prev_index=0, cur_index=1, next_index=2, total_size, dimension;
    vector<double> *curve;

    prev_index=0, cur_index=1, next_index=2;
    total_size=query.get_dimension();
    curve = query.get_curve_address();

    origin_query_curve = query;

    while(next_index < total_size){

        if(abs((*curve)[prev_index] - (*curve)[cur_index]) <= EPSILON &&
            abs((*curve)[cur_index] - (*curve)[next_index]) <= EPSILON){

            query.erase_coordinate(cur_index);
            total_size--;
        }
        else{
            prev_index++;
            cur_index++;
            next_index++;
        }
    }

    query.set_filter_curve(query.get_curve());

    //snapping
    dimension=query.get_dimension();
    for(int j=0;j<dimension;j++){

        y = query.get_curve()[j];
        new_curve = floor((y + t) / delta) * delta;  
        query.set_curve_coordinate(j, new_curve);
    }   


    //minima-maxima
    prev_index=0, cur_index=1, next_index=2;
    total_size=query.get_dimension();
    curve = query.get_curve_address();

    while(next_index < total_size){

        if(((*curve)[prev_index] <= (*curve)[cur_index] && (*curve)[cur_index] <= (*curve)[next_index]) || 
            ((*curve)[prev_index] >= (*curve)[cur_index] && (*curve)[cur_index] >= (*curve)[next_index])){

            query.erase_coordinate(cur_index);
            total_size--;
        }
        else{
            prev_index++;
            cur_index++;
            next_index++;
        }
    }


    //Convert curve to vector
    Point_cont new_point = convert_curve_to_vector(query);
    query_vector = new_point;


    int cur_size=0, size=origin_dimension;

    cur_size = query_vector.get_dimension();
    while(cur_size < size){

        query_vector.add_coordinate(100000);
        cur_size = cur_size + 1;
    }

    return 0;
}




int Continuous::print_query_results(Curve_cont query, string output_file, string algorithm, print_result_info &info){

    ofstream out;
    nn_info nn_distance;
    nn_info bruteforce_distance;
    double time, time_true, current_maf;
    int size;
    unsigned long long int ID;

    convert_query_to_vector(query);

    Curve curve_query(1);
    size = query.get_filtered_curve().size();
    for(int i=0;i<size;i++){

        Point p(1);
        p.set(0, query.get_filtered_curve()[i]);
        curve_query.push_back(p);
    }

    ID = g_hash.ID_compute(query_vector);
    auto time_1 = high_resolution_clock::now();
    find_nearest_neighbor(curve_query, nn_distance, ID);
    auto time_2 = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(time_2 - time_1);
    time = duration.count();
    info.total_time += time;


    auto time_true_1 = high_resolution_clock::now();
    NN_bruteforce(curve_query, bruteforce_distance);
    auto time_true_2 = high_resolution_clock::now();
    auto duration_true = duration_cast<microseconds>(time_true_2 - time_true_1);
    time_true = duration_true.count();
    info.total_time_true += time_true;

    if(bruteforce_distance.distance != 0){

        current_maf = nn_distance.distance / bruteforce_distance.distance;
        if(current_maf > info.maf)
            info.maf = current_maf;
        
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
}