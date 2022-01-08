#include "Hypercube.hpp"

Hypercube::Hypercube(int k, int M, int probes, int N, int R, int w, int d, vector <Point_hc *> points_set){

    this->k=k;
    this->M=M;
    this->w=w;
    this->N=N;
    this->R=R;
    this->probes=probes;
    this->original_dimensions=d;
    this->points_set=points_set;

    for(int i=0;i<k;i++){

        H_hash h(w, original_dimensions);
        h_hash.push_back(h);
        map<int, int> mymap;
        // k f_functions , each for one h 
        f_function.push_back(mymap);
    }
}

void Hypercube::hypercube_insert()
{
    int i=0;
    int j=0;
    int counter=0;
    map<int,int>::iterator it;
    unsigned int final_key=0; //only works with k <= 32
    int long long temp_key=0;
    int size=points_set.size();
    random_device rand_dev;
    mt19937 gen(rand_dev());
    uniform_int_distribution<int> x(0, 1); //for our f function


    for(i=0;i<size;i++) //every element
    {
        for(j=0;j<this->k;j++) //every d'(k) dimension 
        { 
            temp_key=h_hash[j].h_compute(*points_set[i]); //compute h[j]
            it = f_function[j].find(temp_key);  //search in f[j]
            if(it==f_function[j].end()) //does not exist
            {
                f_function[j].insert(pair<int,int>(temp_key,x(gen))); //rand 0,1 

                //cout << "Inserted " << temp_key << " at " << j <<endl; 
            }
            it=f_function[j].find(temp_key); //finds it definitely
            if(it->second==1)
            {
                temp_key=1;
                temp_key = temp_key << (k-j); //left bit shift to represent the appropriate bucket
                final_key=final_key|temp_key; //bitwise or
            }
        }
        hyperc.insert(make_pair(final_key,points_set[i])); //insert in the hypercube
        final_key=0;
    }
}



int Hypercube::range(Point_hc query,int R,vector<string> &range_search_names)
{
    map<int,int>::iterator it;
    int i=0;
    int j=0;
    int final_key=0;
    int temp_key=0;
    random_device rand_dev;
    mt19937 gen(rand_dev());
    uniform_int_distribution<int> x(0, 1); //for our f function

    for(j=0;j<this->k;j++) //every d'(k) dimension 
    { 
        temp_key=h_hash[j].h_compute(query); //compute h[j]
        it = f_function[j].find(temp_key);  //search in f[j]
        if(it==f_function[j].end()) //does not exist
        {
            f_function[j].insert(pair<int,int>(temp_key,x(gen))); //rand 0,1 

                //cout << "Inserted " << temp_key << " at " << j <<endl; 
        }
        it=f_function[j].find(temp_key); //finds it definitely
        if(it->second==1)
        {
            temp_key=1;
            temp_key = temp_key << (k-j); //left bit shift to represent the appropriate bucket
            final_key=final_key|temp_key; //bitwise or
        }
    }
    int counter=0; //for M
    int prob_counter=0; // for
    int cur_nearby_buckets=0;// for the current nearby buckets
    vector <int> nearby_buckets;
    vector <string> nearby_buckets_str;
    int cur_hamming=1;
    nearby_buckets=hamming(final_key,k,nearby_buckets_str,cur_hamming);
    cur_hamming++;
    int min=999999;
    double tmp_distance;
    int flag=0;
    int pos=0;
    auto its = hyperc.equal_range(final_key);

    while(counter<M && prob_counter<probes)//carefull of empty buckets
    {
        if(flag)
        {
            if(cur_nearby_buckets==nearby_buckets.size())
            {//time to increase the hamming distance
                //cout << "Changing distance " << endl;
                cur_nearby_buckets=0;
                nearby_buckets.clear();
                nearby_buckets_str.clear();
                nearby_buckets=hamming(final_key,k,nearby_buckets_str,cur_hamming);
                cur_hamming++;
            }
            //changes to a nearby bucket
            final_key=nearby_buckets[cur_nearby_buckets];
            prob_counter++;
            cur_nearby_buckets++;
            its=hyperc.equal_range(final_key);
        }
        if(flag==0)
        {
            //first iteration
            flag=1;
        }
        for(auto it = its.first; it != its.second; ++it){

            counter++;
            tmp_distance = euclidean_distance(query.get_coordinates(), it->second->get_coordinates());
            if(tmp_distance<R)
            {
                range_search_names.push_back(it->second->get_item_id());
            }
            if(counter>=M)
            {             
                break;
            }
        }

    }
    return 0;
}

int Hypercube::knn(Point_hc query, int N, vector <knn_info_hyper> &point_info)
{
    map<int,int>::iterator it;
    int i=0;
    int j=0;
    int final_key=0;
    int temp_key=0;
    random_device rand_dev;
    mt19937 gen(rand_dev());
    uniform_int_distribution<int> x(0, 1); //for our f function

    for(j=0;j<this->k;j++) //every d'(k) dimension 
    { 
        temp_key=h_hash[j].h_compute(query); //compute h[j]
        it = f_function[j].find(temp_key);  //search in f[j]
        if(it==f_function[j].end()) //does not exist
        {
            f_function[j].insert(pair<int,int>(temp_key,x(gen))); //rand 0,1 

                //cout << "Inserted " << temp_key << " at " << j <<endl; 
        }
        it=f_function[j].find(temp_key); //finds it definitely
        if(it->second==1)
        {
            temp_key=1;
            temp_key = temp_key << (k-j); //left bit shift to represent the appropriate bucket
            final_key=final_key|temp_key; //bitwise or
        }
    }
    //something must be added for hamming distance
    int counter=0; //for M
    int prob_counter=0; // for
    int cur_nearby_buckets=0;// for the current nearby buckets
    vector <int> nearby_buckets;
    vector <string> nearby_buckets_str;
    int cur_hamming=1;
    nearby_buckets=hamming(final_key,k,nearby_buckets_str,cur_hamming);
    cur_hamming++;
    int min=999999;
    knn_info_hyper mypoint;//change name asap
    double tmp_distance;
    int flag=0;
    string point_name;
    std::vector<knn_info_hyper>::iterator insert_it;
    int pos=0;
    auto its = hyperc.equal_range(final_key);

    while(counter<M && prob_counter<probes)//carefull of empty buckets
    {
        if(flag)
        {
            if(cur_nearby_buckets==nearby_buckets.size())
            {//increasing hamming distance
                //cout << "Changing distance " << endl;
                cur_nearby_buckets=0;
                nearby_buckets.clear();
                nearby_buckets_str.clear();
                nearby_buckets=hamming(final_key,k,nearby_buckets_str,cur_hamming);
                cur_hamming++;
            }
            //changing to a nearby bucket
            final_key=nearby_buckets[cur_nearby_buckets];
            prob_counter++;
            cur_nearby_buckets++;
            its=hyperc.equal_range(final_key);
        }
        if(flag==0)
        {
            //first iteration
            flag=1;
        }
        for (auto it = its.first; it != its.second; ++it){
                counter++;
                    
                tmp_distance = euclidean_distance(query.get_coordinates(), it->second->get_coordinates());
                if(counter<=N)
                { //for the first stage where we have not even found N neighbours
                    mypoint.distance=tmp_distance;
                    mypoint.item_id=it->second->get_item_id();
                    if(counter!=1)
                    {
                        for(pos=0;pos<point_info.size();pos++)
                        {   //find the correct place for insertion, so that the vector remains sorted
                            if(point_info[pos].distance>tmp_distance)
                            {
                                break;
                            }
                        }
                            insert_it=point_info.begin();
                            insert_it=point_info.insert(insert_it+pos,mypoint);
                    }
                    else
                    {
                        mypoint.distance=tmp_distance;
                        mypoint.item_id=it->second->get_item_id();
                        point_info.push_back(mypoint);
                            
                    }
                    
                }
                else if(point_info[point_info.size()-1].distance > tmp_distance)
                {
                    //find the appropriate pos to insert it

                    for(pos=0;pos<point_info.size();pos++)
                    {
                        if(point_info[pos].distance>tmp_distance)
                        {
                            break;
                        }
                    }
                    mypoint.distance=tmp_distance;
                    mypoint.item_id=it->second->get_item_id();
                    insert_it=point_info.begin();
                    insert_it=point_info.insert(insert_it+pos,mypoint);
                    if(point_info.size()>N) //delete the last one 
                    {
                    insert_it=point_info.erase(point_info.end());
                    }
                }
            
            if(counter>=M) // break if it reaches the threshold
            {
                break;
            }
        }

    }
    return 0;
}



int Hypercube::find_bit(int num,int pos,int k)
{
    //find the corresponding bit
    int mask =  1 << pos;
    int masked_n = num & mask;
    int thebit = masked_n >> pos;
    return thebit;

}

int Hypercube::reverse_to_int(string mystring)
{
    //turns a string into an int
    int length=mystring.length();
    int i=0;
    int value=0;
    string one="1";
    for(i=length-1;i>=0;i--)
    {
        if( mystring[i]==one[0])
        {
            value=value+1;
        }
        if(i==0)
        {
            break;
        }
        value=value << 1;
    }
    return value;
}
void Hypercube::hamming_distance_strings(string mystring,int cur_height,int hamming_distance,vector<string> &all_comb,int k,int start_point)
{
    //creates the neighbouring buckets to our hypercube based on the given hamming distance as strings
    int i=0;
    string to_add;
    int endpoint=k-hamming_distance+(hamming_distance-cur_height)+1;
    for(i=start_point;i<endpoint;i++)
    {
        to_add=mystring;
        if(to_add[i]=='1')
        {
            to_add[i]='0';
        }
        else
        {
            to_add[i]='1';
        }
        if(cur_height==1)
        {
            all_comb.push_back(to_add);
        }
        else
        {
            hamming_distance_strings(to_add,cur_height-1,hamming_distance,all_comb,k,i+1);
        }
    }

}

vector<int> Hypercube::hamming(int l,int k,vector <string> &everything,int distance)
{
    //creates and returns the neighbouring buckets to our hypercube based on the given hamming distance as ints
    vector <int> hamming_vector;
    int original = l;
    int temp;
    int current;
    string hamming_str;
    int i=0;
    for(i=0;i<k;i++)
    {
        temp=original;
        current=find_bit(original,i, k);
        if(current==0)
        {
            hamming_str.append("0");
        }
        else
        {
            hamming_str.append("1");
        }

    }
    //cout << "This is the hamming_str " << hamming_str <<endl;
    int value;
    value=reverse_to_int(hamming_str);
    //cout << "This is the value " << value << endl;
    hamming_distance_strings(hamming_str,distance,distance,everything,k,0);
    for(i=0;i<everything.size();i++)
    {
        hamming_vector.push_back(reverse_to_int(everything[i]));
    }

    return hamming_vector;
}

Hypercube::~Hypercube(){
    for(vector<Point_hc *>::iterator it=points_set.begin(); it != points_set.end(); it++){
    
        delete *it;
    }

}



int Hypercube::NN_bruteforce(Point_hc query, knn_info_hyper &info){
    
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


void Hypercube::swap(knn_info_hyper &x, knn_info_hyper &y){ 

    knn_info_hyper temp = x; 
    x = y; 
    y = temp; 
}

void Hypercube::bubbleSort(vector<knn_info_hyper> &items, int n){ 

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


int Hypercube::kNN_bruteforce(Point_hc query, int N, vector<knn_info_hyper> &min_distances){

    knn_info_hyper tmp;

    for(int i=0;i<points_set.size();i++){
        
        tmp.distance = euclidean_distance(points_set[i]->get_coordinates(), query.get_coordinates());
        tmp.item_id = points_set[i]->get_item_id();
        min_distances.push_back(tmp);
    }

    bubbleSort(min_distances, min_distances.size());

    return 0;
}



void Hypercube::print_query_results(Point_hc query, string output_file, string algorithm, print_result_info &info){
    
    ofstream out;
    double current_maf;
    vector<knn_info_hyper> nn_distance;
    knn_info_hyper bruteforce_distance;
    double time, time_true;

    auto time_1 = high_resolution_clock::now();
    knn(query, 1, nn_distance);
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

        current_maf = nn_distance[0].distance / bruteforce_distance.distance;
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
        out << "Approximate Nearest neighbor: " << nn_distance[0].item_id << endl;
        out << "True Nearest neighbor: " << bruteforce_distance.item_id << endl;
        out << "distanceApproximate: " << nn_distance[0].distance << endl;
        out << "distanceTrue: " << bruteforce_distance.distance << endl;
        out << endl;
        out.close(); 
    }
}

