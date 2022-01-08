#include "Hypercube.hpp"


Hypercube::Hypercube(int k, int M, int probes, int N, int R, int w, int d,vector <Point_frechet *> points_set){
    
    this->k=k;
    this->M=M;
    this->w=w;
    this->N=N;
    this->R=R;
    this->probes=probes;
    this->original_dimensions=d;
    this->points_set=points_set;
    //cout << "M1 = " << M << endl;
    random_device rand_dev;
    mt19937 gen(rand_dev());
    uniform_int_distribution<int> x(0, 1);

    int i=0;
    int j=0;

    for(i=0;i<k;i++){

        H_hash h(w, original_dimensions);
        h_hash.push_back(h);
        map<int, int> mymap;

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
                temp_key = temp_key << k-j; //left bit shift to represent the appropriate bucket
                final_key=final_key|temp_key; //bitwise or
            }
        }
        //cout << "This is the final key for this one " << final_key << endl;
        hyperc.insert(make_pair(final_key,points_set[i])); //insert in the hypercube
        final_key=0;
    }
}


int Hypercube::hypercube_insert_point(vector<centroid_info> centroids, int total_clusters,int index)
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


    for(j=0;j<this->k;j++) //every d'(k) dimension 
    { 
        temp_key=h_hash[j].h_compute(*centroids[index].centroid); //compute h[j]
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
            temp_key = temp_key << k-j; //left bit shift to represent the appropriate bucket
            final_key=final_key|temp_key; //bitwise or
        }
    }
    
    hyperc.insert(make_pair(final_key,centroids[index].centroid)); //insert in the hypercube
    return 0;
}


int Hypercube::init_points(vector<vector<int>> keys,int total_clusters){

    for(int i=0;i<total_clusters;i++)
    {
        for(int j=0;j<keys[i].size();j++){

            auto its = hyperc.equal_range(keys[i][j]);
            for(auto it = its.first; it != its.second; ++it){

                it->second->set_marked(-1);
                it->second->set_cur_distance(INT32_MAX);
            }
        }
    }
    return 0;
}

int Hypercube::range_search(vector<centroid_info> &centroids,int centroid_num,double R,vector<vector<int>> &keys)
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
        temp_key=h_hash[j].h_compute(*centroids[i].centroid); //compute centroid bucket
        it = f_function[j].find(temp_key);  //search in f[j]
        if(it==f_function[j].end()) //does not exist
        {
            f_function[j].insert(pair<int,int>(temp_key,x(gen))); //rand 0,1 
        }
        it=f_function[j].find(temp_key); //finds it definitely
        if(it->second==1)
        {
            temp_key=1;
            temp_key = temp_key << k-j; //left bit shift to represent the appropriate bucket
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
    knn_info_hyper mypoint;//change name asap
    double tmp_distance;
    int flag=0;
    string point_name;
    int pos=0, cnt_new_points=0;
    auto its = hyperc.equal_range(final_key);

    while(counter<M && prob_counter<probes)//carefull of empty buckets
    {
        if(flag)
        {
            if(cur_nearby_buckets==nearby_buckets.size())
            {
                cout << "Changing distance " << endl;
                cur_nearby_buckets=0;
                nearby_buckets.clear();
                nearby_buckets_str.clear();
                nearby_buckets=hamming(final_key,k,nearby_buckets_str,cur_hamming);
                cur_hamming++;
            }
            final_key=nearby_buckets[cur_nearby_buckets];
            prob_counter++;
            cur_nearby_buckets++;
            its=hyperc.equal_range(final_key);
            keys[centroid_num].push_back(final_key);

        }
        if(flag==0)
        {
            flag=1;
            keys[centroid_num].push_back(final_key);
        }
        for (auto it = its.first; it != its.second; ++it){

            counter++;
            //Skip this iteration. Point has already been assigned to some cluster
            if(it->second->get_marked() == centroid_num){
                
                if(counter>=M)
                {
                    break;
                }
                continue;
            }        
                
            tmp_distance = euclidean_distance(centroids[centroid_num].centroid->get_coordinates(), it->second->get_coordinates());

            if(tmp_distance < R){
    
                if(it->second->get_marked() != -1){
                    
                    if(tmp_distance < it->second->get_cur_distance()){

                        it->second->set_cur_distance(tmp_distance);
                        int tmp_index = it->second->get_marked();
                        it->second->set_marked(centroid_num);
                        centroids[centroid_num].points.push_back(it->second);
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
                    it->second->set_marked(centroid_num);
                    centroids[centroid_num].points.push_back(it->second);
                    cnt_new_points++;
                }

            }

            if(counter>=M)
            {
                break;
            }
        } 

    }

    return cnt_new_points;
}

int Hypercube::range(Point_frechet query,int R,vector <knn_info_hyper> &point_info)
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
            temp_key = temp_key << k-j; //left bit shift to represent the appropriate bucket
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
            {
                cout << "Changing distance " << endl;
                cur_nearby_buckets=0;
                nearby_buckets.clear();
                nearby_buckets_str.clear();
                nearby_buckets=hamming(final_key,k,nearby_buckets_str,cur_hamming);
                cur_hamming++;
            }
            final_key=nearby_buckets[cur_nearby_buckets];
            prob_counter++;
            cur_nearby_buckets++;
            its=hyperc.equal_range(final_key);
        }
        if(flag==0)
        {
            flag=1;
        }
        for (auto it = its.first; it != its.second; ++it){
                counter++;
                    
                tmp_distance = euclidean_distance(query.get_coordinates(), it->second->get_coordinates());
                //cout << "___Distance_ID___: " << tmp_distance << endl;
                if(tmp_distance<R)
                {
                    mypoint.distance=tmp_distance;
                    mypoint.item_id=it->second->get_item_id();
                    point_info.push_back(mypoint);
                }
        }             

    }
    return 0;
}

int Hypercube::knn(Point_frechet query,int N,vector <knn_info_hyper> &point_info)
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
            temp_key = temp_key << k-j; //left bit shift to represent the appropriate bucket
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
    int tmp_distance;
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
            {
                cout << "Changing distance " << endl;
                cur_nearby_buckets=0;
                nearby_buckets.clear();
                nearby_buckets_str.clear();
                nearby_buckets=hamming(final_key,k,nearby_buckets_str,cur_hamming);
                cur_hamming++;
            }
            final_key=nearby_buckets[cur_nearby_buckets];
            prob_counter++;
            cur_nearby_buckets++;
            its=hyperc.equal_range(final_key);
        }
        if(flag==0)
        {
            flag=1;
        }
        for (auto it = its.first; it != its.second; ++it){
                counter++;
                    
                tmp_distance = euclidean_distance(query.get_coordinates(), it->second->get_coordinates());
                //cout << "___Distance_ID___: " << tmp_distance << endl;
                if(counter<N)
                {
                    mypoint.distance=tmp_distance;
                    mypoint.item_id=it->second->get_item_id();
                    for(pos=0;pos<point_info.size();pos++)
                    {
                        if(point_info[pos].distance>tmp_distance)
                        {
                            break;
                        }
                    }
                    insert_it=point_info.begin();
                    insert_it=point_info.insert(insert_it+pos,mypoint);
                    
                }
                else if(point_info[point_info.size()-1].distance > tmp_distance)
                {

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
                    if(point_info.size()>N)
                    {
                    insert_it=point_info.erase(point_info.end());
                    }
                }
            
        }

    }
    return 0;
}


int Hypercube::find_bit(int num,int pos,int k)
{
    int mask =  1 << pos;
    int masked_n = num & mask;
    int thebit = masked_n >> pos;
    return thebit;

}

int Hypercube::reverse_to_int(string mystring)
{
    int length=mystring.length();
    int i=0;
    int value=0;
    string one="1";
    for(i=length-1;i>=0;i--)
    {
        if( mystring[i]==one[0])
        {
            //cout << "Yes" << endl ;
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

    for(vector<Point_frechet *>::iterator it=points_set.begin(); it != points_set.end(); it++){
    
        delete *it;
    }

}
