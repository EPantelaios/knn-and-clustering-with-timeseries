#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include "../Common/utils.hpp"
#include "Clustering.hpp"

Clustering::Clustering(vector<Point_frechet *> dataset, int total_clusters, string method, int L, int k, int w, int max_number_hypercube,
            int hypercube_dimensions, int probes){
    
    centroids.resize(total_clusters);
    this->points_set = dataset;
    this->total_clusters = total_clusters;
    this->method=method;
    this->L=L;
    this->k=k;
    this->w=w;
    this->max_number_hypercube=max_number_hypercube;
    this->hypercube_dimensions=hypercube_dimensions;
    this->probes=probes;
}


//Implement initialization++ (k-Means++)
int Clustering::initialization(){

    int n_clusters=0, rand_index=0, index=0;
    float tmp_distance=0.0, square_result=0.0, max=500.0, normalize=0.0, rand_val=0.0;
    int size_of_points = points_set.size();

    random_device rand_dev;
    mt19937 gen(rand_dev());
    uniform_int_distribution<int> distr(0, size_of_points-1);
    //Choose at random the first centroid from input dataset
    rand_index = distr(gen);
    centroids[n_clusters].centroid = points_set[rand_index];

    //Delete point (first centroid) from input dataset
    points_set.erase(points_set.begin() + rand_index);
    size_of_points = size_of_points - 1;
    n_clusters++;

    //Until all centroids are initialized
    while(n_clusters<total_clusters){

        for(int i=0;i<size_of_points;i++){
            
            points_set[i]->init_sum_square_distance();

            for(int j=n_clusters-1;j<n_clusters;j++){

                tmp_distance = euclidean_distance(points_set[i]->get_coordinates(), centroids[j].centroid->get_coordinates());

                if(tmp_distance < points_set[i]->get_min_distance()){
                    
                    //store min distance for later use
                    points_set[i]->set_min_distance(tmp_distance); 
                    //normalize distance with a constant big integer
                    normalize = points_set[i]->get_min_distance() / max;
                    //D(i)^2 normalize distance
                    square_result = normalize * normalize;
                    //store square distance for later use
                    points_set[i]->set_square_dist(square_result);
                }

                points_set[i]->set_sum_square_dist(points_set[i]->get_square_dist());

                if(i>0){

                    points_set[i]->set_sum_square_dist(points_set[i-1]->get_sum_square_dist());
                }

            }

        }


        uniform_real_distribution<float> distr(0, points_set[size_of_points-1]->get_sum_square_dist());
        rand_val = distr(gen);
        //find lower threshold index with binary search
        index = lower_bound_distance(points_set, size_of_points, rand_val);

        centroids[n_clusters].centroid = points_set[index];
        n_clusters++;
        
        points_set.erase(points_set.begin() + index);
        size_of_points = size_of_points - 1;
    }
    
    return 0;
}


//binary search implementation
int Clustering::lower_bound_distance(vector<Point_frechet *> points_set, int total_size, float value){
    
    int l = 0;
    int h = total_size;

    while(l < h){

        int mid =  l + (h - l) / 2;
        if(value <= points_set[mid]->get_sum_square_dist())
            h = mid;
        else
            l = mid + 1;
    }
    return l;
}


double Clustering::calculate_min_radius(vector<centroid_info> centroids){

    double tmp_distance, min=INT32_MAX;

    for(int i=0;i<total_clusters;i++){

        for(int j=i+1;j<total_clusters;j++){

            tmp_distance = euclidean_distance(centroids[i].centroid->get_coordinates(), centroids[j].centroid->get_coordinates());
            if(tmp_distance<min){

                min = tmp_distance;
            }
        }
    }

    return min;
}




double Clustering::Cluster_Silhouette(vector<double> &s_i, double &stotal)
{
    int i=0, j=0, s=0, pos=0;
    double point_Silhouette=0.0, sum_a=0.0, sum_b=0.0;
    double tmp_distance=0, a_distance=0, b_distance=0;
    double min = INT32_MAX;
    stotal=0.0;
    s_i.resize(total_clusters);

    if(centroids.size()<=1)
    {
        cout << "Invalid amount of clusters" << endl;
        return 1;
    }

    for(i=0;i<centroids.size();i++)
    {//for each cluster
        s_i[i]=0;

        if(centroids[i].points.size()<2)
        {//friendly help
            cout << "Not enough points in the cluster " << endl;
            continue;
        }

        for(j=0;j<centroids[i].points.size();j++)
        {//for each element in the cluster
            min = INT32_MAX;
            sum_a=0;
            for(s=0;s<centroids[i].points.size();s++)
            {//find the average
                if(s==j)
                {
                    continue;
                }
                tmp_distance=euclidean_distance(centroids[i].points[j]->get_coordinates(),centroids[i].points[s]->get_coordinates());
                sum_a=sum_a+tmp_distance;
            }
            a_distance=sum_a/centroids[i].points.size()-1;//a(i)
            //next closest cluster
            for(s=0;s<centroids.size();s++)
            { // find the closest cluster
                if(s==i)
                {
                    continue;
                }
                tmp_distance=euclidean_distance(centroids[i].points[j]->get_coordinates(),centroids[s].centroid->get_coordinates());
                if(min>tmp_distance)
                {
                    min=tmp_distance;
                    pos=s;
                }
            }
            sum_b=0;
            for(s=0;s<centroids[pos].points.size();s++)
            { // average for b
                tmp_distance=euclidean_distance(centroids[i].points[j]->get_coordinates(),centroids[pos].points[s]->get_coordinates());
                sum_b=sum_b+tmp_distance;
            }
            b_distance=sum_b/centroids[pos].points.size();
            //time to find the silhouette
            point_Silhouette=Silhouette(a_distance,b_distance);
            s_i[i]=s_i[i]+point_Silhouette;
        }

        s_i[i]=s_i[i]/centroids[i].points.size();
        cout << i << " Cluster's Silhouette " << s_i[i] << endl;
        stotal=stotal+s_i[i];
    }
    stotal=stotal/centroids.size();
    cout << "Overall Silhouette is " << stotal << endl;

    return 0;
}



double Clustering::Silhouette(double a,double b)
{
    double Silhouette_value;
    if(a>b)
    {
        Silhouette_value =(b-a)/a;
    }
    else
    {
        Silhouette_value=(b-a)/b;
    }
    return Silhouette_value;
}




LSH* Clustering::LSH_reverse_assignment(int L, int k, int w){

    const int max_iter=20;
    const int max_radius_iter=5;
    const float epsilon=1;
    int cnt=0, cur_centroid=0, cnt_new_points=0, sum=0;
    float tmp_distance=0.0, min=INT32_MAX;
    double R, store_R;
    vector<Point_frechet *> means;
    means.resize(total_clusters);
    int size_of_points = points_set.size();

    //store the keys for hashtables that already has visited
    vector<vector<int>> keys;
    vector<int> key;
    for(int i=0;i<total_clusters;i++){

        keys.push_back(key);
    }


    LSH *lsh=new LSH(points_set, L, k, w, points_set[0]->get_dimension());
    lsh->insert_lsh();

    while(cnt<max_iter){
        
        //initialize keys and points at each iteration
        for(int i=0;i<keys.size();i++){

            keys[i].clear();
        }
        clear_points();

        R = calculate_min_radius(centroids) / 2;
        store_R = R;

        for(int i=0;i<total_clusters;i++){
            
            R = store_R;
            for(int j=0;j<max_radius_iter;j++){

                cnt_new_points = lsh->range_search(centroids, i, R, keys);
                //if no new points was added, break
                if(j>1 && cnt_new_points==0){
                    
                    break;
                }
                R = R * 2;
            }
        }


        sum=0;
        for(int i=0;i<total_clusters;i++){

            sum += centroids[i].points.size();
        }


        //Assign each unassigned point with Lloyd approach
        if(sum != points_set.size()){

            for(int i=0;i<size_of_points;i++){

                if(points_set[i]->get_marked() != -1){

                    continue;
                }

                min=INT32_MAX;
                for(int j=0;j<total_clusters;j++){

                    tmp_distance = euclidean_distance(centroids[j].centroid->get_coordinates(), points_set[i]->get_coordinates());

                    if(tmp_distance<min){

                        min = tmp_distance;
                        cur_centroid = j;
                    }
                }

                centroids[cur_centroid].points.push_back(points_set[i]);
            }
        }


        //Store the current centroids
        for(int k=0;k<total_clusters;k++){

            if(cnt==0){
                
                centroids[k].centroid->set_marked(-1);
                centroids[k].centroid->set_cur_distance(INT32_MAX);
                points_set.push_back(centroids[k].centroid);
                size_of_points = points_set.size();

                lsh->insert_K_points(centroids, total_clusters, k);
            }

            means[k] = centroids[k].centroid;
        }
        

        update_means();

        cnt++;

        if(centroids_not_same_position(means, epsilon) == false){

            return lsh;
        }
        //cout << "Here we go again" << endl;

        if(cnt>2)
        {
            for(int i=0;i<total_clusters;i++)
            {
                delete means[i];
            }
        }

        lsh->init_points(keys, total_clusters);
    }

    return lsh;
}




Hypercube* Clustering::Hypercube_reverse_assignment(int k,int M,int probes,int w)
{
    const int max_iter=20;
    const int radius_iter=5;
    const float epsilon=1;
    int cnt=0, cur_centroid=0, cur_point=0, sum=0, cnt_new_points=0;

    Hypercube* xcube=new Hypercube(k,M,probes,0,0,w,points_set[0]->get_dimension(),points_set);
    float tmp_distance=0.0, min=INT32_MAX;
    double R, store_R;

    vector<Point_frechet *> means;
    means.resize(total_clusters);
    vector<vector<int>> keys;
    vector<int> key;
    for(int i=0;i<total_clusters;i++){

        keys.push_back(key);
    }
    int size_of_points = points_set.size();
    xcube->hypercube_insert();
    while(cnt<max_iter){
        
        for(int i=0;i<keys.size();i++){

            keys[i].clear();
        }

        clear_points();

        R = calculate_min_radius(centroids) / 2;
        store_R=R;
        for(int i=0;i<total_clusters;i++){
            
            R=store_R;
            for(int j=0;j<radius_iter;j++){

                cnt_new_points = xcube->range_search(centroids, i, R, keys);
                if(j>1 && cnt_new_points==0){
                    
                    break;
                }
                R = R * 2;
            }
        }

        sum=0;
        for(int i=0;i<total_clusters;i++){

            sum += centroids[i].points.size();
        }
        //Assign each unassigned point with Lloyd approach
        if(sum != points_set.size()){

            for(int i=0;i<size_of_points;i++){

                if(points_set[i]->get_marked() != -1){

                    continue;
                }

                min=INT32_MAX;
                for(int j=0;j<total_clusters;j++){

                    tmp_distance = euclidean_distance(centroids[j].centroid->get_coordinates(), points_set[i]->get_coordinates());

                    if(tmp_distance<min){

                        min = tmp_distance;
                        cur_centroid = j;
                    }
                }

                centroids[cur_centroid].points.push_back(points_set[i]);
            }
        }


        //Store the current centroids
        for(int k=0;k<total_clusters;k++){

            if(cnt==0){
                
                centroids[k].centroid->set_marked(-1);
                centroids[k].centroid->set_cur_distance(INT32_MAX);
                points_set.push_back(centroids[k].centroid);
                size_of_points = points_set.size();

                xcube->hypercube_insert_point(centroids, total_clusters, k);
            }
            means[k] = centroids[k].centroid;
        }
        

        update_means();

        cnt++;

        if(centroids_not_same_position(means, epsilon) == false){

            return xcube;
        }
        cout << "This is the counter " << cnt << endl;
        if(cnt>2)
        {
            for(int i=0;i<total_clusters;i++)
            {
                delete means[i];
            }
        }
        xcube->init_points(keys,total_clusters);
    }
    return xcube;
}



int Clustering::Lloyd_assignment(){

    const int max_iter=20;
    const float epsilon=1;
    int cnt=0, cur_centroid=0;
    float tmp_distance=0.0, min=INT32_MAX;
    vector<Point_frechet *> means;
    means.resize(total_clusters);
    int size_of_points = points_set.size();

    //break only if can overcome max iterations or all the centoirds are almost at the same position
    //with the previous iteration.
    while(cnt<max_iter){

        clear_points();

        for(int i=0;i<size_of_points;i++){

            min=INT32_MAX;
            for(int j=0;j<total_clusters;j++){

                tmp_distance = euclidean_distance(centroids[j].centroid->get_coordinates(), points_set[i]->get_coordinates());

                if(tmp_distance<min){

                    min = tmp_distance;
                    cur_centroid = j;
                }
            }

            centroids[cur_centroid].points.push_back(points_set[i]);
        }


        //Store the current centroids
        for(int k=0;k<total_clusters;k++){

            if(cnt==0){
                points_set.push_back(centroids[k].centroid);
                size_of_points = points_set.size();
            }
            means[k] = centroids[k].centroid;
        }


        update_means();
        //cout << "It does update "<< cnt <<endl ;

        if(centroids_not_same_position(means, epsilon) == false){

            return 0;
        }

        cnt++;
        if(cnt>2)
        {
            for(int i=0;i<total_clusters;i++)
            {
                delete means[i];
            }
        }

    }

    return 0;
}





int Clustering::update_means(){

    int size_of_points=0, i, j;
    int dimensions = points_set[0]->get_dimension();
    Point_frechet *new_point;
    double sum=0.0;

    for(i=0;i<total_clusters;i++){

        new_point = new Point_frechet("Centroid_" + to_string(i+1));
        size_of_points = centroids[i].points.size();

        if(size_of_points==0)
            continue;

        for(j=0;j<dimensions;j++){
            
            sum=0.0;
            for(int k=0;k<size_of_points;k++){

                sum += centroids[i].points[k]->get_a_coordinate(j);
            }

            sum = sum / size_of_points;
            new_point->add_coordinate(sum);
        }

        centroids[i].centroid = new_point;
    }

    return 0;
}


//For each centroid if all coordinates have distance less than epsilon from the previous iteration
//return false and stop the process
bool Clustering::centroids_not_same_position(vector<Point_frechet *> means, float epsilon){

    int dimensions = points_set[0]->get_dimension();
    for(int i=0;i<total_clusters;i++){

        for(int j=0;j<dimensions;j++){

            if(abs(means[i]->get_a_coordinate(j) - centroids[i].centroid->get_a_coordinate(j)) > epsilon)
                return true;
        }
    }

    return false;
}




int Clustering::clear_points(){

    for(int i=0;i<total_clusters;i++){

        centroids[i].points.clear();
    }

    return 0;
}



void Clustering::print_results(int complete_param, string output_file){

    ofstream out;
    double time, stotal;
    vector<double> s_i;
    clock_t start=0, end=0;
    LSH *mylsh;
    Hypercube *xcube;

    if(method=="Classic"){
        
        auto time_1 = high_resolution_clock::now();
        Lloyd_assignment();
        auto time_2 = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(time_2 - time_1);
        time = duration.count();
    }
    else if(method=="LSH"){
        
        auto time_1 = high_resolution_clock::now();
        mylsh=LSH_reverse_assignment(L, k, w);
        auto time_2 = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(time_2 - time_1);
        time = duration.count();
    }
    else if(method=="Hypercube"){

        auto time_1 = high_resolution_clock::now();
        xcube=Hypercube_reverse_assignment(hypercube_dimensions, max_number_hypercube, probes, w);
        auto time_2 = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(time_2 - time_1);
        time = duration.count();
    }


    Cluster_Silhouette(s_i, stotal);

    out.open(output_file);
    if (!out) {
        cout << "File not created!";
    }
    else {

        out << "Algorithm: " << method << endl << endl;

        for(int i=0; i<total_clusters;i++){

            out << "CLUSTER-" << i+1 << " {size: " << centroids[i].points.size() << ", centroid: ";
            centroids[i].centroid->print_coordinates_output(out);
            if(complete_param==1){

                out << endl;
                out << "\nCLUSTER-" << i+1 << " Items_id: {" << endl;
                for(int j=0;j<centroids[i].points.size();j++){
                    
                    if(j==centroids[i].points.size()-1){
                        out << centroids[i].points[j]->get_item_id();
                        break;
                    }
                    out << centroids[i].points[j]->get_item_id() << ", ";
                }
            }
            out << "}" << endl << endl;
        }
        out << "clustering_time: " << time/1000.0 << " seconds" << endl << endl;
        out << "Silhouette: [";
        for(int i=0;i<s_i.size();i++)
            out << s_i[i] << ", ";
        out << "stotal: " << stotal << "]" << endl; 

        out.close();

        if(method=="Classic"){
            
            for(vector<Point_frechet *>::iterator it=points_set.begin(); it != points_set.end(); it++){
        
                delete *it;
            }
        }
        else if(method=="LSH"){
            
            mylsh->~LSH();
        }
        else if(method=="Hypercube"){

           xcube->~Hypercube();
        }
    }
}
