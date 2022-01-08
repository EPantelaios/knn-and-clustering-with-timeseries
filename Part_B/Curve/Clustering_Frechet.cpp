#include <iostream>
#include <vector>
#include <cmath>
#include <random>
//#include "../Common/utils.hpp"
#include "Clustering_Frechet.hpp"

Clustering_Frechet::Clustering_Frechet(vector<Point_frechet *> origin_curves, int total_clusters, string update_method, 
                                        string assignment_method, int L, int k, int max_number_hypercube, 
                                        int hypercube_dimensions, int probes, int w, double delta){
    
    centroids.resize(total_clusters);
    this->origin_curves = origin_curves;
    this->total_clusters = total_clusters;
    this->update_method=update_method;
    this->assignment_method=assignment_method;
    this->max_number_hypercube=max_number_hypercube;
    this->hypercube_dimensions=hypercube_dimensions;
    this->probes=probes;
    this->L=L;
    this->k=k;
    this->w=w;
    this->delta=delta;
    int C_size=origin_curves[0]->get_dimension();
    dimensions=C_size;
    C_size=2*C_size;
    C = new double*[C_size];
    for(int i=0;i<C_size;i++)
    {

        C[i] = new double[C_size];
    }   }


inline double euclidean_distance_curve(curve_struct p, curve_struct q){

    //sqrt( (x1-x2)^2 + (y1-y2)^2 )
    double x = p.x - q.x;
    double y = p.y - q.y;
    return sqrt(x*x + y*y);
}




double Clustering_Frechet::frechet_distance_cluster(vector<curve_struct> *p, vector<curve_struct> *q){

    int m1 = p->size(), m2 = q->size();

    double euclidean_distane, current_index, left, up, diagonal;
    double final_value;


    for(int i=0;i<m1;i++){

        for(int j=0;j<m2;j++){

            if(i==0 && j==0){
                
                C[i][j] = euclidean_distance_curve(p->at(i), q->at(j));
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

            euclidean_distane = euclidean_distance_curve(p->at(i), q->at(j));

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








//Implement initialization++ (k-Means++)
int Clustering_Frechet::initialization(){

    int n_clusters=0, rand_index=0, index=0;
    float tmp_distance=0.0, square_result=0.0, max=500.0, normalize=0.0, rand_val=0.0;
    int size_of_points = origin_curves.size();

    random_device rand_dev;
    mt19937 gen(rand_dev());
    uniform_int_distribution<int> distr(0, size_of_points-1);
    //Choose at random the first centroid from input dataset
    rand_index = distr(gen);
    centroids[n_clusters].centroid = origin_curves[rand_index];

    //Delete point (first centroid) from input dataset
    origin_curves.erase(origin_curves.begin() + rand_index);
    size_of_points = size_of_points - 1;
    n_clusters++;
    vector <curve_struct> *curve1;
    vector <curve_struct> *curve2;

    //Until all centroids are initialized
    while(n_clusters<total_clusters){

        for(int i=0;i<size_of_points;i++){
            
            origin_curves[i]->init_sum_square_distance();

            for(int j=n_clusters-1;j<n_clusters;j++){
                //cout << "This is j "<< j << endl;
                curve1=origin_curves[i]->get_curve_address();
                curve2=centroids[j].centroid->get_curve_address();
                tmp_distance = frechet_distance_cluster(curve1,curve2);

                if(tmp_distance < origin_curves[i]->get_min_distance()){
                    
                    //store min distance for later use
                    origin_curves[i]->set_min_distance(tmp_distance); 
                    //normalize distance with a constant big integer
                    normalize = origin_curves[i]->get_min_distance() / max;
                    //D(i)^2 normalize distance
                    square_result = normalize * normalize;
                    //store square distance for later use
                    origin_curves[i]->set_square_dist(square_result);
                }

                origin_curves[i]->set_sum_square_dist(origin_curves[i]->get_square_dist());

                if(i>0){

                    origin_curves[i]->set_sum_square_dist(origin_curves[i-1]->get_sum_square_dist());
                }

            }

        }

        uniform_real_distribution<float> distr(0, origin_curves[size_of_points-1]->get_sum_square_dist());
        rand_val = distr(gen);
        //find lower threshold index with binary search
        index = lower_bound_distance(origin_curves, size_of_points, rand_val);

        centroids[n_clusters].centroid = origin_curves[index];
        n_clusters++;
        
        origin_curves.erase(origin_curves.begin() + index);
        size_of_points = size_of_points - 1;
    }

    return 0;
}



//binary search implementation
int Clustering_Frechet::lower_bound_distance(vector<Point_frechet *> points_set, int total_size, float value){
    
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


double Clustering_Frechet::calculate_min_radius(vector<centroid_info> centroids){

    double tmp_distance, min=INT32_MAX;
    vector <curve_struct> *curve1;
    vector <curve_struct> *curve2;


    for(int i=0;i<total_clusters;i++){

        for(int j=i+1;j<total_clusters;j++){
            curve1=centroids[i].centroid->get_curve_address();
            curve2=centroids[j].centroid->get_curve_address();

            tmp_distance = frechet_distance_cluster(curve1, curve2);

            if(tmp_distance<min){

                min = tmp_distance;
            }
        }
    }

    return min;
}



double Clustering_Frechet::Cluster_Silhouette(vector<double> &s_i, double &stotal)
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
                tmp_distance=frechet_distance_cluster(centroids[i].points[j]->get_curve_address(), centroids[i].points[s]->get_curve_address());
                sum_a=sum_a+tmp_distance;
            }
            if(centroids[i].points.size()<=1)
                a_distance=sum_a;
            else
                a_distance=sum_a/centroids[i].points.size()-1;//a(i)
            //next closest cluster
            for(s=0;s<centroids.size();s++)
            { // find the closest cluster
                if(s==i)
                {
                    continue;
                }
                tmp_distance=frechet_distance_cluster(centroids[i].points[j]->get_curve_address(), centroids[s].centroid->get_curve_address());
                if(min>tmp_distance)
                {
                    min=tmp_distance;
                    pos=s;
                }
            }
            sum_b=0;
            for(s=0;s<centroids[pos].points.size();s++)
            { // average for b
                tmp_distance=frechet_distance_cluster(centroids[i].points[j]->get_curve_address(), centroids[pos].points[s]->get_curve_address());
                sum_b=sum_b+tmp_distance;
            }
            if(centroids[pos].points.size()<=1)
                b_distance=sum_b;
            else
                b_distance=sum_b/centroids[pos].points.size();
            //time to find the silhouette
            point_Silhouette=Silhouette(a_distance,b_distance);
            s_i[i]=s_i[i]+point_Silhouette;
        }
        if(centroids[i].points.size()<=0)
            s_i[i]=0;
        else
            s_i[i]=s_i[i]/centroids[i].points.size();
        cout << i << " Cluster's Silhouette " << s_i[i] << endl;
        stotal=stotal+s_i[i];
    }
    stotal=stotal/centroids.size();
    cout << "Overall Silhouette is " << stotal << endl;

    return 0;
}



double Clustering_Frechet::Silhouette(double a,double b)
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




Discrete* Clustering_Frechet::LSH_Frechet(int L, int k, int w){

    const int max_iter=10;
    const int max_radius_iter=3; 
    const float epsilon=10;
    int cnt=0, cur_centroid=0, cnt_new_points=0, sum=0;
    float tmp_distance=0.0, min=INT32_MAX;
    double R, store_R;
    vector<Point_frechet *> means;
    means.resize(total_clusters);
    int size_of_points = origin_curves.size();

    //store the keys for hashtables that already has visited
    vector<vector<int>> keys;
    vector<int> key;
    for(int i=0;i<total_clusters;i++){

        keys.push_back(key);
    }

    Discrete *lsh = new Discrete(origin_curves, L, k, w, origin_curves[0]->get_dimension(), delta);

    lsh->hash_all_input();
    int number_of_acceptables_coordinates=lsh->get_dimension();


    while(cnt<max_iter){
        
        //initialize keys and points at each iteration
        for(int i=0;i<keys.size();i++){

            keys[i].clear();
        }

        clear_points();

        R = calculate_min_radius(centroids) / 2;
        store_R = R;

        for(int i=0;i<total_clusters;i++){
            
            lsh->convert_query_to_vector(*centroids[i].centroid);

            R = store_R;
            for(int j=0;j<max_radius_iter;j++){
                
                cnt_new_points = lsh->range_search_discrete(centroids, i, R, keys);
                //if no new points was added, break
                if(j>1 && cnt_new_points==0){
                    
                    break;
                }
                R = R * 2;

            }
            //cout << "This is the current centroid size: " << centroids[i].points.size() << endl;
        }   

        sum=0;
        for(int i=0;i<total_clusters;i++){

            sum += centroids[i].points.size();
        }

        //Assign each unassigned point with Lloyd approach
        if(sum != origin_curves.size()){

            for(int i=0;i<size_of_points;i++){

                if(lsh->get_points_set()[0][i]->get_marked() != -1){

                    continue;
                }

                min=INT32_MAX;
                for(int j=0;j<total_clusters;j++){

                    tmp_distance = frechet_distance_cluster(centroids[j].centroid->get_curve_address(), 
                                                            lsh->get_points_set()[0][i]->get_curve_address());

                    if(tmp_distance<min){

                        min = tmp_distance;
                        cur_centroid = j;
                    }
                }

                centroids[cur_centroid].points.push_back(lsh->get_points_set()[0][i]);
            }
        }   


        //Store the current centroids
        for(int k=0;k<total_clusters;k++){

            if(cnt==0){
                
                centroids[k].centroid->set_marked(-1);
                centroids[k].centroid->set_cur_distance(INT32_MAX);
                lsh->get_points_set_address()->push_back(centroids[k].centroid);
                size_of_points = lsh->get_points_set()[0].size();
                lsh->insert_K_points(centroids, total_clusters, k);
            }

            means[k] = centroids[k].centroid;
        }


        sum=0;
        for(int i=0;i<total_clusters;i++){

            sum += centroids[i].points.size();
        }

        //update_means_frechet();
        update_means();

        cnt++;

        if(centroids_not_same_position(means, epsilon) == false){

            return lsh;
        }
        


        if(cnt>2)
        {
            for(int i=0;i<total_clusters;i++)
            {
                delete means[i];
            }
           // exit(0);

        }

        lsh->init_points(keys, total_clusters);

        cout << "Iteration: " << cnt << endl;
    }


    return lsh;
}






int Clustering_Frechet::Lloyd_assignment(){

    const int max_iter=10;
    const float epsilon=1;
    int cnt=0, cur_centroid=0;
    float tmp_distance=0.0, min=INT32_MAX;
    vector<Point_frechet *> means;
    means.resize(total_clusters);
    int size_of_points = origin_curves.size();

    //break only if can overcome max iterations or all the centoirds are almost at the same position
    //with the previous iteration.
    while(cnt<max_iter){

        clear_points();

        for(int i=0;i<size_of_points;i++){

            min=INT32_MAX;
            for(int j=0;j<total_clusters;j++){

                tmp_distance = frechet_distance_cluster(centroids[j].centroid->get_curve_address(), origin_curves[i]->get_curve_address());

                if(tmp_distance<min){

                    min = tmp_distance;
                    cur_centroid = j;
                }
            }

            centroids[cur_centroid].points.push_back(origin_curves[i]);
        }


        //Store the current centroids
        for(int k=0;k<total_clusters;k++){

            if(cnt==0){
                origin_curves.push_back(centroids[k].centroid);
                size_of_points = origin_curves.size();
            }
            means[k] = centroids[k].centroid;
        }


        update_means();

        //exit(0);
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
        cout << "Iteration " << cnt << endl;
    }

    return 0;
}





int Clustering_Frechet::update_means(){

    int size_of_points=0, i, j;
    int dimensions = origin_curves[0]->get_dimension();
    //cout << "These are the dimensions " << dimensions << endl;
    Point_frechet *new_curve;
    double sum=0.0;
    comp_bin_tree *update_tree;

    for(i=0;i<total_clusters;i++){

        new_curve = new Point_frechet("Centroid_" + to_string(i+1));
        size_of_points = centroids[i].points.size();

        if(size_of_points==0)
        {
            cout << "Unfortunately this cluster doesnt have any points: " << new_curve->get_item_id() << endl;
            continue;
        }

        //here is where the updating happens
        //for starters calculate the tree
        //cout << "Centroid : " << i << endl;
        //cout << "That's how many points currently  in the centroid " << size_of_points<< endl;

        update_tree=new comp_bin_tree(&centroids[i].points);
        new_curve->set_curve(update_tree->comp_bin_tree_compute(update_tree->head,dimensions));
        update_tree->delete_tree(update_tree->head);
        delete update_tree;
        // maybe not necessary
        //filtering(&new_curve->get_curve(),dimensions);

        centroids[i].centroid = new_curve;
        //cout << "First value of the new centroid " << centroids[i].centroid->get_curve()[0].y << endl;

    }

    return 0;
}


//For each centroid if all coordinates have distance less than epsilon from the previous iteration
//return false and stop the process
bool Clustering_Frechet::centroids_not_same_position(vector<Point_frechet *> means, float epsilon){

    int dimensions = origin_curves[0]->get_dimension();
    curve_struct curve_mean, curve_centroid;

    for(int i=0;i<total_clusters;i++){

        for(int j=0;j<dimensions;j++){
            
            curve_mean = means[i]->get_curve_coordinate(j);
            curve_centroid = centroids[i].centroid->get_curve_coordinate(j);

            if(abs(curve_mean.x - curve_centroid.x) > epsilon && abs(curve_mean.y - curve_centroid.y) > epsilon)
                return true;
        }
    }

    return false;
}




int Clustering_Frechet::clear_points(){

    for(int i=0;i<total_clusters;i++){

        centroids[i].points.clear();
    }

    return 0;
}



void Clustering_Frechet::print_results(int complete_param, int silhouette_param, string output_file){

    ofstream out;
    double time, stotal;
    vector<double> s_i;
    clock_t start=0, end=0;
    Discrete *discrete_lsh;

    if(assignment_method=="Classic" && update_method=="mean_frechet"){
        
        auto time_1 = high_resolution_clock::now();
        Lloyd_assignment();
        auto time_2 = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(time_2 - time_1);
        time = duration.count();
    }
    else if(assignment_method=="LSH" && update_method=="mean_frechet"){
        
        auto time_1 = high_resolution_clock::now();
        discrete_lsh=LSH_Frechet(L, k, w);
        auto time_2 = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(time_2 - time_1);
        time = duration.count();
    }

    if(silhouette_param==1){

        Cluster_Silhouette(s_i, stotal);
    }

    out.open(output_file);
    if (!out) {
        cout << "File not created!";
    }
    else {

        out << "Algorithm: " << assignment_method << endl << endl;

        for(int i=0; i<total_clusters;i++){

            out << "CLUSTER-" << i+1 << " {size: " << centroids[i].points.size() << ", centroid: ";
            centroids[i].centroid->print_curves_output(out);
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

        if(silhouette_param==1){

            out << "Silhouette: [";
            for(int i=0;i<s_i.size();i++)
                out << s_i[i] << ", ";
            out << "stotal: " << stotal << "]" << endl; 
        }

        out.close();

        if(update_method=="mean_frechet")
        {
            if(assignment_method=="Classic")
            {
                
                for(vector<Point_frechet *>::iterator it=origin_curves.begin(); it != origin_curves.end(); it++){
            
                    delete *it;
                }
            }
            else if(assignment_method=="LSH")
            {
                
                discrete_lsh->~Discrete();
            }
        }
    }
}








/////////////////

int filtering(vector <curve_struct> *mean_curve,int dimensions){
    
    double epsilon=1;
    int prev_index=0, cur_index=1, next_index=2, total_size;
    vector<curve_struct> *curve;

    total_size=mean_curve->size();
    while (total_size>(dimensions))
    {
        prev_index=0, cur_index=1, next_index=2, 

        epsilon=2*epsilon;
        while(next_index < total_size){
            if(total_size<=dimensions)
            {
                break;
            }
            // will check it out later
            if(abs((*mean_curve)[prev_index].y - (*mean_curve)[cur_index].y) <= epsilon &&
                abs((*mean_curve)[cur_index].y - (*mean_curve)[next_index].y) <= epsilon){

                mean_curve->erase(mean_curve->begin()+cur_index);
                total_size--;
            }
            else{
                prev_index++;
                cur_index++;
                next_index++;
            }
        }
        //total_size=mean_curve->size();

    }
    
    
    return 0;
}

/*inline double euclidean_distance_curve(double p, double q){

    return sqrt((p-q) * (p-q));
}*/

back_track_path *  frechet_distance_and_path(vector <curve_struct> *p,vector <curve_struct> *q ){
//calculates the frechet distance and finds an optimal traversal
    vector<curve_struct> *m1_curve;
    vector<curve_struct> *m2_curve;
    m1_curve=p;
    m2_curve=q;
    int m1 = m1_curve->size(), m2 = m2_curve->size();
    double euclidean_distane, current_index, left, up, diagonal;

    back_track_path *path;
    double **C;
    if(m2==0 || m1==0)
    {
        path=new back_track_path[1];
        path->t=0;
        path->value=0;
        path->path=new int[m1+m2];
        return path;
    }
    C = new double*[m1];
    for(int i=0;i<m1;i++)
    {

        C[i] = new double[m2];
    }   


    for(int i=0;i<m1;i++){

        for(int j=0;j<m2;j++){

            if(i==0 && j==0){
                
                C[i][j] = euclidean_distance_curve(m1_curve->at(i), m2_curve->at(j));
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

            euclidean_distane = euclidean_distance_curve(m1_curve->at(i), m2_curve->at(j));
            if(euclidean_distane > current_index){

                C[i][j]=euclidean_distane;
            }
            else{

                C[i][j]=current_index;
            }
        }
    }
    path = new back_track_path[1];
    path->value=C[m1-1][m2-1];
    path->path=new int[m1+m2]; // allocating enough space
    //time for the backtracking
    int posp=m1-1;
    int posq=m2-1;
    int min=0;
    int counter=0;
    while(posp>0 && posq>0)
    {//starting from the last position and going backwards to the smallest value
        if(C[posp-1][posq-1]<=C[posp][posq-1] && C[posp-1][posq-1]<=C[posp-1][posq])//double step
        {
            path->path[counter]=2;
            posp--;
            posq--;
            //cout << "P is : " << p[posp] << "Q is : " << q[posq] << endl; 
        }
        else if(C[posp-1][posq]<=C[posp-1][posq-1] && C[posp-1][posq]<=C[posp][posq-1]) // P only
        {
            posp--;
            path->path[counter]=1;
            //cout << "P is : " << p[posp] << "Q is : " << q[posq] << endl; 

        }
        else // Q only
        {
            posq--;
            path->path[counter]=0;
            //cout << "P is : " << p[posp] << "Q is : " << q[posq] << endl; 
        }
        counter++;
    }
    while(posp>0)
    {
        posp--;
        //cout << "P is : " << p[posp] << "Q is : " << q[posq] << endl; 
        path->path[counter]=1;
        counter++;
    }
    while(posq>0)
    {
        posq--;
        //cout << "P is : " << p[posp] << "Q is : " << q[posq] << endl; 
        path->path[counter]=0;

        counter++;
    }
    for(int i=0;i<m1;i++){
        delete [] C[i];
    }
    delete [] C;


    path->t=counter;
    return path;
}

vector <curve_struct> find_mean_curve(back_track_path *the_path,vector <curve_struct> p,vector <curve_struct> q)
{
    //we have to consider the two dimensions now
    //traverse it backwards
    vector <curve_struct> mean_curve;
    double mean_value;
    int m1=0;
    int m2=0;
    curve_struct temp;
    if(p.size()==0 && q.size()!=0)
    {
        return q;   
    }
    if(q.size()==0 && p.size()!=0)
    {
        return p;
    }
    if(q.size()==0 && p.size()==0)
    {
        return p;
    }

    temp.x=(p[m1].x+q[m2].x)/2;
    temp.y=(p[m1].y+q[m2].y)/2;
    mean_curve.push_back(temp);
    for(int i=the_path->t-1; i >= 0;i--)
    {//for each state in the traversal , calculate the appropriate value
        if(the_path->path[i]==2) //both
        {
            m1++;
            m2++;
            temp.x=(p[m1].x+q[m2].x)/2;
            temp.y=(p[m1].y+q[m2].y)/2;
            mean_curve.push_back(temp);
        
        }
        else if(the_path->path[i]==1) //p
        {
            m1++;
            temp.x=(p[m1].x+q[m2].x)/2;
            temp.y=(p[m1].y+q[m2].y)/2;
            mean_curve.push_back(temp);

        }
        else//q
        {
            m2++;
            temp.x=(p[m1].x+q[m2].x)/2;
            temp.y=(p[m1].y+q[m2].y)/2;
            mean_curve.push_back(temp);
        }
       // cout << "Let's go" << endl;

    }
    return mean_curve;
}



int *find_path(int tree_counter,int &counter)
{
    int mydiv=tree_counter;
    int *path;
    int num=ceil(log2(tree_counter));
    path = new  int[num+1];
    int next_step;
    while(mydiv>1)
    {
        next_step=mydiv%2;
        mydiv=mydiv/2;
        path[counter]=next_step;
        //cout << "This is the next step " << next_step << " "<< tree_counter << endl;
        counter++;
    }
    return path;
}

comp_bin_node * comp_bin_node_init(Point_frechet *value,int height)
{
    comp_bin_node *temp;
    temp = new comp_bin_node;
    temp->value=value;
    temp->left_child=NULL;
    temp->right_child=NULL;
    temp->cur_height=height;
    return temp;

}

void comp_bin_add_node(int tree_counter,Point_frechet *value,comp_bin_node * head,int height)
{
    comp_bin_node *temp;
    comp_bin_node *temp_parent;
    temp_parent=head;
    temp=comp_bin_node_init(value,height);
    
    int counter=0;
    int *path;
    path = find_path(tree_counter,counter); //find where the node is supposed to be placed
    //cout << counter << " is the counter" << endl;
    for(int i=counter-1;i>0;i--)//from the last element(which is the first one), to the second (which is the parent of the position)
    {
        if(path[i]==1) //uses the path constructed to traverse the tree in logn time
        {
            temp_parent=temp_parent->right_child;
        }
        else
        {
            temp_parent=temp_parent->left_child;
        }
    }
    if(path[0])
    {
        //cout << "RIGHT " << path[0] << endl;
        temp_parent->right_child=temp;
        temp->parent=temp_parent;
    }
    else
    {
        //cout << "LEFT" << path[0] << endl;
        temp_parent->left_child=temp;
        temp->parent=temp_parent;
    }
    
    
};
// really risky, lets hope it works

comp_bin_tree::comp_bin_tree(vector <Point_frechet *> *values)
    {
        tree_counter=0;
        int counter=0;
        height=0;
        int no_nodes;
        int n=values->size();

        height=ceil(log2(n));
        //cout << "this is the height "<< height << endl;
        if(n==1)
        { 
                head = new comp_bin_node;
                head->left_child=NULL;
                head->right_child=NULL;
                head->value=values->at(0);
                head->cur_height=0;
                tree_counter++;
        }
        else
        {
            while(counter<=height)
            {
                if(counter==0)//root
                {
                    head = new comp_bin_node;
                    head->left_child=NULL;
                    head->right_child=NULL;
                    head->value=NULL;
                    head->cur_height=0;
                    tree_counter++;
                }
                else if(counter<height)
                {                    //inner nodes
                    no_nodes=pow(2,counter);
                    for(int i=0;i<no_nodes;i++)
                    {
                        tree_counter++;
                        comp_bin_add_node(tree_counter,NULL,head,counter); //for each level 2^counter nodes
                    }
                }
                else
                {   //children
                    //lets create n nodes
                    for(int i=0;i<n;i++)
                    {
                        tree_counter++;
                        comp_bin_add_node(tree_counter,values->at(i),head,counter);
                    }
                }
                counter++;
            }
        }
};

vector <curve_struct> comp_bin_tree::comp_bin_tree_compute(comp_bin_node *the_node,int dimensions)
{//uses post order tree traversal to recursively calculate the mean approximate mean curve 
    vector <curve_struct> value_left,value_right;
    int flag=0;
    back_track_path *mypath;
    vector <curve_struct> mean_curve;
    if(the_node->cur_height==height)
    {
        return the_node->value->get_curve();
    }
    else
    {
        if(the_node->left_child!=NULL) // in case we have many empty nodes
        {
            value_left=comp_bin_tree_compute(the_node->left_child,dimensions);
            if(the_node->right_child!=NULL)
            {
                value_right=comp_bin_tree_compute(the_node->right_child,dimensions);
            }
            else
            {
                
                flag=1;
            }
            if(!flag)
            {
                if(value_left.size()==0)
                {

                    return value_left;
                }
                mypath=frechet_distance_and_path(&value_left,&value_right); //find path between the two nodes (left,right)
                mean_curve = find_mean_curve(mypath,value_left,value_right);//find the mean curve between them
                filtering(&mean_curve,dimensions);//filter it to have reasonable size
                delete[] mypath->path;
                delete[] mypath;
                //cout << "This is the current size " << mean_curve.size() << endl;
                return mean_curve;//return it
            }
            else
            {
                //cout << "This plays a role value right non existent " << endl;
                return value_left;
            }
        }
        else
        {//tricky
        // this is the error. damn.
        // i have to fix this somehow
            //cout << "This plays a role, value left non existstent " << endl;
            return value_left;
        }
    }
}

void comp_bin_tree::delete_tree(comp_bin_node *the_node)
{
    //now its time to start deleting the tree,like a destructor
    int flag=0;
    if(the_node->cur_height==height)
    {
        delete the_node;
        return;
    }
    else
    {
        if(the_node->left_child!=NULL) // in case we have many empty nodes
        {
            delete_tree(the_node->left_child);
            if(the_node->right_child!=NULL)
            {
                delete_tree(the_node->right_child);
                return;
            }
            else
            {
                flag=1;
            }
            if(!flag)
            {
                delete the_node;
                return;
            }
            else
            {
                //cout << "This plays a role value right non existent " << endl;
                delete the_node;
                return;
            }
        }
        else
        {//tricky
        // this is the error. damn.
        // i have to fix this somehow
            //cout << "This plays a role, value left non existstent " << endl;
            delete the_node;
        }
    }

}

    Clustering_Frechet::~Clustering_Frechet()
    {
        //cout << "This is happening" << endl;

        for(int i=0;i<2*dimensions;i++){

            delete [] C[i];
        }
        delete [] C;


}     
