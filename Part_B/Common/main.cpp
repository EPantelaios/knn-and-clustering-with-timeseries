#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include "../Curve/Clustering_Frechet.hpp"
//#include "../Vector/Clustering.hpp"
//#include "Clustering.hpp"

int main(int argc, char *argv[]){

    int number_of_clusters=0, hash_tables_lsh=3, hash_functions_lsh=4;
    int max_number_hypercube=10; //M hypercube
    int hypercube_dimensions=3; // k=(d') hypercube
    int probes=2, cnt_input=0;
    int complete_param=0, silhouette_param=0; //default value=false
    string update_method, assignment_method;
    string input_file, conf_file, output_file="output.txt";
    curve_struct curve;
    print_result_info info;
    double delta=1;

    check_arguments(argc, argv, input_file, conf_file, output_file, update_method, assignment_method, complete_param, silhouette_param);

    if(file_exist(input_file)==false){

        cout << "Cannot open the input file.\nExit.\n";
        return 1;
    }
    if(file_exist(conf_file)==false){

        cout << "Cannot open the conf file.\nExit.\n";
        return 1;
    }

    if(update_method!="mean_frechet" && update_method!="mean_vector"){

        cout << "Wrong method. Only 'mean_frechet' 'mean_vector' are acceptable for update.\nExit\n";
        return 1;
    }

    if(assignment_method!="Classic" && assignment_method!="LSH" && assignment_method!="Hypercube"){

        cout << "Wrong method. Only 'Classic' 'LSH' and 'Hypercube' are acceptable for assignment.\nExit\n";
        return 1;
    }



    ifstream open_conf(conf_file);
    string line, str;
    int token=0;
    //open and read cluster.conf file
    while(getline(open_conf, line)){

        istringstream ss(line);
        ss >> str;
        ss >> token;

        if(str=="number_of_clusters:"){
            number_of_clusters = token;
        }
        else if(str=="number_of_vector_hash_tables:"){
            hash_tables_lsh = token;
        }
        else if(str=="number_of_vector_hash_functions:"){
            hash_functions_lsh = token;
        }
        else if(str=="max_number_M_hypercube:"){
            max_number_hypercube = token;
        }
        else if(str=="number_of_hypercube_dimensions:"){
            hypercube_dimensions = token;
        }
        else if(str=="number_of_probes:"){
            probes = token;
        }
        else{
            cout << "Cluster conf: Invalid value " << str << "Exit." << endl;
            return 1;
        }
    }

    if(number_of_clusters<1){

        cout << "There must be at least one cluster.\n";
    }

    open_conf.close();

    ifstream open_input(input_file);
    double cur_y, time_x=1.0;
    vector<Point_frechet *> origin_curves;
    Point_frechet *cur_curve;

    //open and read input file
    while(getline(open_input, line)){

        time_x=1;
        istringstream ss(line);
        ss >> str;
        cur_curve = new Point_frechet(str);
        if(update_method=="mean_frechet")
        {
            while(ss >> cur_y){
                
                curve.x = time_x;
                curve.y = cur_y;
                cur_curve->add_curve_coordinate(curve);
                time_x++;
            }
        }
        else
        {
            while(ss >> cur_y)
            {
                cur_curve->add_coordinate(cur_y);
            }

        }
        origin_curves.push_back(cur_curve);

        cnt_input++;

    }
    open_input.close();

    //W parameter (hard-coded). Explained in README
    int w=100;
    if(update_method=="mean_frechet")
    {
        cout << "This is the update_method " << update_method  << endl;
        Clustering_Frechet cluster(origin_curves, number_of_clusters, update_method, assignment_method, hash_tables_lsh, 
                                hash_functions_lsh, w, max_number_hypercube, hypercube_dimensions, probes, delta);
        cluster.initialization();


        cluster.print_results(complete_param, silhouette_param, output_file);
    }
    else
    {
        cout << "This is the update_method " << update_method  << endl;
        Clustering mycluster(origin_curves,number_of_clusters,assignment_method,hash_tables_lsh,hash_functions_lsh,w,max_number_hypercube,hypercube_dimensions,probes);
        mycluster.initialization();
        mycluster.print_results(complete_param,output_file);
    }

    return 0;
}