#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include "LSHCurve.hpp"


int main(int argc, char *argv[]){

    int k=4, L=5, dd=14, M=10, probes=2, time_x=1;
    double delta=-1.0;
    string input_file, query_file, output_file="output.txt", algorithm, metric;
    curve curve_struct;
    print_result_info info;
    info.maf=INT32_MIN;

    random_device rand_dev;
    mt19937 gen(rand_dev());
    uniform_int_distribution<int> distr(2, 6);
    int w = distr(gen);

    check_arguments(argc, argv, input_file, query_file, output_file, k, L, M, probes, algorithm, metric, delta);


    cout << input_file << "    "<< query_file << "    "<< k << "    "<< L << "    "<< M << "    "<< probes << "    "<< algorithm <<\
        "    "<<metric << "    "<< delta << endl << endl;


    if(file_exist(input_file)==false){

        cout << "Cannot open the input file.\nExit.\n";
        return 1;
    }
    if(file_exist(query_file)==false){

        cout << "Cannot open the query file.\nExit.\n";
        return 1;
    }

    if(algorithm!="LSH" && algorithm!="Hypercube" && algorithm!="Frechet"){
        
        cout << "Wrong algorithm. Only 'LSH', 'Hypercube' and 'Frechet' are acceptable.\nExit\n";
        return 1;
    }

    if(algorithm=="Frechet"){

        if(delta==-1.0){

            cout <<"Delta must be included.\nExit\n";
            return 1;
        }

        if(metric!="discrete" && metric!="continuous"){

            cout << "'metric' parameter: Only 'discrete' and 'continuous' are acceptable.\nExit\n";
            return 1;
        }
    }

    if(algorithm=="LSH")
        algorithm = "LSH_Vector";
    else if(algorithm=="Frechet" && metric=="discrete")
        algorithm = "LSH_Frechet_Discrete";
    else if(algorithm=="Frechet" && metric=="continuous")
        algorithm = "LSH_Frechet_Continuous";


    //create output file
    ofstream file;
    file.open(output_file);
    file.close(); 
    cout << "Output file named: '" << output_file << "' just created." << endl;


    ifstream open_input(input_file);
    string line;
    double cur_y;
    string str;

    vector<Curve *> curves_set;
    Curve *cur_curve;

    int cnt_input=0;

    while(getline(open_input, line)){

        time_x=1;

        istringstream ss(line);
        ss >> str;
        cur_curve = new Curve(str);
        while(ss >> cur_y){
            
            curve_struct.x = time_x;
            curve_struct.y = cur_y;
            cur_curve->add_curve_coordinate(curve_struct);
            time_x++;
        }

        curves_set.push_back(cur_curve);

        cnt_input++;

        //if(cnt_input==5){

         //   break;
       // }
        
    }
    open_input.close();


    LSHCurve lsh(curves_set, L, k, w, cur_curve->get_dimension(), delta);
    lsh.hash_all_input();


    int cnt_query=0;
    ifstream open_query(query_file);
    while(getline(open_query, line)){
        
        time_x=1;

        istringstream ss(line);
        ss >> str;
        cur_curve = new Curve(str);
        while(ss >> cur_y){

            curve_struct.x = time_x;
            curve_struct.y = cur_y;
            cur_curve->add_curve_coordinate(curve_struct);
            time_x++;
        }

        lsh.print_query_results(*cur_curve, output_file, algorithm, info);
        
        //nn_info info;
        //lsh.find_nearest_neighbor(*cur_curve, info, cnt_query);
        //cout << endl;
        //lsh.NN_bruteforce(*cur_curve, info, 0);  

        cout << cnt_query+1 << " queries have been executed." << endl;
        delete cur_curve;
        cout << endl << endl << endl;

        cnt_query++;
        if(cnt_query==200){
            break;
        }

    }

    open_query.close();

    ofstream out;
    out.open(output_file, ios::app);
    out << "tApproximateAverage: " << info.total_time / cnt_query << " ms" << endl;
    out << "tTrueAverage: " << info.total_time_true / cnt_query << " ms" << endl;
    out << "MAF: " << info.maf << endl;

    return 0;
}
