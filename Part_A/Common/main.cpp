#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include "../LSH/LSH.hpp"
#include "../Hypercube/Hypercube.hpp"
#include "../Discrete_Frechet/Discrete.hpp"
#include "../Continuous_Frechet/Continuous.hpp"


int main(int argc, char *argv[]){

    int k=4, L=4, dd=3, M=20, probes=3, N=1;
    double delta=-1.0, cur_y, cur_coordinate;
    string input_file, query_file, output_file="output.txt", algorithm, metric;
    string line, str;
    int cnt_input=0, cnt_query=0, R=10;
    double time_x=1.0;

    print_result_info info;
    info.maf=INT32_MIN;
    info.total_af=0;

    random_device rand_dev;
    mt19937 gen(rand_dev());
    //W parameter (hard-coded). Explained in README
    int w=100;

    check_arguments(argc, argv, input_file, query_file, output_file, k, L, M, probes, algorithm, metric, delta);

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

    //create output file
    ofstream file;
    file.open(output_file);
    file.close(); 
    cout << "Output file named: '" << output_file << "' just created." << endl;
    ifstream open_input(input_file);


    if(algorithm=="LSH"){

        algorithm = "LSH_Vector";

        vector<Point_lsh *> points_set;
        Point_lsh *cur_point;

        while(getline(open_input, line)){

            istringstream ss(line);
            ss >> str;
            cur_point = new Point_lsh(str);
            while(ss >> cur_coordinate){

                cur_point->add_coordinate(cur_coordinate);
            }

            points_set.push_back(cur_point);
        }
        open_input.close();


        LSH lsh(points_set, L, k, w, cur_point->get_dimension(), R);
        lsh.insert_lsh();

        ifstream open_query(query_file);
        while(getline(open_query, line)){

            istringstream ss(line);
            ss >> str;
            cur_point = new Point_lsh(str);
            while(ss >> cur_coordinate){

                cur_point->add_coordinate(cur_coordinate);
            }

            lsh.print_query_results(*cur_point, output_file, algorithm, info);


            delete cur_point;

            cnt_query++;

            cout << cnt_query << " queries have been executed." << endl;
        }

        open_query.close();

    }
    else if(algorithm=="Hypercube"){

        algorithm = "Hypercube";

        vector<Point_hc *> points_set;
        Point_hc *cur_point;

        while(getline(open_input, line)){

            istringstream ss(line);
            ss >> str;
            cur_point = new Point_hc(str);
            while(ss >> cur_coordinate){

                cur_point->add_coordinate(cur_coordinate);
            }

            points_set.push_back(cur_point);
        }
        open_input.close();


        //initialising the hypercube
        Hypercube mycube=Hypercube(dd, M, probes, N, R, w, cur_point->get_dimension(), points_set);
        mycube.hypercube_insert();

        ifstream open_query(query_file);
        while(getline(open_query, line)){

            istringstream ss(line);
            ss >> str;
            cur_point = new Point_hc(str);
            while(ss >> cur_coordinate){

                cur_point->add_coordinate(cur_coordinate);
            }

            mycube.print_query_results(*cur_point, output_file, algorithm, info);

            delete cur_point;

            cnt_query++;

            cout << cnt_query << " queries have been executed." << endl;
        }

        open_query.close();
    }
    else if(algorithm=="Frechet" && metric=="discrete"){

        algorithm = "LSH_Frechet_Discrete";

        vector<Curve_discrete *> curves_set;
        Curve_discrete *cur_curve;
        curve_struct curve;

        while(getline(open_input, line)){

            time_x=1;

            istringstream ss(line);
            ss >> str;
            cur_curve = new Curve_discrete(str);
            while(ss >> cur_y){
                
                curve.x = time_x;
                curve.y = cur_y;
                cur_curve->add_curve_coordinate(curve);
                time_x++;
            }

            curves_set.push_back(cur_curve);

            cnt_input++;
        }
        open_input.close();


        Discrete lsh(curves_set, L, k, w, cur_curve->get_dimension(), delta);
        lsh.hash_all_input();


        ifstream open_query(query_file);
        while(getline(open_query, line)){
            
            time_x=1;

            istringstream ss(line);
            ss >> str;
            cur_curve = new Curve_discrete(str);
            while(ss >> cur_y){

                curve.x = time_x;
                curve.y = cur_y;
                cur_curve->add_curve_coordinate(curve);
                time_x++;
            }

            lsh.print_query_results(*cur_curve, output_file, algorithm, info);
            
            delete cur_curve;

            cnt_query++;

            cout << cnt_query << " queries have been executed." << endl;
        }

        open_query.close();

    }
    else if(algorithm=="Frechet" && metric=="continuous"){

        algorithm = "LSH_Frechet_Continuous";

        vector<Curve_cont *> curves_set;
        Curve_cont *cur_curve;

        while(getline(open_input, line)){

            istringstream ss(line);
            ss >> str;
            cur_curve = new Curve_cont(str);
            while(ss >> cur_y){
                
                cur_curve->add_curve_coordinate(cur_y);
            }

            curves_set.push_back(cur_curve);
            cnt_input++;
 
        }
        
        open_input.close();

        Continuous lsh(curves_set, k, w, cur_curve->get_dimension(), delta);
        lsh.hash_all_input();

        ifstream open_query(query_file);
        while(getline(open_query, line)){
            
            istringstream ss(line);
            ss >> str;
            cur_curve = new Curve_cont(str);
            while(ss >> cur_y){

                cur_curve->add_curve_coordinate(cur_y);
            }

            lsh.print_query_results(*cur_curve, output_file, algorithm, info);

            cout << cnt_query+1 << " queries have been executed." << endl;
        
            delete cur_curve;

            cnt_query++;
        }
        open_query.close();
    }

    ofstream out;
    out.open(output_file, ios::app);

    //convert microseconds (μs) to milliseconds (ms)
    info.total_time = info.total_time / 1000.0;
    info.total_time_true = info.total_time_true / 1000.0;
    out << "tApproximateAverage: " << info.total_time / (double)cnt_query << " ms" << endl;
    out << "tTrueAverage: " << info.total_time_true / (double)cnt_query << " ms" << endl;
    out << "MAF: " << info.maf << endl;
    out << "[Extra] Average Approximation Factor: " << info.total_af / cnt_query << endl;

    return 0;
}