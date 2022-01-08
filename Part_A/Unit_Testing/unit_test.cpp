#include <vector>
#include <fstream>
#include <sstream>
#include "acutest/include/acutest.h"
#include "../Discrete_Frechet/Discrete.hpp"
#include "../Continuous_Frechet/Continuous.hpp"
#include "../Common/utils.hpp"

void test_discrete_frechet_distance(void){

    int k=4, L=5, w=50;
    double delta, cur_y;
    string input_file, query_file;
    string line, str;
    int cnt_input=0, cnt_query=0;
    double time_x=1.0, tmp_distance=0;

    input_file="nasd_input.csv";
    query_file="nasd_query.csv";
    delta = 1;

    ifstream open_input(input_file);
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

        if(cnt_input==10){

            break;
        }
    }
    open_input.close();

    TEST_CHECK(curves_set.size() == (long unsigned int)cnt_input);

    Discrete lsh(curves_set, L, k, w, cur_curve->get_dimension(), delta);

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
        cnt_query++;

        if(cnt_query==1){
            break;
        }
    }

    open_query.close();

    for(int i=0;i<cnt_input;i++){

        tmp_distance = lsh.frechet_distance(cur_curve->get_curve(), curves_set[i]->get_curve());

        TEST_CHECK(tmp_distance = 16.4265);
        TEST_CHECK(tmp_distance = 26.8497);
        TEST_CHECK(tmp_distance = 19.4428);
        TEST_CHECK(tmp_distance = 44.6719);
        TEST_CHECK(tmp_distance = 13.5);
        TEST_CHECK(tmp_distance = 427.329);
        TEST_CHECK(tmp_distance = 49.7201);
        TEST_CHECK(tmp_distance = 21.3231);
        TEST_CHECK(tmp_distance = 49.5318);
        TEST_CHECK(tmp_distance = 26.72);
    }
 
   delete cur_curve;
}




void test_discrete_frechet_insert(void){

    int k=4, L=5, w=50;
    double delta, cur_y;
    string input_file, query_file, output_file="output.txt", algorithm;
    string line, str;
    int cnt_input=0, cnt_query=0;
    double time_x=1.0;
    print_result_info info;
    info.maf=INT32_MIN;
    info.total_af=0;

    input_file="nasd_input.csv";
    query_file="nasd_query.csv";
    algorithm="testing_process";
    delta = 1;

    ifstream open_input(input_file);
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

    TEST_CHECK(curves_set.size() == (long unsigned int)cnt_input);

    Discrete lsh(curves_set, L, k, w, cur_curve->get_dimension(), delta);

    int size = cur_curve->get_dimension() * 2;

    TEST_CHECK(lsh.hash_all_input() == 0);

    for(int i=0;i<L;i++){

        vector<Point_discrete>  points = lsh.get_points_set(i);

        for(vector<Point_discrete>::iterator it = points.begin(); it != points.end(); it++){
            
            TEST_CHECK(it->get_dimension() == size);
        }
    }

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

        TEST_CHECK(lsh.print_query_results(*cur_curve, output_file, algorithm, info) == 1);

        delete cur_curve;

        cnt_query++;

        if(cnt_query==1){
            break;
        }
    }

    open_query.close();
    TEST_CHECK(info.maf <= info.total_af);
}




void test_number_of_points_continuous(void){

    double cur_y;
    string input_file, query_file;
    string line, str;
    int cnt_input=0;

    input_file="nasd_input.csv";
    query_file="nasd_query.csv";

    ifstream open_input(input_file);
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

        TEST_CHECK(curves_set[cnt_input]->get_curve().size() == 730);

        cnt_input++;
    }

    open_input.close();

    TEST_CHECK(curves_set.size() == (long unsigned int)cnt_input);

    for(vector<Curve_cont *>::iterator it=curves_set.begin(); it != curves_set.end(); it++){
    
        delete *it;
    }
}





TEST_LIST = {

    {"frechet_distance", test_discrete_frechet_distance},
    {"discrete_frechet", test_discrete_frechet_insert},
    {"number_of_points_continuous_frechet", test_number_of_points_continuous},
    {NULL, NULL}     /* zeroed record marking the end of the list */
};