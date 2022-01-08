#include <vector>
#include <fstream>
#include <sstream>
#include "acutest/include/acutest.h"
#include "./Curve/Clustering_Frechet.hpp"
#include "./Common/utils.hpp"

void test_back_track_path()
{
    back_track_path *the_path;
    //check path_length
    //check path values
    //3 tests
    string input_file="test_input1.csv";
    ifstream open_input(input_file);
    double cur_y, time_x=1.0;
    vector<Point_frechet *> origin_curves;
    Point_frechet *cur_curve;
    string line, str;
    curve_struct curve;
    //open and read input file
    while(getline(open_input, line)){

        time_x=1;
        istringstream ss(line);
        ss >> str;
        cur_curve = new Point_frechet(str);
        
        while(ss >> cur_y){
                
            curve.x = time_x;
            curve.y = cur_y;
            cur_curve->add_curve_coordinate(curve);
            time_x++;
        }
        
        origin_curves.push_back(cur_curve);

    }
    open_input.close();
    int i=0;
    int result[4]={4,2,7,3};
    int result_optimal[4][8]={{2,0,2,1},{2,0},{1,1,1,2,2,2,0},{2,2,1}};
    int j=0;
    while(i<=6)
    {
        the_path=frechet_distance_and_path(origin_curves[i]->get_curve_address(),origin_curves[i+1]->get_curve_address());
        TEST_CHECK(the_path->t==result[j]);
        //cout << "The path:" << the_path->t <<" result" << result[j] << endl;
        for(int k=0;k<the_path->t;k++)
        {
            TEST_CHECK(the_path->path[k]==result_optimal[j][k]); //tests if it has an optimal path
        }
        
        j++;
        i=i+2;

    }
}

void test_mean_curve()
{//tests if the mean curve is returned
    back_track_path *the_path;
    //check path_length
    //check path values
    //3 tests
    string input_file="test_input1.csv";
    ifstream open_input(input_file);
    double cur_y, time_x=1.0;
    vector<Point_frechet *> origin_curves;
    Point_frechet *cur_curve;
    string line, str;
    curve_struct curve;
    //open and read input file
    while(getline(open_input, line)){

        time_x=1;
        istringstream ss(line);
        ss >> str;
        cur_curve = new Point_frechet(str);
        
        while(ss >> cur_y){
                
            curve.x = time_x;
            curve.y = cur_y;
            cur_curve->add_curve_coordinate(curve);
            time_x++;
        }
        
        origin_curves.push_back(cur_curve);

    }
    open_input.close();
    int i=0;
    int j=0;
    vector <curve_struct> mean_curve;
    double values_x[4][8]={{1,1.5,2.5,3,4},{1,1.5,2.5},{1,1.5,2.5,3.5,4.5,5,5.5,6},{1,1.5,2.5,3.5}};
    double values_y[4][8]={{5,6.5,12.5,12,9},{5.5,6,2.5},{5.5,6.5,16.5,30,17,20.5,19,18},{10.35,9.6,6.9,9.7}};
    while(i<=6)
    {
        the_path=frechet_distance_and_path(origin_curves[i]->get_curve_address(),origin_curves[i+1]->get_curve_address());
        mean_curve=find_mean_curve(the_path,origin_curves[i]->get_curve(),origin_curves[i+1]->get_curve());   
        for(int k=0;k<mean_curve.size();k++)
        {
            TEST_CHECK(values_x[j][k]==mean_curve[k].x);
            TEST_CHECK(values_y[j][k]==mean_curve[k].y);
        }  
        j++;
        i=i+2;

    }

}


TEST_LIST = {

    {"back_track_path", test_back_track_path},
    {"mean_curve", test_mean_curve},
    {NULL, NULL}     /* zeroed record marking the end of the list */
};