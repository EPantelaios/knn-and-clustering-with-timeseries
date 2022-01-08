#include "Point_hc.hpp"

using namespace std;


Point_hc::Point_hc(string id){

    this->item_id = id;
    this->dimension=0;
}

Point_hc::~Point_hc(){}

int Point_hc::add_coordinate(double x){

    coordinates.push_back(x);
    dimension++;
    return 0;
}

vector<double> Point_hc::get_coordinates(){

    return coordinates;
}

double Point_hc::get_a_coordinate(int i){

    return coordinates[i];
}


string Point_hc::get_item_id(){

    return item_id;
}


void Point_hc::print_coordinates(){

    for(int i=0;i<this->coordinates.size();i++)
        cout << this->coordinates[i] << " ";
    cout << endl << endl;
}


void Point_hc::set_dimension(int d){

     this->dimension = d;
}

int Point_hc::get_dimension(){

    return this->dimension;
}