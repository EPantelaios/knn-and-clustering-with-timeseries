#include "Point_lsh.hpp"


Point_lsh::Point_lsh(string id){

    this->item_id = id;
    this->dimension=0;
}
Point_lsh::~Point_lsh(){}

int Point_lsh::add_coordinate(double x){

    coordinates.push_back(x);
    dimension++;
    return 0;
}

vector<double> Point_lsh::get_coordinates(){

    return coordinates;
}


string Point_lsh::get_item_id(){

    return item_id;
}


void Point_lsh::print_coordinates(Point_lsh x){

    for(int i=0;i<x.coordinates.size();i++)
        cout << x.coordinates[i] << " ";
    cout << endl << endl;
}


void Point_lsh::set_dimension(int d){

     this->dimension = d;
}

int Point_lsh::get_dimension(){

    return this->dimension;
}

void Point_lsh::add_ID(unsigned long long int id){

    this->ID.push_back(id);
}

unsigned long long int Point_lsh::get_ID(int i){

    return this->ID.at(i);
}