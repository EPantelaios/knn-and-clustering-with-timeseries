#include "Point_cont.hpp"

Point_cont::Point_cont(){}

Point_cont::Point_cont(string id){

    this->item_id = id;
    this->dimension=0;
}

Point_cont::~Point_cont(){}

int Point_cont::add_coordinate(double x){

    coordinates.push_back(x);
    dimension++;
    return 0;
}


vector<double> Point_cont::get_coordinates(){

    return coordinates;
}


string Point_cont::get_item_id(){

    return item_id;
}


void Point_cont::print_coordinates(){

    for(int i=0;i<this->coordinates.size();i++)
        cout << i+1 << "    " << this->coordinates[i] << endl;
    cout << endl << endl;
}


void Point_cont::print_filtered_curve(){

    for(int i=0;i<this->filtered_curve.size();i++){

        cout << i+1 << "    " << this->filtered_curve[i] << endl;
    }
    cout << endl << endl;
}


void Point_cont::set_dimension(int d){

     this->dimension = d;
}

int Point_cont::get_dimension(){

    return this->dimension;
}

void Point_cont::set_ID(unsigned long long int id){

    this->ID = id;
}

unsigned long long int Point_cont::get_ID(){

    return this->ID;
}


vector<double> Point_cont::get_filtered_curve(){

    return this->filtered_curve;
}


void Point_cont::set_filtered_curve(vector<double> curve){

     this->filtered_curve = curve;
}

