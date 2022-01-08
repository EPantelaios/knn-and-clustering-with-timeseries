#include "Point_discrete.hpp"


Point_discrete::Point_discrete(string id){

    this->item_id = id;
    this->dimension=0;
}

Point_discrete::~Point_discrete(){}

int Point_discrete::add_coordinate(double x){

    coordinates.push_back(x);
    dimension++;
    return 0;
}


vector<double> Point_discrete::get_coordinates(){

    return coordinates;
}


string Point_discrete::get_item_id(){

    return item_id;
}


void Point_discrete::print_coordinates(){

    for(int i=0;i<this->coordinates.size();i++)
        cout << this->coordinates[i] << " ";
    cout << endl << endl;
}



void Point_discrete::print_origin_curve(){

    for(int i=0;i<origin_curve.size()-1;i++){

        cout << "(" << origin_curve[i].x << ", " << origin_curve[i].y << "),\t";
    }

    int index=origin_curve.size() - 1;
    cout << "(" << origin_curve[index].x << ", " << origin_curve[index].y << ")";
    cout << endl << endl;
}



void Point_discrete::set_dimension(int d){

     this->dimension = d;
}

int Point_discrete::get_dimension(){

    return this->dimension;
}

void Point_discrete::set_ID(unsigned long long int id){

    this->ID = id;
}

unsigned long long int Point_discrete::get_ID(){

    return this->ID;
}


vector<curve_struct> Point_discrete::get_origin_curve(){

    return this->origin_curve;
}


void Point_discrete::set_origin_curve(vector<curve_struct> curve){

     this->origin_curve = curve;
}

