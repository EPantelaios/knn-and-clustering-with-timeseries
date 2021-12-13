#include "Point.hpp"


Point::Point(string id){

    this->item_id = id;
    this->dimension=0;
}

Point::~Point(){}

int Point::add_coordinate(double x){

    coordinates.push_back(x);
    dimension++;
    return 0;
}


vector<double> Point::get_coordinates(){

    return coordinates;
}


string Point::get_item_id(){

    return item_id;
}


void Point::print_coordinates(){

    for(int i=0;i<this->coordinates.size();i++)
        cout << this->coordinates[i] << " ";
    cout << endl << endl;
}



void Point::print_origin_curve(){

    for(int i=0;i<origin_curve.size()-1;i++){

        cout << "(" << origin_curve[i].x << ", " << origin_curve[i].y << "),\t";
    }

    int index=origin_curve.size() - 1;
    cout << "(" << origin_curve[index].x << ", " << origin_curve[index].y << ")";
    cout << endl << endl;
}





void Point::set_dimension(int d){

     this->dimension = d;
}

int Point::get_dimension(){

    return this->dimension;
}

void Point::set_ID(unsigned long long int id){

    this->ID = id;
}

unsigned long long int Point::get_ID(){

    return this->ID;
}


vector<curve> Point::get_origin_curve(){

    return this->origin_curve;
}


void Point::set_origin_curve(vector<curve> curve){

     this->origin_curve = curve;
}

