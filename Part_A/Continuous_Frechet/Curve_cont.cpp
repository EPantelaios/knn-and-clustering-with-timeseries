#include "Curve_cont.hpp"

Curve_cont::Curve_cont(){}


Curve_cont::Curve_cont(string id){

    this->item_id = id;
    this->dimension=0;
}

Curve_cont::~Curve_cont(){}

int Curve_cont::add_curve_coordinate(double c){

    polygonal_curve.push_back(c);
    dimension = dimension + 1;
    return 0;
}

int Curve_cont::set_curve_coordinate(int index, double c){

    polygonal_curve[index] = c;
    return 0;
}

int Curve_cont::erase_coordinate(int index){

    polygonal_curve.erase(polygonal_curve.begin()+index);
    dimension = dimension - 1;
    return 0;
}

vector<double> Curve_cont::get_curve(){

    return polygonal_curve;
}

vector<double> *Curve_cont::get_curve_address(){

    return &polygonal_curve;
}


string Curve_cont::get_item_id(){

    return item_id;
}


void Curve_cont::print_curves(){

    for(int i=0;i<polygonal_curve.size();i++){

        cout << i+1 << "   " << polygonal_curve[i] << endl;
    }

    cout << endl << endl;
}

void Curve_cont::print_filtered_curves(){

    for(int i=0;i<filter_curve.size();i++){

        cout << i+1 << "   " << filter_curve[i] << endl;
    }

    cout << endl << endl;
}


void Curve_cont::set_filter_curve(vector<double> curve){

     this->filter_curve = curve;
}

vector<double> Curve_cont::get_filtered_curve(){

    return this->filter_curve;
}

void Curve_cont::set_dimension(int d){

     this->dimension = d;
}

int Curve_cont::get_dimension(){

    return this->dimension;
}