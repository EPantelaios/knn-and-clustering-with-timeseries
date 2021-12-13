#include "Curve.hpp"


Curve::Curve(string id){

    this->item_id = id;
    this->dimension=0;
}
Curve::~Curve(){}


int Curve::add_curve_coordinate(curve c){

    polygonal_curve.push_back(c);
    dimension = dimension + 1;
    return 0;
}

int Curve::erase_coordinate(int index){

    polygonal_curve.erase(polygonal_curve.begin()+index);
    dimension = dimension - 1;
    return 0;
}

vector<curve> Curve::get_curve(){

    return polygonal_curve;
}


string Curve::get_item_id(){

    return item_id;
}


void Curve::print_curves(){

    for(int i=0;i<polygonal_curve.size()-1;i++){

        cout << "(" << polygonal_curve[i].x << ", " << polygonal_curve[i].y << "),\t";
    }

    int index=polygonal_curve.size() - 1;
    cout << "(" << polygonal_curve[index].x << ", " << polygonal_curve[index].y << ")";
    cout << endl << endl;
}


void Curve::set_init_curve(vector<curve> curve){

     this->init_curve = curve;
}

vector<curve> Curve::get_init_curve(){

    return this->init_curve;
}

void Curve::set_dimension(int d){

     this->dimension = d;
}

int Curve::get_dimension(){

    return this->dimension;
}