#include "Curve_discrete.hpp"

Curve_discrete::Curve_discrete(string id){

    this->item_id = id;
    this->dimension=0;
}
Curve_discrete::~Curve_discrete(){}


int Curve_discrete::add_curve_coordinate(curve_struct c){

    polygonal_curve.push_back(c);
    dimension = dimension + 1;
    return 0;
}

int Curve_discrete::erase_coordinate(int index){

    polygonal_curve.erase(polygonal_curve.begin()+index);
    dimension = dimension - 1;
    return 0;
}

vector<curve_struct> Curve_discrete::get_curve(){

    return polygonal_curve;
}


string Curve_discrete::get_item_id(){

    return item_id;
}


void Curve_discrete::print_curves(){

    for(int i=0;i<polygonal_curve.size()-1;i++){

        cout << "(" << polygonal_curve[i].x << ", " << polygonal_curve[i].y << "),\t";
    }

    int index=polygonal_curve.size() - 1;
    cout << "(" << polygonal_curve[index].x << ", " << polygonal_curve[index].y << ")";
    cout << endl << endl;
}


void Curve_discrete::set_init_curve(vector<curve_struct> curve){

     this->init_curve = curve;
}

vector<curve_struct> Curve_discrete::get_init_curve(){

    return this->init_curve;
}

void Curve_discrete::set_dimension(int d){

     this->dimension = d;
}

int Curve_discrete::get_dimension(){

    return this->dimension;
}