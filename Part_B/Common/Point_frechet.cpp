#include "Point_frechet.hpp"

Point_frechet::Point_frechet(string id){

    this->item_id = id;
    this->dimension=0;
    
    this->min_distance_from_centroids=INT32_MAX;
    this->square_dist=0.0;
    this->sum_square_dist=0.0;

    this->marked=-1;
    this->cur_distance=INT32_MAX;
}


int Point_frechet::add_coordinate(double x){

    coordinates.push_back(x);
    dimension++;
    return 0;
}


vector<double> Point_frechet::get_coordinates(){

    return coordinates;
}


double Point_frechet::get_a_coordinate(int i){

    return coordinates[i];
}


string Point_frechet::get_item_id(){

    return item_id;
}


void Point_frechet::print_coordinates(){

    for(int i=0;i<this->coordinates.size();i++)
        cout << this->coordinates[i] << " ";
    cout << endl;
}


void Point_frechet::print_coordinates_output(ofstream &out){

    for(int i=0;i<this->coordinates.size();i++)
        out << this->coordinates[i] << " ";
}



void Point_frechet::set_dimension(int d){

     this->dimension = d;
}


int Point_frechet::get_dimension(){

    return this->dimension;
}


void Point_frechet::set_min_distance(float i){

    this->min_distance_from_centroids = i;
}


float Point_frechet::get_min_distance(){

    return this->min_distance_from_centroids;
}



void Point_frechet::set_square_dist(float i){

    this->square_dist = i;
}


float Point_frechet::get_square_dist(){

    return this->square_dist;
}


void Point_frechet::set_sum_square_dist(float i){

    this->sum_square_dist += i;
}


float Point_frechet::get_sum_square_dist(){

    return this->sum_square_dist;
}


float Point_frechet::init_sum_square_distance(){

    return this->sum_square_dist = 0.0;
}


void Point_frechet::set_cur_distance(double i){

    this->cur_distance = i;
}


double Point_frechet::get_cur_distance(){

    return this->cur_distance;
}


void Point_frechet::set_marked(int i){

    this->marked = i;
}


int Point_frechet::get_marked(){

    return this->marked;
}



curve_struct Point_frechet::get_curve_coordinate(int i){

    return this->origin_curve[i];
}


vector<curve_struct> Point_frechet::get_curve(){

    return this->origin_curve;
}


vector<curve_struct> * Point_frechet::get_curve_address(){

    return &this->origin_curve;
}



void Point_frechet::set_curve(vector<curve_struct> curve){

     this->origin_curve = curve;
}



int Point_frechet::add_curve_coordinate(curve_struct c){

    origin_curve.push_back(c);
    dimension = dimension + 1;
    return 0;
}



int Point_frechet::erase_curve_coordinate(int index){

    origin_curve.erase(origin_curve.begin()+index);
    dimension = dimension - 1;
    return 0;
}



void Point_frechet::print_curves(){

    for(int i=0;i<origin_curve.size()-1;i++){

        cout << "(" << origin_curve[i].x << ", " << origin_curve[i].y << "),\t";
    }

    int index=origin_curve.size() - 1;
    cout << "(" << origin_curve[index].x << ", " << origin_curve[index].y << ")";
    cout << endl << endl;
}



void Point_frechet::print_curves_output(ofstream &out){

    for(int i=0;i<origin_curve.size()-1;i++){

        out << "(" << origin_curve[i].x << ", " << origin_curve[i].y << "),\t";
    }

    int index=origin_curve.size() - 1;
    out << "(" << origin_curve[index].x << ", " << origin_curve[index].y << ")";
    cout << endl;
}
