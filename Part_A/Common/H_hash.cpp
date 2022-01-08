#include "H_hash.hpp"

H_hash::H_hash(int w, int d){
    
    this->w=w;
    this->dimension=d;

    random_device rand_dev;
    mt19937 gen(rand_dev());
    uniform_real_distribution<float> x(0,w);
    t = x(gen);

    normal_distribution<double> n_d(0, 1);

    for(int i=0;i<d;i++){

        v.push_back(n_d(rand_dev));
    }
}

H_hash::~H_hash(){}


int H_hash::h_compute(Point_lsh p){

    double inner_product=0;

    for(int i=0;i<dimension;i++){

        inner_product += p.get_coordinates().at(i) * v[i];
    }

    return floor((inner_product + t)/w);
}

int H_hash::h_compute(Point_hc p){

    double inner_product=0;

    for(int i=0;i<dimension;i++){

        inner_product += p.get_coordinates().at(i) * v[i];
    }

    return floor((inner_product + t)/w);
}

int H_hash::h_compute(Point_discrete p){

    double inner_product=0;

    for(int i=0;i<dimension;i++){

        inner_product += p.get_coordinates().at(i) * v[i];
    }

    return floor((inner_product + t)/w);
}

int H_hash::h_compute(Point_cont p){

    double inner_product=0;

    for(int i=0;i<dimension;i++){

        inner_product += p.get_coordinates().at(i) * v[i];
    }

    return floor((inner_product + t)/w);
}