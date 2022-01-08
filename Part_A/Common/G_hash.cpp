#include "G_hash.hpp"

G_hash::G_hash(){}

G_hash::G_hash(int k, int w, int d, int table_size){

    this->k=k;
    this->w=w;
    this->dimension=d;
    this->M=pow(2, 32) - 5;
    this->table_size = table_size;
    this->ID=0;

    random_device rand_dev;
    mt19937 gen(rand_dev());
    uniform_int_distribution<int> x(0, INT32_MAX);

    for(int i=0;i<k;i++){

        H_hash h(w, dimension);
        h_hash.push_back(h);
        r.push_back(x(gen));
    }
}


unsigned long long int G_hash::ID_compute(Point_lsh p){

    unsigned long long int sum=0;
    unsigned long long int tmp=0;

    for(int i=0;i<k;i++){

        tmp = h_hash[i].h_compute(p) * r[i];
        sum += euclidean_mod(tmp, M);
    }

    this->ID = euclidean_mod(sum, M);

    return ID;
}


unsigned long long int G_hash::ID_compute(Point_hc p){

    unsigned long long int sum=0;
    unsigned long long int tmp=0;

    for(int i=0;i<k;i++){

        tmp = h_hash[i].h_compute(p) * r[i];
        sum += euclidean_mod(tmp, M);
    }

    this->ID = euclidean_mod(sum, M);

    return ID;
}


unsigned long long int G_hash::ID_compute(Point_discrete p){

    unsigned long long int sum=0;
    unsigned long long int tmp=0;

    for(int i=0;i<k;i++){

        tmp = h_hash[i].h_compute(p) * r[i];
        sum += euclidean_mod(tmp, M);
    }

    this->ID = euclidean_mod(sum, M);

    return ID;
}


unsigned long long int G_hash::ID_compute(Point_cont p){

    unsigned long long int sum=0;
    unsigned long long int tmp=0;

    for(int i=0;i<k;i++){

        tmp = h_hash[i].h_compute(p) * r[i];
        sum += euclidean_mod(tmp, M);
    }

    this->ID = euclidean_mod(sum, M);

    return ID;
}

unsigned long long int G_hash::get_ID(){

    return this->ID;
}