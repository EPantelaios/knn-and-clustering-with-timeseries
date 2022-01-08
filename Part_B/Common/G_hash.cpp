#include "G_hash.hpp"
#include "utils.hpp"

G_hash::G_hash(int k, int w, int d, int table_size){

    this->k=k;
    this->w=w;
    this->dimension=d;
    this->M=pow(2, 32) - 5;
    this->table_size = table_size;

    random_device rand_dev;
    mt19937 gen(rand_dev());
    uniform_int_distribution<int> x(0, INT32_MAX);

    for(int i=0;i<k;i++){

        H_hash h(w, dimension);
        h_hash.push_back(h);

        r.push_back(x(gen));
    }
}


unsigned long long int G_hash::G_compute(Point_frechet p){

    unsigned long long int sum=0;
    unsigned long long int tmp=0;
    
    for(int i=0;i<k;i++){

        tmp = h_hash[i].h_compute(p) * r[i];

        sum += euclidean_mod(tmp, M);
    }

    sum = euclidean_mod(sum, M);

    return sum % table_size;
}