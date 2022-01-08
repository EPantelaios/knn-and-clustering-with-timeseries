#pragma once
#include "H_hash.hpp"
#include "utils.hpp"

using namespace std;

class G_hash{

    private:
        vector<H_hash> h_hash;
        vector<unsigned int> r;
        int k;
        int w;
        int dimension;
        int table_size;
        unsigned long long int M; //constant value 2^32 -5
        unsigned long long int ID;//Object ID (slides page 22)

    public:
        G_hash();
        G_hash(int k, int w, int d, int table_size);
        unsigned long long int ID_compute(Point_lsh p);
        unsigned long long int ID_compute(Point_hc p);
        unsigned long long int ID_compute(Point_discrete p);
        unsigned long long int ID_compute(Point_cont p);
        unsigned long long int get_ID();

};