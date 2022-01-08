#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include <cstring>
#include <random>
#include "../LSH/Point_lsh.hpp"
#include "../Hypercube/Point_hc.hpp"
#include "../Discrete_Frechet/Point_discrete.hpp"
#include "../Continuous_Frechet/Point_cont.hpp"


using namespace std;

class H_hash{

    private:
        vector<double> v;
        float t;
        int w;
        int dimension;

    public:
        H_hash(int w, int d);
        ~H_hash();
        int h_compute(Point_lsh p);
        int h_compute(Point_hc p);
        int h_compute(Point_discrete p);
        int h_compute(Point_cont p);
};