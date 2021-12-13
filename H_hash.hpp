#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include <cstring>
#include <random>
#include "Point.hpp"

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
        int h_compute(Point p);
};