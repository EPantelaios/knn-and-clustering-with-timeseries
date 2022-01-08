#include <iostream>
#include <vector>
#include <cmath>
#include <cstring>
#include <random>
#include "Point_frechet.hpp"

using namespace std;

class H_hash{

    private:
        vector<double> v;
        float t;
        int w;
        int dimension;

    public:
        H_hash(int w, int d);

        int h_compute(Point_frechet p);
};