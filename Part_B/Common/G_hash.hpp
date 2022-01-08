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

    public:
        G_hash(int k, int w, int d, int table_size);
        
        unsigned long long int G_compute(Point_frechet p);
};