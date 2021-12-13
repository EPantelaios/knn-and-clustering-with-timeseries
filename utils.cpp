#include "utils.hpp"


int check_arguments(int argc, char **argv, string &input_file, string &query_file, string &output_file, int &k,
                    int &L, int &M, int &probes, string &algorithm, string &metric, double &delta){

    char *arg=NULL;
    int cnt=0;

    if(argc<7){

        cout << "'input_file', 'query_file' and 'algorithm_method' must be included\nExit.\n";
        exit(1);
    }

    //Edit the parameters given at program
    while(--argc){

        arg = *++argv;
        
        if(!strcmp(arg, "-i")){	

            if(argc > 1 && --argc) 
                input_file = *++argv;
            
            cnt++;
        }
        else if(!strcmp(arg, "-q")){	

            if(argc > 1 && --argc) 
                query_file = *++argv;

            cnt++;
        }
        else if(!strcmp(arg, "-o")){	

            if(argc > 1 && --argc) 
                output_file = *++argv;
        }
        else if(!strcmp(arg, "-k")){	

            if(argc > 1 && --argc) 
                k = atoi(*++argv);
        }
        else if(!strcmp(arg, "-L")){	

            if(argc > 1 && --argc) 
                L = atoi(*++argv);
        }
        else if(!strcmp(arg, "-M")){	

            if(argc > 1 && --argc) 
                M = atoi(*++argv);
        }
        else if(!strcmp(arg, "-probes")){	

            if(argc > 1 && --argc) 
                probes = atoi(*++argv);
        }
        else if(!strcmp(arg, "-algorithm")){	

            if(argc > 1 && --argc) 
                algorithm = *++argv;
            
            cnt++;
        }
        else if(!strcmp(arg, "-metric")){	

            if(argc > 1 && --argc) 
                metric = *++argv;
            
            cnt++;
        }
        else if(!strcmp(arg, "-delta")){	

            if(argc > 1 && --argc) 
                delta = stod(*++argv);
            
            cnt++;
        }
        else{

            cout << "FALSE arguments. Please enter the correct arguments (-i <input file> -q <query file>\
                     -algorithm <LSH or Hypercube or Frechet)\nExit.\n";
            exit(1);
        }
    }


    if(cnt!=3 && cnt!=5){

        cout << "Please enter the required arguments.\n \
        '-i <input_file>', '-c <cluster.conf>' and '-m <method>' must be included\n.Exit.\n";

        exit(1);
    }
    
    return 0;

    
}


bool file_exist(const string &filename){

    if (FILE *file = fopen(filename.c_str(), "r")){
        fclose(file);
        return true;
    }else{
        return false;
    }
}

unsigned long long int euclidean_mod(long long int a, unsigned long long int b){ //for us b is always positive
    
    assert(b != 0);
    unsigned long long int r = a % b;

    if(r>=0)
        return r;
    else
        return r + b;
}


double euclidean_distance(vector<double> x, vector<double> y){

    double res=0.0;

    if(x.size() != y.size()){

        assert(1==0);
        return -1.0;
    }

    for(int i=0;i<x.size();i++){

        res += (x[i] - y[i]) * (x[i] - y[i]);
    }

    return sqrt(res);
}