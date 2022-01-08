#include "utils.hpp"
#include <cmath>


int check_arguments(int argc, char **argv, string &input_file, string &conf_file, string &output_file,
                    string &update_method, string &assignment_method, int &complete_param, int &silhouette_param){

    char *arg=NULL;
    int cnt=0;

    if(argc<7){

        cout << "'-i <input_file>', '-c <cluster.conf>' and '-m <method>' must be included\nExit.\n";
        exit(1);
    }

    if(argc>14){

        cout << "Please enter the right arguments\nExit.\n";
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
        else if(!strcmp(arg, "-c")){	

            if(argc > 1 && --argc) 
                conf_file = *++argv;

            cnt++;
        }
        else if(!strcmp(arg, "-o")){	

            if(argc > 1 && --argc) 
                output_file = *++argv;
        }
        else if(!strcmp(arg, "-complete")){	
            
            complete_param = 1;
        }
        else if(!strcmp(arg, "-silhouette")){	
            
            silhouette_param = 1;
        }
        else if(!strcmp(arg, "-update")){	

            if(argc > 1 && --argc) 
                update_method = *++argv;
            cnt++;
        }
        else if(!strcmp(arg, "-assignment")){	

            if(argc > 1 && --argc) 
                assignment_method = *++argv;
            cnt++;
        }
        else{

            cout << "FALSE arguments. Please enter the correct arguments \
                    ('-i <input_file>', '-c <cluster.conf>', '-update <Mean Frechet or Mean Vector>' and \
                    '-assignment <Classic or LSH or Hypercube>').\nExit.\n";
            exit(1);
        }
    }

    if(cnt!=4){

        cout << "Please enter the required arguments.\n \
                 '-i <input_file>', '-c <cluster.conf>', '-update <Mean Frechet or Mean Vector>' and \
                 '-assignment <Classic or LSH or Hypercube>' must be included\n.Exit.\n";

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

    if(x.size() != y.size())
        return -1.0;

    for(int i=0;i<x.size();i++){

        res += (x[i] - y[i]) * (x[i] - y[i]);
    }

    return sqrt(res);
}