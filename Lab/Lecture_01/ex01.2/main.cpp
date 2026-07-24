#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include "random.h"

using namespace std;

int main(int argc, char *argv[]){
    Random rnd;
    int seed[4];
    int p1, p2;

    // Files
    ifstream Primes("Primes");
    if (Primes.is_open()){
        Primes >> p1 >> p2;
    } else cerr << "\nUnable to open 'Primes' file.\n" << endl;
    Primes.close();

    ifstream input("seed.in");
    string property;
    if (input.is_open()){
        while(!input.eof()){
            input >> property;
            if (property == "RANDOMSEED"){
                input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
                rnd.SetRandom(seed,p1,p2);
            }
        }
        input.close();
    } else cerr << "\nUnable to open 'seed.in' file.\n" << endl;

    // Parameters
    const int M=10000;
    const double lambda=1.;
    const double gamma=1.;
    const double mu=0.;

    vector <int> N = {1,2,10,100};

    vector <ofstream> out;
    vector <string> filename = {"out_N1.txt", "out_N2.txt", "out_N10.txt", "out_N100.txt"};

    for(size_t j = 0; j<N.size(); j++){
        out.emplace_back(filename[j]);
        if(!out[j].is_open()){
            cerr << "\nUnable to open '" << filename[j] << "'\n" << endl;
            return -1;
        }
    }

    // Calculation
    for(size_t j=0; j<N.size(); j++){
        for(int i=0; i<M; i++){
            double sum_unif=0., sum_exp=0., sum_cauchy_lorentz=0.;

            // Generation and averages computation
            for(int k=0; k<N[j]; k++){
                sum_unif += rnd.Rannyu();
                sum_exp += rnd.Exponential(lambda);
                sum_cauchy_lorentz += rnd.CauchyLorentz(gamma, mu);
            }
            sum_unif /= N[j];
            sum_exp /= N[j];
            sum_cauchy_lorentz /= N[j];

            out[j] << sum_unif << " " << sum_exp << " " << sum_cauchy_lorentz << "\n";
        }
        out[j].close();
    }

    rnd.SaveSeed();
    return 0;
}