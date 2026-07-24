#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include "random.h"

using namespace std;

double error(const vector<double>& AV, const vector<double>& AV2, int n) {
    // This function calculates the statistical error of the mean
    if (n == 0) return 0.0;
    else return sqrt((AV2[n] - AV[n] * AV[n]) / n);
}

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
    int M=10000;        // throws
    int N=100;          // blocks
    int Nt=M/N;         // number of throws for each block
    double L=0.8;       // length of the needle
    double d=1.0;       // spacing between the lines
    double alpha, y0;   // orientation and height of the center

    vector <double> avg, avg2;
    vector <double> sum, sum2, err;

    // Buffon's Simulation
    for(int i=0; i<N; i++){
        int N_hit=0;
        for(int j=0; j<Nt; j++){
            alpha = rnd.Angle();
            y0 = rnd.Rannyu(0,1);

            double y_up = y0 + (L/2.) * sin(alpha);
            double y_down = y0 - (L/2.) * sin(alpha);

            // Hit condition
            if (y_up >= d || y_down <= 0.) N_hit++;
        }
        if (N_hit > 0) {
            avg.push_back((2 * L * Nt) / (d * N_hit));
        } else {
            avg.push_back(0); 
        }       
        avg2.push_back(avg[i]*avg[i]);    
    }

    for(int i=0; i<N; i++){
        sum.push_back(0);
        sum2.push_back(0);

        for(int j=0; j<i+1; j++){
            sum[i] += avg[j];
            sum2[i] += avg2[j];
        }
        sum[i] /= (i+1);
        sum2[i] /= (i+1);
        err.push_back(error(sum,sum2,i));
    }

    // Output File
    ofstream out_avg;
    out_avg.open("out_avg.data");
    for(int i=0; i<N; i++){
        out_avg << (i+1)*Nt << "\t" << sum[i] << "\t" << err[i] << endl;
    }
    out_avg.close();

    // Results on Terminal
    cout << "\n\nNumber of blocks = " << N << endl;
    cout << "Number of throws per block = " << Nt << endl;
    cout << "Number of total throws = " << M << endl;
    cout << "Best value of pi = " << sum[N-1] << " +-" << err[N-1] << "\n\n" << endl;

    rnd.SaveSeed();
    return 0;
    
}