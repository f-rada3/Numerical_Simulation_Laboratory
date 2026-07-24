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

// Asset price funtion at time t
double S(double S0, double r, double sigma, double Z, double t){
    return S0 * exp((r-0.5*sigma*sigma)* t + sigma * Z * sqrt(t));
}

// Europian Call
double C(double r, double t, double S, double K){
    return exp(-r*t) * max(0., S-K);
}

// Europian Put 
double P(double r, double t, double S, double K){
    return exp(-r*t) * max(0., K-S);
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
    int M=100000;           // throws
    int N=100;              // blocks
    int L=M/N;              // number of throws for each block
    double S0=100;          // initial asset price
    double K=100;           // strike price
    double T=1.;            // delivery time
    double r=0.1;           // risk-free interest rate
    double sigma=0.25;      // volatility

    // Direct
    vector <double> avg1_call, avg1_put, avg2_call, avg2_put;
    vector <double> sum_call, sum2_call, err_call, sum_put, sum2_put, err_put;

    // Discretize
    vector <double> avg1_call_disc, avg1_put_disc, avg2_call_disc, avg2_put_disc;
    vector <double> sum_call_disc, sum2_call_disc, err_call_disc, sum_put_disc, sum2_put_disc, err_put_disc;

    // Calculation
    for(int i=0; i<N; i++){
        double s=0, c=0, p=0, c_disc=0, p_disc=0;

        for(int j=0; j<L; j++){
            double Z=rnd.Gauss(0,1);
            s = S(S0,r,sigma,Z,T);
            c += C(r,T,s,K);
            p += P(r,T,s,K);

            double s_disc=S0;
            for(int k=0; k<100;k++){
                double Z_disc = rnd.Gauss(0,1);
                s_disc = S(s_disc, r, sigma, Z_disc, T/100.);
            }
            c_disc += C(r,T,s_disc,K);
            p_disc += P(r,T,s_disc,K);
        }

        avg1_call.push_back(c/L);
        avg1_put.push_back(p/L);
        avg2_call.push_back(avg1_call[i] * avg1_call[i]);
        avg2_put.push_back(avg1_put[i] * avg1_put[i]);

        avg1_call_disc.push_back(c_disc/L);
        avg1_put_disc.push_back(p_disc/L);
        avg2_call_disc.push_back(avg1_call_disc[i] * avg1_call_disc[i]);
        avg2_put_disc.push_back(avg1_put_disc[i] * avg1_put_disc[i]);
    }

    for(int i=0; i<N; i++){
        // Initialize to zero (no first step calculation)
        sum_call.push_back(0);
        sum2_call.push_back(0);
        sum_put.push_back(0);
        sum2_put.push_back(0);

        sum_call_disc.push_back(0);
        sum2_call_disc.push_back(0);
        sum_put_disc.push_back(0);
        sum2_put_disc.push_back(0);

        for(int j=0; j<i+1; j++){
            sum_call[i] += avg1_call[j];
            sum_put[i] += avg1_put[j];
            sum2_call[i] += avg2_call[j];
            sum2_put[i] += avg2_put[j];

            sum_call_disc[i] += avg1_call_disc[j];
            sum_put_disc[i] += avg1_put_disc[j];
            sum2_call_disc[i] += avg2_call_disc[j];
            sum2_put_disc[i] += avg2_put_disc[j];
        }
        sum_call[i] /= (i+1);
        sum2_call[i] /= (i+1);
        sum_put[i] /= (i+1);
        sum2_put[i] /= (i+1);

        sum_call_disc[i] /= (i+1);
        sum2_call_disc[i] /= (i+1);
        sum_put_disc[i] /= (i+1);
        sum2_put_disc[i] /= (i+1);

        err_call.push_back(error(sum_call,sum2_call,i));
        err_put.push_back(error(sum_put,sum2_put,i));

        err_call_disc.push_back(error(sum_call_disc,sum2_call_disc,i));
        err_put_disc.push_back(error(sum_put_disc,sum2_put_disc,i));
    }

    // Output File
    ofstream out_call;
    ofstream out_put;
    ofstream out_call_disc;
    ofstream out_put_disc;

    out_call.open("out_call.data");
    out_put.open("out_put.data");
    out_call_disc.open("out_call_disc.data");
    out_put_disc.open("out_put_disc.data");

    for(int i=0; i<N; i++){
        out_call << (i+1)*L << "\t" << sum_call[i] << "\t" << err_call[i] << endl;
        out_put << (i+1)*L << "\t" << sum_put[i] << "\t" << err_put[i] << endl;

        out_call_disc << (i+1)*L << "\t" << sum_call_disc[i] << "\t" << err_call_disc[i] << endl;
        out_put_disc << (i+1)*L << "\t" << sum_put_disc[i] << "\t" << err_put_disc[i] << endl;
    }
    out_call.close();
    out_put.close();
    out_call_disc.close();
    out_put_disc.close();

    rnd.SaveSeed();
    return 0;
}