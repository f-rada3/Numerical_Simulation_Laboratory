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

double f(double x){
    return (M_PI/2.) * cos((M_PI * x)/2.);
}

double p(double x){
    return -2*x + 2;
}

double g(double x){
    return (20./3.)*cos((M_PI*x)/2.) / (8-pow(M_PI*x , 2));
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
    int M=100000;       // throws
    int N=100;          // blocks
    int L=M/N;          // number of throws for each block
    
    vector <double> avg1, avg2, var1, var2, avg1_imp_samp, avg2_imp_samp;
    vector <double> sum, sum2, err;
    vector <double> sum_var, sum_var2, err_var;
    vector <double> sum_imp_samp, sum_imp_samp2, err_imp_samp;

    // Integral evaluation
    for(int i=0; i<N; i++){
        double temp=0.;
        double temp_var=0.;
        double temp_imp_samp=0.;

        for(int j=0; j<L; j++){
            double x_uniform = rnd.Rannyu();
            temp += f(x_uniform);
            temp_var += pow(x_uniform-1, 2);

            double x_sampl = rnd.Linear();
            temp_imp_samp += f(x_sampl)/p(x_sampl);
        }
        avg1.push_back(temp/L);
        avg2.push_back(avg1[i]*avg1[i]);  

        var1.push_back(temp_var/L);
        var2.push_back(var1[i]*var1[i]);

        avg1_imp_samp.push_back(temp_imp_samp/L);
        avg2_imp_samp.push_back(avg1_imp_samp[i]*avg1_imp_samp[i]);
    }

    for(int i=0; i<N; i++){
        sum.push_back(0);
        sum2.push_back(0);

        sum_var.push_back(0);
        sum_var2.push_back(0);

        sum_imp_samp.push_back(0);
        sum_imp_samp2.push_back(0);

        for(int j=0; j<i+1; j++){
            sum[i] += avg1[j];
            sum2[i] += avg2[j];
            sum_var[i] += var1[j];
            sum_var2[i] += var2[j];
            sum_imp_samp[i] += avg1_imp_samp[j];
            sum_imp_samp2[i] += avg2_imp_samp[j];
        }
        sum[i] /= (i+1);
        sum2[i] /= (i+1);
        sum_var[i] /= (i+1);
        sum_var2[i] /= (i+1);
        sum_imp_samp[i] /= (i+1);
        sum_imp_samp2[i] /= (i+1);

        err.push_back(error(sum,sum2,i));
        err_var.push_back(error(sum_var,sum_var2,i));
        err_imp_samp.push_back(error(sum_imp_samp,sum_imp_samp2,i));
    }

    // Output File
    ofstream out_avg;
    ofstream out_var;
    ofstream out_imp_samp;

    out_avg.open("out_avg.data");
    out_var.open("out_var.data");
    out_imp_samp.open("out_imp_samp.data");

    for(int i=0; i<N; i++){
        out_avg << (i+1)*L << "\t" << sum[i] << "\t" << err[i] << endl;
        out_var << (i+1)*L << "\t" << sum_var[i] << "\t" << err_var[i] << endl;
        out_imp_samp << (i+1)*L << "\t" << sum_imp_samp[i] << "\t" << err_imp_samp[i] << endl;
    }
    out_avg.close();
    out_var.close();
    out_imp_samp.close();

    // Results on Terminal
    cout << "\n\nNumber of blocks = " << N << endl;
    cout << "Number of throws per block = " << L << endl;
    cout << "Number of total throws = " << M << endl;
    cout << "Best integral value = " << sum[N-1] << " +-" << err[N-1] << endl;
    cout << "Best integral value (importance sampling) = " << sum_imp_samp[N-1] << " +-" << err_imp_samp[N-1] << "\n\n" << endl;

    rnd.SaveSeed();

    return 0;
}