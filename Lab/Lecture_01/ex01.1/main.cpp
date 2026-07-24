#include <fstream>
#include <iostream>
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

double chi2(const vector<int>& observed, double expected) {
    // This function calculates the chi-squared statistic
    double chi2 = 0.0;
    for (int i=0; i<observed.size(); i++){
        chi2 += pow((observed[i] - expected),2)/expected;
    }
    return chi2;
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
    int M=10000; // throws
    int N=100;  // blocks
    int L=M/N;  // number of throws in each block

    vector <double> avg1, avg2, var1, var2;
    vector <double> sum, sum2, err;
    vector <double> sum_var, sum_var2, err_var;

    // Integral evaluation
    for(int i=0; i<N; i++){
        double sum=0;
        double sum_var=0;
        for(int j=0; j<L; j++){
            sum += rnd.Rannyu();
            sum_var += pow(rnd.Rannyu()-0.5,2);
        }
        avg1.push_back(sum/L);
        avg2.push_back(avg1[i]*avg1[i]);
        var1.push_back(sum_var/L);
        var2.push_back(var1[i]*var1[i]);
    }  

    for(int i=0; i<N; i++){
        sum.push_back(0);
        sum2.push_back(0);
        sum_var.push_back(0);
        sum_var2.push_back(0);

        for(int j=0; j<i+1; j++){
            sum[i] += avg1[j];
            sum2[i] += avg2[j];
            sum_var[i] += var1[j];
            sum_var2[i] += var2[j];
        }        
        sum[i] /= (i+1);
        sum2[i] /= (i+1);
        sum_var[i] /= (i+1);
        sum_var2[i] /= (i+1);
        err.push_back(error(sum, sum2, i));
        err_var.push_back(error(sum_var, sum_var2, i));
    }

    // Output file
    ofstream out_avg;
    ofstream out_var;
    out_avg.open("out_avg.data");
    out_var.open("out_var.data");
    for(int i=0; i<N; i++){
        out_avg << (i+1)*L <<"\t" << sum[i] << "\t" << err[i] << endl;
        out_var << (i+1)*L <<"\t" << sum_var[i] << "\t" << err_var[i] << endl;
    }
    out_avg.close();
    out_var.close();

    // Results on terminal
    cout << "\n\nBlocks number = " << N << endl;
    cout << "Number of throws for each block = " << L << endl;
    cout << "\nBest integral value = " << sum[N-1] << " +-" << err[N-1] << "\n\n" << endl;

    // Chi2 test
    int n = 10000;          // throws
    int m=100;              // number of intervals in [0,1)
    // the number of expected events in each sub-interval, 
    // according to a normal distribution is np = n x 1/m
    int exp_events= n/m;    // expected events in each sub-interval 
    vector <int> test_cases{100,1000,10000};

    for (int num_test : test_cases){
        vector <double> chi2_values;
        string out_chi_file = "out_chi_" + to_string(num_test) + ".data";
        
        for(int i=0; i< num_test; i++){
            vector <int> counts(m,0.);
            for (int j=0; j<n; j++){
                int index = static_cast<int>(rnd.Rannyu()*m);    // interval index which contains rnd.Rannyu()
                counts[index]++;
            }
            chi2_values.push_back(chi2(counts, exp_events));
        }

        ofstream out_chi(out_chi_file);
        for(double chi_val : chi2_values){
            out_chi << chi_val << endl;
        }
        out_chi.close();
    }

    rnd.SaveSeed();
    
    return 0;
}