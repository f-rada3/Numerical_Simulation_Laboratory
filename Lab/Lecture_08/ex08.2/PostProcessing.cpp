#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include "random.h"

using namespace std;

// Functions

double error(const vector<double> &sum_prog, const vector<double> &sum_prog2, int n) {
    if (n == 0) return 0;
    return sqrt((sum_prog2[n] - pow(sum_prog[n], 2)) / n);
}

double acceptance_ratio(double prob_old, double prob_new) {
   return min(1., prob_new / prob_old);
}

double acceptance_ratio_SA(double energy_old, double energy_new, double T) {
   return min(1., exp(-1./T*(energy_new-energy_old)));
}

double wave_function(double x, double sigma, double mu){
    return (exp(-pow(x - mu, 2) / (2*pow(sigma, 2))) + 
            exp(-pow(x + mu, 2) / (2*pow(sigma, 2))));
}

double sqrd_mod_wave_function(double x, double sigma, double mu){
    return exp(-pow(x - mu, 2) / (pow(sigma, 2))) + 
            exp(-pow(x + mu, 2) / (pow(sigma, 2))) + 
            2 * exp( - (pow(x + mu, 2) + pow(x - mu, 2))/ (2 * pow(sigma, 2)));
}

double kinetic_energy(double x, double sigma, double mu){
    // prefactor
    double f0 = -1. / (sigma * sigma);
    double f1 = pow((x + mu) / sigma, 2);
    double f2 = pow((x - mu) / sigma, 2);
    double f3 = exp(-f1 / 2.);
    double f4 = exp(-f2 / 2.);
    double dderiv_wave_function = f0 * ((1. - f1) * f3 + (1. - f2) * f4);

    // K = -0.5 * dderiv_wave_function / wave_function
    return -0.5 *  dderiv_wave_function / (wave_function(x,sigma,mu));
}

double potential_energy(double x){
    return pow(x, 4) - 2.5 * x * x;
}

double energy(double x, double sigma, double mu){
    return kinetic_energy(x,sigma,mu) + potential_energy(x);
}

double temperature(int step){
    return 3. / (1. + 4.*step);
}

double energy_Metropolis(Random &rnd, double sigma, double mu){
    int M=100000;
    double avg_energy=0.;
    double x=0.;
    double step=2.;

    for(int i=0; i<M; i++){
        double x1=rnd.Rannyu(-step,step);
        double alpha = min(1., sqrd_mod_wave_function(x+x1,sigma,mu)/sqrd_mod_wave_function(x,sigma,mu));

        if (rnd.Rannyu() <= alpha){
            x = x+x1;
        }
        avg_energy += energy(x,sigma,mu);
    }
    return avg_energy/M;
}

void data_blocking_sqrd_mod_wave_function(int M, int N, Random &rnd, double sigma, double mu){
    double x_start = 0;
    double step=3.;
    int L=M/N;

    vector <double> avg1, avg2;
    vector <double> sum1, sum2, err;

    ofstream steps_data("steps.dat");
    ofstream avg_energy_data("avg_energy.dat");

    for(int i=0; i<N; i++){
        double sum=0;
        for(int j=0; j<L; j++){
            double x=rnd.Rannyu(-step,step);
            double acceptance = acceptance_ratio(sqrd_mod_wave_function(x_start,sigma,mu),
                                                sqrd_mod_wave_function(x_start + x,sigma,mu));
            
            if (rnd.Rannyu() <= acceptance){
                x_start += x;
            }
            steps_data << x_start << endl;
            sum += energy(x_start,sigma,mu);
        }
        avg1.push_back(sum/L);
        avg2.push_back(avg1[i]*avg1[i]);
    }
    // Data Blocking
    for (int i=0; i<N; i++){
        sum1.push_back(0);
        sum2.push_back(0);
        for (int j=0; j<i+1; j++){
            sum1[i] += avg1[j];
            sum2[i] += avg2[j];
        }
        sum1[i] /= (i+1);
        sum2[i] /= (i+1);
        err.push_back(error(sum1,sum2,i));
        avg_energy_data << (i+1)*L << "\t" << sum1[i] << "\t" << err[i] << endl;
    }
}


int main(int argc, char *argv[]){
    
    // Files
    Random rnd;
    int seed[4], p1, p2;
    ifstream Primes("Primes");
    if (Primes.is_open()) Primes >> p1 >> p2;
    else cerr << "PROBLEM: Unable to open Primes" << endl;
    Primes.close();

    ifstream input("seed.in");
    string property;
    if (input.is_open()) {
        while (!input.eof()) {
            input >> property;
            if (property == "RANDOMSEED") {
                input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
                rnd.SetRandom(seed, p1, p2);
            }
        }
        input.close();
    } else {
        cerr << "PROBLEM: Unable to open seed.in" << endl;
        return 1;
    }

    ifstream in_params("SimAnn_accepted.dat");
    if(!in_params.is_open()){
        cerr << "PROBLEM: Unable to open SimAnn_accepted.dat" << endl;
        return 1;
    }

    ofstream out_params_and_energy("energy_vs_SA_steps.dat");
    if(!out_params_and_energy.is_open()){
        cerr << "PROBLEM: Unable to open energy_vs_SA_steps.dat" << endl;
        return 1;
    }
    out_params_and_energy << "# mu\t\tsigma\t\tenergy\t\terror\n";

    string line;
    int nline = 0;
    while(getline(in_params,line)){
        istringstream iss(line);
        double temp, mu, sigma;
        if(!(iss >> temp >> mu >> sigma)) continue;

        // Parameters
        int M=10000, N=100, L=M/N;
        //int accept = 0;
        double step = 3.;
        double x_start = 0.;

        vector <double> avg1, avg2;
        vector <double> sum1, sum2, err;

        for(int i=0; i<N; i++){
            double sum=0;
            for (int j=0; j<L; j++){
                double x = rnd.Rannyu(-step,step);
                double acceptance = acceptance_ratio(sqrd_mod_wave_function(x_start,sigma,mu), 
                                                    sqrd_mod_wave_function(x_start + x,sigma,mu));

                if(rnd.Rannyu() <= acceptance){
                    x_start += x;
                    //accept++;
                }
                sum += energy(x_start,sigma,mu);
            }
            avg1.push_back(sum/L);
            avg2.push_back(avg1[i] * avg1[i]);
        }
        for(int i=0; i<N; i++){
            double s1=0, s2=0;
            for (int j=0; j<i+1; j++){
                s1 += avg1[j];
                s2 += avg2[j];
            }
            sum1.push_back(s1/(i+1));
            sum2.push_back(s2/(i+1));
            err.push_back(error(sum1, sum2, i));
        }
        nline++;
        double ENERGY = sum1[N-1];
        double ERROR = err[N-1];
        out_params_and_energy << nline << "\t" << mu << "\t" << sigma << "\t" 
                        << ENERGY << "\t" << ERROR << endl; 
    }

    // Data blocking for |Psi|^2
    double sigma_best = 0.618567;
    double mu_best = 0.803221;
    // nota: potrei implementare pezzo di codice per leggerlo in automatico
    data_blocking_sqrd_mod_wave_function(1000000, 100, rnd, sigma_best, mu_best);

    in_params.close();
    out_params_and_energy.close();
    rnd.SaveSeed();
    return 0;
}