#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdlib>
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
    int seed[4];
    int p1, p2;
    ifstream Primes("Primes");
    if (Primes.is_open()){
        Primes >> p1 >> p2 ;
    } else cerr << "PROBLEM: Unable to open Primes" << endl;
    Primes.close();

    ifstream input("seed.in");
    string property;
    if (input.is_open()){
        while ( !input.eof() ){
            input >> property;
            if( property == "RANDOMSEED" ){
                input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
                rnd.SetRandom(seed,p1,p2);
            }
        }
        input.close();
    } else cerr << "PROBLEM: Unable to open seed.in" << endl;

    ofstream output_SA("SimAnn.dat");                   // File containing every step of SA 
    if (!output_SA){
        cerr << "PROBLEM: Unable to open SimmAnn.dat" << endl;
    }
    ofstream output_SA_accepted("SimAnn_accepted.dat"); // FIle containig just accepted steps
    if (!output_SA){
        cerr << "PROBLEM: Unable to open SimmAnn_accepted.dat" << endl;
    }

    // Parameters
    int n_T=1000;
    int n_cycles=20;
    double step=0.05, sigma0=1., mu0=1.;

    double energy_prev = energy_Metropolis(rnd, sigma0, mu0);
    double energy_inst = 0.; 

    // Simulated Annealing
    double best_sigma=0., best_mu=0.;
    for(int i=0; i<n_T; i++){
        cout << "STEP " << i+1 << "\t" << "step size = " << step << endl;

        int n_accept=0;

        double T=temperature(i);
        for(int j=0; j<n_cycles; j++){
            double sigma = sigma0 + rnd.Rannyu(-step,step);
            double mu = mu0 + rnd.Rannyu(-step,step);
            energy_inst = energy_Metropolis(rnd, sigma, mu);

            if(rnd.Rannyu() < acceptance_ratio_SA(energy_prev,energy_inst,T)){
                sigma0 = sigma;
                mu0 = mu;
                energy_prev = energy_inst;
                n_accept++;
                cout << "sigma = " << sigma << "\tmu = " << mu << "\tenergy = " << energy_prev << endl;
                output_SA_accepted << T << "\t" << mu0 << "\t" << sigma0 << "\t" << energy_prev << endl;
            }
            output_SA << T << "\t" << mu0 << "\t" << sigma0 << "\t" << energy_prev << endl;
            if (i==1000 && j==n_cycles-1){
                best_mu = mu0;
                best_sigma = sigma0;
            }
        }
        cout << "T = " << T << "\tAcceptance = " << double((n_accept)/n_cycles*100) << "%" << endl;

        if (double(n_accept)/n_cycles<0.40 && step > 0.01) step *= 0.9;
        else if (double(n_accept)/n_cycles>0.60 && step < 2.) step *= 1.1;
    }

    // Show and save results
    cout << "========================================" << endl;
    cout << "Simulated Annealing completed." << endl;
    cout << "Best parameters found:" << endl;
    cout << "Best mu: " << setprecision(10) << best_mu << endl;
    cout << "Best sigma: " << setprecision(10) << best_sigma << endl;
    cout << "========================================" << endl;

    ofstream output_best_params("best_params.dat");
    if (!output_best_params) {
        cerr << "PROBLEM: Unable to open best_params.dat" << endl;
    } else {
        output_best_params << "# Best parameters found by Simulated Annealing" << endl;
        output_best_params << fixed << setprecision(10); // Imposta precisione per il file
        output_best_params << "Best mu: " << best_mu << endl;
        output_best_params << "Best sigma: " << best_sigma << endl;
        output_best_params.close();
    }

    output_SA.close();
    output_SA_accepted.close();
    rnd.SaveSeed();
    return 0;
}