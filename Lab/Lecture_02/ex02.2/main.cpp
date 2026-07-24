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

void DiscreteRW(int M, int N, int a, int Nsteps){

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
    int L=M/N;
    vector <double> global_avg(Nsteps, 0.), global_avg2(Nsteps, 0.), prog_avg(Nsteps, 0.), prog_avg2(Nsteps, 0.), prog_err(Nsteps, 0.);
    ofstream traj("discrete_traj.data");

    // Implementation
    for(int i=0; i<N; i++){
        vector <double> block_avg(Nsteps,0.);

        for(int j=0; j<L; j++){
            double x=0.,y=0.,z=0.;

            for(int step=0; step<Nsteps; step++){
                // Save trajectory
                if(j==0 && i==0) 
                traj << step+1 << "\t" << x << "\t" << y << "\t" << z << endl;

                int dir = int(rnd.Rannyu(0,3));             // Chooses a direction (x,y,z)
                int sign = (rnd.Rannyu(-1,1) < 0) ? : -1;   // Chooses positive or negative direction 

                if(dir==0) x += a*sign;
                else if(dir==1) y += a*sign;
                else z += a*sign;

                block_avg[step] += x*x + y*y + z*z;
            }
        }
        for(int step=0; step<Nsteps; step++){
            block_avg[step] /= L;
            global_avg[step] += block_avg[step];
            global_avg2[step] += block_avg[step] * block_avg[step];
        }
    }
    traj.close();

    // Save Averages
    ofstream out_avg("out_avg_discrete.data");
    for(int step=0; step<Nsteps; step++){
        prog_avg[step] = global_avg[step]/N;
        prog_avg2[step] = global_avg2[step]/N;
        prog_err[step] = error(prog_avg, prog_avg2, step);

        out_avg << step+1 << "\t" << sqrt(prog_avg[step]) << "\t" << prog_err[step] / (2*sqrt(prog_avg[step])) << endl; // Correct error on sqrt(r^2)!
    }

    rnd.SaveSeed(); 
}

void ContinuosRW(int M, int N, int a, int Nsteps){

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
    int L=M/N;
    vector <double> global_avg(Nsteps, 0.), global_avg2(Nsteps, 0.), prog_avg(Nsteps, 0.), prog_avg2(Nsteps, 0.), prog_err(Nsteps, 0.);
    ofstream traj("continuous_traj.data");

    // Implementation
    for(int i=0; i<N; i++){
        vector <double> block_avg(Nsteps,0.);

        for(int j=0; j<L; j++){
            double x=0.,y=0.,z=0.;

            for(int step=0; step<Nsteps; step++){
                double costheta = rnd.Rannyu(-1,1);
                double theta = acos(costheta);
                double phi = rnd.Rannyu(0, 2*M_PI);

                x += a*sin(theta)*cos(phi);
                y += a*sin(theta)*sin(phi);
                z += a*cos(theta);

                block_avg[step] += x*x + y*y + z*z;

                // Save trajectory
                if(j==0 && i==0) 
                traj << step+1 << "\t" << x << "\t" << y << "\t" << z << endl;
            }
        }
        for(int step=0; step<Nsteps; step++){
            block_avg[step] /= L;
            global_avg[step] += block_avg[step];
            global_avg2[step] += block_avg[step] * block_avg[step];
        }
    }
    traj.close();

    // Save Averages
    ofstream out_avg("out_avg_continuous.data");
    for(int step=0; step<Nsteps; step++){
        prog_avg[step] = global_avg[step]/N;
        prog_avg2[step] = global_avg2[step]/N;
        prog_err[step] = error(prog_avg, prog_avg2, step);

        out_avg << step+1 << "\t" << sqrt(prog_avg[step]) << "\t" << prog_err[step] / (2*sqrt(prog_avg[step])) << endl; // Correct error on sqrt(r^2)!
    }

    rnd.SaveSeed(); 
}

int main(int argc, char *argv[]){

    // Parameters
    int M = 10000;      // Number of Total Random Walks
    int N = 100;        // Number of Blocks
    int a=1;            // Step Size
    int Nsteps=100;     // Maximum Number of Steps for Each Random Walk

    DiscreteRW(M,N,a,Nsteps);
    ContinuosRW(M,N,a,Nsteps);

    return 0;
}