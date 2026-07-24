#include <iostream>
#include <fstream>

#include "population.h"
#include "random.h"

using namespace std;

int main(){

    Random rnd;
    int seed[4], p1, p2;
    ifstream Primes("Primes");
    if (Primes.is_open()) Primes >> p1 >> p2;
    else { cerr << "Error: Primes file missing" << endl; return 1; }
    Primes.close();

    ifstream input("seed.in");
    if (input.is_open()) {
        input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
        rnd.SetRandom(seed, p1, p2);
        input.close();
    } else { cerr << "Error: seed.in missing" << endl; return 1; }

    // Parameters
    int num_of_cities=34;
    int npop=200;
    int ngen=600;
    double mut_prob=0.2;
    double p_exponent=3.;

    // Circle
    Population circle;
    circle.Initialize(rnd, "circle", num_of_cities, npop, ngen, mut_prob, p_exponent);

    for(int i=0; i<circle.GetNumGen(); i++){
        circle.NewGeneration();
    }
    circle.Closure();

    // Square
    Population square;
    square.Initialize(rnd, "square", num_of_cities, npop, ngen, mut_prob, p_exponent);

    for(int i=0; i<square.GetNumGen(); i++){
        square.NewGeneration();
    }
    square.Closure();

    cout << "Simulations for circle and square are completed.\n" << endl;
    rnd.SaveSeed();

    return 0;
}
