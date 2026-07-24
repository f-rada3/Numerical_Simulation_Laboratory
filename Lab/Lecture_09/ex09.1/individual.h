#ifndef __INDIVIDUAL_H__
#define __INDIVIDUAL_H__

#include <armadillo>
#include "random.h"

using namespace std;
using namespace arma;


class Individual{
    private:
        uvec _path;             // vettore di interi di armadillo
        double _fitness;
        int _num_of_cities;
        mat _distance_matrix;   // matrici armadillo

    public:
        Individual(int num_of_cities, const mat &distance_matrix);

        // Main functions
        void Initialize(Random & rnd);
        void Check();
        void Swap(int i, int j);
        void Fitness();

        // Mutations
        void PairPermutation(Random & rnd);
        void AdiacentPermutation(Random & rnd);
        void Shift(Random & rnd);
        void Inversion(Random & rnd);

        // Get & Set functions
        double GetFitness() const { return _fitness; }
        uvec GetPath() const { return _path; }
        void SetPath( const uvec &path );
}; 

#endif