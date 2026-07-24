#ifndef __POPULATION_H__
#define __POPULATION_H__

#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <cmath>

#include "individual.h"

using namespace std;
using namespace arma;

class Population {
    private:
        int _num_of_cities;
        int _npop;
        int _ngen;
        int _gen_index;
        double _mut_prob;
        double _p_exponent;
        string _type;

        mat _pos;
        mat _distance_matrix;
        vector<Individual> _individuals;
        Random* _rnd;

        ofstream _loss;

        void InitializeCircle();
        void InitializeSquare();
        void CreateDistanceMatrix();
        void FitnessSorting();
        Individual SelectRank();
        
        void Crossover(Individual &parent1, Individual &parent2, Individual &child1, Individual &child2);
        double GetBestHalfAvg();

    public:
        Population();
        ~Population();

        void Initialize(Random & rnd, const string & type, int num_of_cities, int npop, int ngen, double mut_prob, double p_exp);
        void NewGeneration();
        int GetNumGen();
        void Closure();

        // metodi parallelizzazione
        double BestFitness();
        void InitializeByFile(Random & rnd, const string & filename, int num_of_cities, int npop, int ngen, double mut_prob, double p_exp, int rank);
        arma::uvec GetBestPath() { return _individuals[0].GetPath();}
        void EvaluateImmigrant(const arma::uvec & path);
        void Closure(int rank);
};


#endif