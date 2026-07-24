#include <iostream>
#include <fstream>
#include <vector>
#include <numeric> 
#include <algorithm>
#include <random> 
#include <mpi.h>
#include "population.h"
#include "random.h"

using namespace std;

int main(int argc, char * argv[]){

    // inizializzazione MPI
    int rank, size;
    MPI_Init(&argc, &argv);              
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Get the ID of the current process (0, 1, 2...)
    MPI_Comm_size(MPI_COMM_WORLD, &size); // Get total number of processes running
    MPI_Barrier(MPI_COMM_WORLD); 
    double start = MPI_Wtime();

    Random rnd;
    int seed[4], p1, p2;
    ifstream Primes("Primes");
    if (Primes.is_open()) {
        for(int i = 0; i <= rank; i++) {
            if(!(Primes >> p1 >> p2)){
                cerr << "Error reading Primes for rank " << rank << " at step " << i << endl;
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }
    } else { 
        if(rank == 0) cerr << "Error: Primes file missing" << endl; 
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    Primes.close();

    ifstream input("seed.in");
    if (input.is_open()) {
        string property;
        if(!(input >> property >> seed[0] >> seed[1] >> seed[2] >> seed[3])){
            cerr << "Error reading seed.in for rank " << rank << endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        rnd.SetRandom(seed, p1, p2);

        cerr << "Rank " << rank 
            << " p1=" << p1 
            << " p2=" << p2
            << " seed=(" << seed[0] << "," << seed[1] << "," << seed[2] << "," << seed[3] << ")"
            << endl;

        for(int t = 0; t < 5; t++){
            double r = rnd.Rannyu();
            cerr << "Rank " << rank << " rnd[" << t << "] = " << r << endl;
        }

        input.close();
    } else {
        cerr << "Error: seed.in missing for rank " << rank << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // parametri
    int num_of_cities = 110;
    int npop = 700;
    int ngen = 4000;
    int nmig = 50; // numero di generazioni tra migrazioni
    double mut_prob = 0.1;
    double p_exp = 3.0;
    bool migration = false;

    Population population;
    population.InitializeByFile(rnd, "cap_prov_ita.dat", num_of_cities, npop, ngen, mut_prob, p_exp, rank);

    // loop principale per l'evoluzione
    for(int gen=0; gen<ngen; gen++){
        population.NewGeneration();

        // migrazione ogni nmig generazioni
        if(migration && gen>0 && (gen % nmig==0)){

            vector<int> destinations(size);
            iota(destinations.begin(), destinations.end(), 0); // Riempimento con 0, 1, ..., size-1

            if(rank==0){
                // Rank 0 sceglie casualmente le destinazioni per la migrazione
                shuffle(destinations.begin(), destinations.end(), default_random_engine(rnd.Rannyu()));
            }
            // Broadcast delle destinazioni a tutti i processi
            MPI_Bcast(destinations.data(), size, MPI_INT, 0, MPI_COMM_WORLD);

            arma::uvec best_path = population.GetBestPath();
            arma::uword* send_buffer = best_path.memptr(); // Puntatore al buffer dei dati da inviare
            arma::uword* recv_buffer = new arma::uword[num_of_cities]; // Buffer per ricevere i dati

            // determinazione destinazione e provenienza
            int target_rank = destinations[rank];
            int source_rank = -1;
            for(int i=0; i<size; i++){
                if(destinations[i] == rank){
                    source_rank = i;
                    break;
                }
            }

            // invio e ricezione del percorso migliore
            MPI_Status status;
            MPI_Sendrecv(send_buffer, num_of_cities, MPI_UNSIGNED_LONG, target_rank, 1,
                         recv_buffer, num_of_cities, MPI_UNSIGNED_LONG, source_rank, 1,
                         MPI_COMM_WORLD, &status);

            // conversione del buffer ricevuto in arma::uvec
            arma::uvec received_path(num_of_cities);
            for(int i=0; i<num_of_cities; i++){
                received_path(i) = recv_buffer[i];
            }
            population.EvaluateImmigrant(received_path);
            delete[] recv_buffer; // Deallocazione del buffer di ricezione
        }
    }

    // conclusione
    double end = MPI_Wtime();
    double elapsed = end - start;
    double max_elapsed;
    MPI_Reduce(&elapsed, &max_elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    // risultati finali
    double local_best = population.BestFitness();
    double global_best;
    MPI_Reduce(&local_best, &global_best, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);

    if(rank==0){
        cout << "\nBest fitness found: " << global_best << endl;
        cout << "Total elapsed time: " << max_elapsed << " seconds" << endl;
    }

    population.Closure(rank);
    MPI_Finalize();
    return 0;
}