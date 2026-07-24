#include <armadillo>
#include "population.h"

using namespace std;
using namespace arma;

Population::Population() {}
Population::~Population() {
    if(_loss.is_open()) _loss.close();
}

void Population::Initialize(Random & rnd, const string & type, int num_of_cities, int npop, int ngen, double mut_prob, double p_exp) {
    _rnd = &rnd;
    _type = type;
    _num_of_cities = num_of_cities;
    _npop = npop;
    _ngen = ngen;
    _mut_prob = mut_prob;
    _p_exponent = p_exp;
    _gen_index = 0;
    _pos.set_size(_num_of_cities , 2);

    if(_type=="circle") InitializeCircle();
    else if(_type=="square") InitializeSquare();

    // calcolo delle distanze
    CreateDistanceMatrix();

    // loop per popolazione iniziale di npop individui
    for(int i=0; i<_npop; i++){
        Individual ind(_num_of_cities, _distance_matrix);
        // inizializzazione random
        ind.Initialize(*_rnd);
        // aggiunta del nuovo individuo alla popolazione
        _individuals.push_back(ind);
    }
    FitnessSorting();

    // creazione file output
    _loss.open(_type + "_losses.dat");
    cout << "\nOptimization: " << _type << endl; 
}

void Population::InitializeCircle(){
    for(int i=0; i<_num_of_cities; ++i){        //pre increment
        double alpha = _rnd->Rannyu(0, 2*M_PI); //angolo random
        _pos(i, 0) = cos(alpha);                //coor x
        _pos(i, 1) = sin(alpha);                //coor y
    }
}

void Population::InitializeSquare(){
    for(int i=0; i<_num_of_cities; ++i){
        //distribuisco unif tra -1,1 sia coor x che y
        _pos(i, 0) = _rnd->Rannyu(-1.0, 1.0);   
        _pos(i, 1) = _rnd->Rannyu(-1.0, 1.0);   
    }
}

void Population::CreateDistanceMatrix(){
    // Modifica _distance_matrix
    // Inizializzo matrice per conterere dist. num_of_cities x num_of_cities
    _distance_matrix.set_size(_num_of_cities, _num_of_cities);
    // loop su riga corrispondente a prima città
    for(int i=0; i<_num_of_cities; ++i){
        // loop su indice colonna per la seconda città
        for(int j=i+1; j<_num_of_cities; ++j){
            // devo partire dal successivo di i per non fare doppi calcoli
            _distance_matrix(i,j) = norm(_pos.row(i) - _pos.row(j));
            // simmetria
            _distance_matrix(j,i) = _distance_matrix(i,j);
        }
    }
}

// Riordino popolazione da migliore a peggiore in base a fitness
void Population::FitnessSorting(){
    sort(_individuals.begin(), _individuals.end(), [](const Individual &low, const Individual &high){
        return low.GetFitness() < high.GetFitness();
    }
    );
}

// Selezione di un individuo per la riproduzione
Individual Population::SelectRank(){
    // calcolo indice usando rnd number elevato a _p_exponent
    int index = static_cast<int>(_npop * pow(_rnd->Rannyu(), _p_exponent));
    // return individual selezionato tra popolazione riordinata
    return _individuals[index];
}

void Population::Crossover(Individual &parent1, Individual &parent2, Individual &child1, Individual &child2){
    uvec x1 = parent1.GetPath();
    uvec x2 = parent2.GetPath();
    
    // scelta random di pto di crossover (no primo o ultimo elemento)
    int division = static_cast<int>(_rnd->Rannyu(1, _num_of_cities-1));
    // indici per individuare la parte di genoma dei figli oltre la divisione
    int i1=division+1, i2=division+1;   
    
    // vettori per trasporto informazione genetica da parents a children
    uvec sec1(_num_of_cities);
    uvec sec2(_num_of_cities);

    // copia della sequenza dei parents fino a division
    for(int i=0; i<=division; i++){
        sec1(i) = x1(i);
        sec2(i) = x2(i);
    }

    // Iterazione su tutta la sequenza
    for(int i=0; i<_num_of_cities; i++){
        // Flags per tracciare se una città esiste già nel genoma del figlio
        bool flag1=false, flag2=false;  // inizializzo false
        // loop entro 'division' per controllare la parte già ereditata 
        for(int j=0; j<=division; j++){
            if(x2(i) == sec1(j)){
                // flag1=true se la città da parent2 è già in child1
                flag1 = true;
            }
            if(x1(i) == sec2(j)){
                // flag2=true se la città da parent 1 è già in child2
                flag2=true;
            }
        }
        // se la città non è in child1, aggiungo e incremento l'indice di child1
        if(!flag1) sec1(i1++) = x2(i);
        // se la città non è in child2, aggiungo e incremento l'indice di child2
        if(!flag2) sec2(i2++) = x1(i);
    } 
    child1.SetPath(sec1); // assegnazione nuovo genoma a child1
    child2.SetPath(sec2); // assegnazione nuovo genoma a child2
}

void Population::NewGeneration(){
    vector<Individual> future_generation;   // temp vec
    future_generation.reserve(_npop);       // pre-allocazione memoria per nuova popolazione

    // loop per far raggiungere le dimensioni attese per la popolazione
    while(future_generation.size() < static_cast<size_t>(_npop)){
        // selezione parents usando operatore di selezione rango (biased)
        Individual parent1 = SelectRank();
        Individual parent2 = SelectRank();

        // inizializzazione figli 
        Individual child1(_num_of_cities, _distance_matrix);
        Individual child2(_num_of_cities, _distance_matrix);

        // eredità
        if(_rnd->Rannyu()<0.7){
            // con prob 70% avviene crossover
            Crossover(parent1, parent2, child1, child2);
        } else {
            // altrimenti viene copiato il genoma dei genitori
            child1.SetPath(parent1.GetPath());
            child2.SetPath(parent2.GetPath());
        }

        if(_rnd->Rannyu() < _mut_prob) child1.PairPermutation(*_rnd);
        if(_rnd->Rannyu() < _mut_prob) child1.Shift(*_rnd);
        if(_rnd->Rannyu() < _mut_prob) child1.AdiacentPermutation(*_rnd);
        if(_rnd->Rannyu() < _mut_prob) child1.Inversion(*_rnd);

        if(_rnd->Rannyu() < _mut_prob) child2.PairPermutation(*_rnd);
        if(_rnd->Rannyu() < _mut_prob) child2.Shift(*_rnd);
        if(_rnd->Rannyu() < _mut_prob) child2.AdiacentPermutation(*_rnd);
        if(_rnd->Rannyu() < _mut_prob) child2.Inversion(*_rnd);

        child1.Fitness();
        child2.Fitness();

        // aggiunta primo figlio
        future_generation.push_back(child1); 
        // se non è stato raggiunto il limite, aggiunta secondo figlio
        if(future_generation.size() < static_cast<size_t>(_npop)){
            future_generation.push_back(child2);
        }
    }    
    // salvo nuova popolazione
    _individuals=future_generation;
    FitnessSorting();

    // calcolo loss del best individual nella nuova generazione
    double best_loss = _individuals[0].GetFitness();
    // calcolo loss media del 50% della popolazione per vedere convergenza
    double avg_loss = GetBestHalfAvg();
    _loss << _gen_index << "\t" << best_loss << "\t" << avg_loss << endl;

    if(_gen_index % 100 == 0){
        cout << _type << " generation " << _gen_index << "\tbest loss = " << best_loss << endl;
    } 
    _gen_index++;
}

double Population::GetBestHalfAvg(){
    double sum=0;
    // loop sulla prima metà + aggiunta rispettivo fitness
    for(int i=0; i<_npop/2; i++){
        sum += _individuals[i].GetFitness();
    }
    // divido sum per metà della grandezza della popolazione
    return sum/(0.5*double(_npop));
}

int Population::GetNumGen(){
    // numero totale popolazioni da runnare
    return _ngen;
}

void Population::Closure(){
    if(_loss.is_open()) _loss.close();

    // estrazione genoma del miglior individuo
    uvec best_path = _individuals[0].GetPath();

    ofstream path(_type + "_best_path.dat");

    for(int i=0; i<_num_of_cities; i++){
        int j=best_path(i);
        path << _pos(j,0) << "\t" << _pos(j,1) << endl;
    }
    // in chiusura riporto anche la prima città per chiudere il cerchio
    path << _pos(best_path(0), 0) << "\t" << _pos(best_path(0), 1) << endl;
    path.close();
    cout << "\nBest path found.\n" << endl;
}