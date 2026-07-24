#include "individual.h"
#include <iostream>

using namespace std;
using namespace arma;

//-----------------------------------------------------
// ho problemi per cat utilizzando Rannyu(1,_num_of_cities) 
// che dovrebbe produrre un numero in [1, _num_of_cities) ma usando poi
// arma::uvec e Swap, se anche solo una volta j==_num_of_cities, allora _path(j) non esiste e il programma crasha
// devo proteggere inizializzazione e mutazioni

// funzione di sicurezza per generazione casuale interi per intervalli inclusivi
int RandInt(Random& rnd, int min, int max_included){
    if(max_included < min){
        cerr << "Invalid random range: [" << min << ", " << max_included << "]" << endl;
        exit(1);
    }
    return min + int(rnd.Rannyu() * (max_included - min + 1));
}
//-----------------------------------------------------

Individual::Individual(int num_of_cities, const mat &distance_matrix){
    _num_of_cities = num_of_cities;
    _distance_matrix = distance_matrix;
    _path.set_size(_num_of_cities);
}

// Setto il primo path e lo cambio
void Individual::Initialize(Random & rnd){
    _path.set_size(_num_of_cities);
    // esclusione prima città che deve rimaner fissa
    for(int i=0; i<_num_of_cities; i++){
        _path(i) = i;
    }
    for(int i=0; i<_num_of_cities; i++){
        //int j=static_cast<int>(rnd.Rannyu(1,_num_of_cities));
        //int k=static_cast<int>(rnd.Rannyu(1,_num_of_cities));

        // sostituisco con qualcosa di più protetto, produzione valori in [1, _num_of_cities-1]
        // mai = num_of_cities
        int j = RandInt(rnd, 1, _num_of_cities - 1);
        int k = RandInt(rnd, 1, _num_of_cities - 1);
        Swap(j, k, "Initialize");
    }
    Check();
    Fitness();
}

// Controllo che la prima città stia in zero e che tutte le città compaiano una sola volta
void Individual::Check(){
    if(_path(0) !=0){
        cerr << "PROBLEM: the first city isn't zero." << endl;
        exit(1);
    }
    uvec singles = unique(_path);
    if(singles.n_elem != _num_of_cities){
        cerr << "PROBLEM: some cities have been visited multiple times.";
        exit(1);
    }
}

// Distanza percorsa/cammino <-> fitness
void Individual::Fitness(){
    double dist=0.;
    for(int i=0; i<(_num_of_cities - 1); i++){
        dist += _distance_matrix(_path(i), _path(i+1));
    }
    // periodic boundary (inclusione ritorno a città iniziale)
    dist += _distance_matrix(_path(_num_of_cities-1), _path(0)); 
    _fitness = dist;
}

// Swap di due elementi
void Individual::Swap(int i, int j, const std::string& caller){
    if(i < 0 || i >= _num_of_cities || j < 0 || j >= _num_of_cities){
        cerr << "PROBLEM: Swap indices out of bounds from [" << caller << "] "
             << "i = " << i << ", j = " << j
             << ", num_of_cities = " << _num_of_cities << endl;
        cerr << "_path size = " << _path.n_elem << endl;
        exit(1);
    }
    int temp = _path(i);
    _path(i) = _path(j);
    _path(j) = temp;
}

// Swap di due città random
void Individual::PairPermutation(Random & rnd){
    // devo escludere la prima
    //int i=static_cast<int>(rnd.Rannyu(1,_num_of_cities));
    //int j=static_cast<int>(rnd.Rannyu(1,_num_of_cities));
    //while (i==j) {
    //    j=static_cast<int>(rnd.Rannyu(1,_num_of_cities));
    //}
    int i = RandInt(rnd, 1, _num_of_cities - 1);
    int j = RandInt(rnd, 1, _num_of_cities - 1);

    while(i == j){
        j = RandInt(rnd, 1, _num_of_cities - 1);
    }

    Swap(i, j, "PairPermutation");
    Check();
}

// Shift di un blocco intero di "blk" città di l posizini 
void Individual::Shift(Random & rnd){
    // il blocco blk deve essere composto da almeno 2 città e più corto della sequenza completa
    //int blk=static_cast<int>(rnd.Rannyu(2,_num_of_cities-2));

    // la lunghezza dello shift non deve essere maggiore del numeero di città - blk
    //int l=static_cast<int>(rnd.Rannyu(1,_num_of_cities-blk));

    // nel riordinamento immagino i come nuovo pto di partenza + blk + l a seguire
    // l'ultima città può essere al massimo _num_of_cities - blk - l
    //int i=static_cast<int>(rnd.Rannyu(1,_num_of_cities - blk - l));
    int max_blk = _num_of_cities - 2;  
    int blk = RandInt(rnd, 2, max_blk);

    int max_l = _num_of_cities - 1 - blk;
    int l = RandInt(rnd, 1, max_l);

    int max_i = _num_of_cities - blk - l;
    int i = RandInt(rnd, 1, max_i);

    uvec mut = _path;
    for (int k=0; k<blk; k++){
        _path(i + l + k) = mut(i + k);
    }
    for(int k=0; k<l; k++){
        _path(i + k) = mut(i + k + blk);
    }
    Check();
}

// Swap di due blocchi di grandezza blk vicini
void Individual::AdiacentPermutation(Random & rnd){
    // i blocchi min di due città e max di metà del massimo numero di città
    //int blk=static_cast<int>(rnd.Rannyu(2, _num_of_cities / 2));

    // pto iniziale dopo prima città (fissata)
    // non oltre (_num_of_cities - 2*blk)  
    //int i=static_cast<int>(rnd.Rannyu(1, _num_of_cities - 2*blk));

    // chiaramente il secondo blocco parte dopo il primo e deve essere max (_num_of_cities - blk+1)
    //int j=static_cast<int>(rnd.Rannyu(i + blk, _num_of_cities - blk + 1));
    int max_blk = (_num_of_cities - 1) / 2;
    int blk = RandInt(rnd, 2, max_blk);

    int max_i = _num_of_cities - 2*blk;
    int i = RandInt(rnd, 1, max_i);

    int j = i + blk;

    if(i < 1 || j < 1 || j + blk - 1 >= _num_of_cities){
        cerr << "AdiacentPermutation invalid:"
             << " blk=" << blk
             << " i=" << i
             << " j=" << j
             << " num=" << _num_of_cities << endl;
        exit(1);
    }

    for(int k = 0; k < blk; k++){
        Swap(i + k, j + k, "AdiacentPermutation");
    }
    Check();
}

// Inversione ordine città entro blocco random
void Individual::Inversion(Random & rnd){
    //int i=static_cast<int>(rnd.Rannyu(1, _num_of_cities - 3));
    //int j=static_cast<int>(rnd.Rannyu(2, _num_of_cities - i));

    int i = RandInt(rnd, 1, _num_of_cities - 3);
    int max_len = _num_of_cities - i;
    int j = RandInt(rnd, 2, max_len);

    int k = i + j - 1;

    if(i < 1 || k >= _num_of_cities){
        cerr << "Inversion invalid:"
             << " i=" << i
             << " j=" << j
             << " k=" << k
             << " num=" << _num_of_cities << endl;
        exit(1);
    }

    for(int l = 0; l < j/2; l++){
        Swap(i + l, k - l, "Inversion");
    }
    Check();
}

void Individual::SetPath(const uvec & path){
    _path = path;
    this->Check();
    this->Fitness();
}