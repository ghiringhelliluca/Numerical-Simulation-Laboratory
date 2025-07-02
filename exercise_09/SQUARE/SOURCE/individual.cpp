/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <math.h>
#include <set>
#include <vector>
#include "individual.h"

using namespace std;
using namespace arma;



void Individual :: initialize(Random &_rnd, int sim_type){

   if(sim_type == 0){
      _ndim = 1; // circle
      _individual.set_size(_ncity);
      for(unsigned int i=0; i<_ncity; i++){
         _individual(i).resize(_ndim + 1);
         _individual(i)(0) = i + 1; // city id
         double angle = _rnd.Rannyu(0., 2 * M_PI); // random angle
         _individual(i)(1) = angle ; // random angles
      }
   } else if(sim_type == 1){
      _ndim = 2; // square
      _individual.set_size(_ncity );
      for(unsigned int i=0; i<_ncity; i++){
         _individual(i).resize(_ndim + 1);
         _individual(i)(0) = i + 1; // city id
         _individual(i)(1) = _rnd.Rannyu(0.0, _l);
         _individual(i)(2) = _rnd.Rannyu(0.0, _l);
      } 
   } else{
      cerr << "PROBLEM: unknown simulation type" << endl;
      exit(EXIT_FAILURE);
   }

   return;
}


// permutation function, make a pair permutation fix first city
void Individual :: permutation(Random &_rnd){
   for(int i=1; i<_ncity; i++){
      int j = int(_rnd.Rannyu(1, _ncity));
      if(i != j){
         vec temp = _individual(i);
         _individual(i) = _individual(j);
         _individual(j) = temp;
      }
   }
   return;
}


void Individual:: permutation_pair(Random &_rnd){
   int i = int(_rnd.Rannyu(2, _ncity));
   int j = int(_rnd.Rannyu(2, _ncity));
   if(i != j){
      vec temp = _individual(i);
      _individual(i) = _individual(j);
      _individual(j) = temp;
   }
   return;
}



double Individual :: loss(int sim_type){
   if(sim_type == 0){
      return loss_circle();
   } else if(sim_type == 1){
      return loss_square();
   } else{
      cerr << "PROBLEM: unknown simulation type" << endl;
      exit(EXIT_FAILURE);
   }
}

double Individual::loss_circle(){
   double loss = 0.0;
   for(unsigned int i=0; i<_ncity; i++){

      if(i == _ncity - 1){
         double dtheta = _individual(i)(1) - _individual(0)(1);
         dtheta = fmod(dtheta, 2*M_PI);
         if (dtheta > M_PI) dtheta = 2*M_PI - dtheta;
         if (dtheta < 0) dtheta += 2 * M_PI;
         loss += _r * dtheta;  
      } else{
         double dtheta = _individual(i)(1) - _individual(i+1)(1);
         dtheta = fmod(dtheta, 2*M_PI);
         if (dtheta > M_PI) dtheta = 2*M_PI - dtheta;
         if (dtheta < 0) dtheta += 2 * M_PI;         
         loss += _r * dtheta;   
      }

   }
   
   return loss;
}

double Individual::loss_square(){
   double loss = 0.0;
   for(unsigned int i=0; i<_ncity; i++){
      if(i == _ncity - 1){
         double dx = _individual(i)(1) - _individual(0)(1);
         double dy = _individual(i)(2) - _individual(0)(2);
         loss += sqrt(dx*dx + dy*dy);
      } else{
         double dx = _individual(i)(1) - _individual(i+1)(1);
         double dy = _individual(i)(2) - _individual(i+1)(2);
         loss += sqrt(dx*dx + dy*dy);
      }
   }
   return loss;
}



pair<Individual, Individual> Individual::crossover_pair(const Individual& partner, Random& _rnd) {
   Individual mom = *this;
   Individual dad = partner;

   vector<int> order1, order2;
   int cut = int(_rnd.Rannyu(1, _ncity));

   // mather: this
   for (int i = 0; i < cut; ++i){
      order1.push_back(mom._individual(i)(0));
   }

   for (int i = 0; i < _ncity; ++i) {
      int id = dad._individual(i)(0);
      // check if id is not already in order1
      if (find(order1.begin(), order1.end(), id) == order1.end()){ // find returns iterator to the first element equal to id, if not found returns end
         order1.push_back(id);
      }
   }

   // father: partner
   for (int i = 0; i < cut; ++i)
      order2.push_back(dad._individual(i)(0));

   for (int i = 0; i < _ncity; ++i) {
      int id = mom._individual(i)(0);
      if (find(order2.begin(), order2.end(), id) == order2.end())
         order2.push_back(id);
   }

   Individual child1, child2;
   child1._individual.set_size(_ncity);
   child2._individual.set_size(_ncity);

   // build child1._individual from order1
   // for each city id in order1, find the corresponding city in _individual and copy it to child1
   for (int i = 0; i < _ncity; ++i) {
      int city_id = order1[i];
      for (int j = 0; j < _ncity; ++j) {
         if (_individual(j)(0) == city_id) {
            child1._individual(i) = mom._individual(j);
            break;
         }
      }
   }

   // build child2._individual from order2
   for (int i = 0; i < _ncity; ++i) {
      int city_id = order2[i];
      for (int j = 0; j < _ncity; ++j) {
         if (_individual(j)(0) == city_id) {
            child2._individual(i) = dad._individual(j);
            break;
         }
      }
   }

   return make_pair(child1, child2);
}


void Individual::mutation_shift(Random &_rnd) {
   int max_start = _ncity - 2; // max start index (second last real city)
   if (max_start < 1) return;  // troppo poche città per shifare

   int m = int(_rnd.Rannyu(1, max_start + 1));  // inizio del blocco
   int max_block_len = min(4, _ncity - m);  // max blocco = 4 città (puoi aumentare qui!)
   int block_len = int(_rnd.Rannyu(1, max_block_len + 1));

   int max_shift = _ncity - m - block_len;
   if (max_shift < 1) return;

   int shift = int(_rnd.Rannyu(1, max_shift + 1));

   // --- Copia il percorso originale tranne il blocco da spostare ---
   vector<vec> new_path;

   // Prima parte (prima del blocco)
   for (int i = 0; i < m; ++i){
      new_path.push_back(_individual(i));
   }
   // Parte intermedia (tra blocco e shift)
   for (int i = m + block_len; i < m + block_len + shift; ++i){
      new_path.push_back(_individual(i));
   }
      

   // Inserisci il blocco spostato
   for (int i = 0; i < block_len; ++i){
       new_path.push_back(_individual(m + i));
   }
   // Parte finale
   for (int i = m + block_len + shift; i < _ncity; ++i){
       new_path.push_back(_individual(i));
   }
   // new_path.push_back(new_path[0]);

   // --- Copia nel field ---
   for (int i = 0; i < _ncity; ++i){
      _individual(i) = new_path[i];
   }
}



//inversion of the order in which they appear in the path of $m$ cities (except for the first city and $m \le N$), e.g. $\left[ 1, 2, 3, 4, 5 \right] \to \left[ 1, 4, 3, 2, 5 \right]$ for the inversion of the cities from 2 to 4.
void Individual::mutation_inversion(Random &_rnd) {

   int max_start = _ncity - 3; // max index start to invert at least 2 cities
   int start = int(_rnd.Rannyu(1, max_start + 1)); // 

   // end at least start+1, max _ncity - 1 
   int end = int(_rnd.Rannyu(start + 1, _ncity));

   // invert subvector
   while (start < end) {
      vec tmp = _individual(start);
      _individual(start) = _individual(end);
      _individual(end) = tmp;
      ++start;
      --end;
   }
   return ;
}


//permutation among $m$ contiguous cities (except for the first city) with other (different!) $m$ contiguous cities ($m<N/2$), e.g. $\left[ 1, 2, 3, 4, 5 \right] \to \left[ 1, 4, 5, 2, 3 \right] $ for a permutation of the second and third cities with the last 2.
void Individual::mutation_swap_blocks(Random &_rnd) {
   int m = int(_rnd.Rannyu(1, int( _ncity / 2) )); // block length
   int a = int(_rnd.Rannyu(1, _ncity - 2 * m)); // index start block 1
   int b = int(_rnd.Rannyu(a + m, _ncity - m + 1)); // index start block 2

   // swap blocks
   for(int i = 0; i < m; ++i){
      vec tmp = _individual(a + i);
      _individual(a + i) = _individual(b + i);
      _individual(b + i) = tmp;
   }
   return;
}


bool Individual::check(){
   // Check if first city is ID 1
   if(_individual(0)(0) != 1){
       cout << "Error: start and end city are not the same!" << endl;
       return false;
   }

   // Check all cities are only one time
   set<int> seen;
   for(int i = 0; i < _ncity; ++i){
       int city_id = _individual(i)(0);
       if(seen.find(city_id) != seen.end()){
           cout << "Error: city " << city_id << " appears more than once!" << endl;
           return false;
       }
       seen.insert(city_id);
   }

   // Check che tutte le città da 1 a n siano presenti
   for(int i = 1; i <= _ncity; ++i){
       if(seen.find(i) == seen.end()){
           cout << "Error: missing city " << i << endl;
           return false;
       }
   }

   return true;
}


void Individual::print_path(){

   for(int i = 0; i < _ncity; ++i){
      cout << _individual(i)(0) << " ";
   }
   cout << endl;
   return;
}

void Individual::save_path(ofstream& path){
   for(int i = 0; i < _ncity; ++i){
      path << _individual(i)(0) << " ";
   } 
   path << endl;
   return;
}

void Individual::save_individual(){
   ofstream pos_individual("position_cities.dat");

   pos_individual << "CITY_ID      ANGLE" << endl;
 
   for(int i =0; i < _ncity; i++){
     pos_individual << _individual(i)(0) << "             " << _individual(i)(1) << endl;
   }
}



/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
