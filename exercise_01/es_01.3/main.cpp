/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

// GIUSTO

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include "random.h"

using namespace std;

struct BlockResult {
    vector<double> mean;
    vector<double> error;
};

Random random_gen();
double error(double AV, double AV2, int n);
BlockResult blocking_method(const vector<double>& data, int nBlocks);
double random_angle(Random &rnd);

// main
int main (int argc, char *argv[]){

   Random rnd = random_gen();

   double L = 1.0;      // Lunghezza dell'ago
   double d = 1.5;      // Distanza tra le linee 
   int M = 100000;      // Numero totale di lanci
   int N = 100;         // Numero di blocchi
   int throws_per_block = M / N;  // Numero di lanci per blocco

   vector<int> hits(M, 0);

   for (int i = 0; i < M; i++){
      // x: distanza del centro dell'ago dalla linea più vicina, uniforme in [0, d/2]
      double x = rnd.Rannyu(0, d/2.0);
      // theta: angolo uniforme in [0, pi]
      double theta = random_angle(rnd);  
      // ago interseca una linea se x <= (L/2) * sin(theta)
      if( x <= (L/2.0) * sin(theta) )
         hits[i] = 1;
      else
         hits[i] = 0;
   }
   
   vector<double> block_pi(N, 0.0);
   for (int i = 0; i < N; i++){
      int block_hits = 0;
      for (int j = 0; j < throws_per_block; j++){
         int index = i * throws_per_block + j;
         block_hits += hits[index];
      }
      if(block_hits == 0)
         block_pi[i] = 0.0;
      else
         block_pi[i] = (2.0 * L * throws_per_block) / (d * block_hits);
   }

   BlockResult result = blocking_method(block_pi, N);
   
   ofstream out("pi_estimate.dat");
   for (int i = 0; i < N; i++){
      out << throws_per_block * (i+1) << " " << result.mean[i] << " " << result.error[i] << endl;
   }
   out.close();

   return 0;
}


// function declaration

Random random_gen(){
   Random rnd;
   int seed[4];

   // prendo i primi due numeri primi dal file Primes
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

   return rnd;
}


// random angle from 0 to pi without using pi
double random_angle(Random &rnd) {
   double x, y;
   do{
      x = rnd.Rannyu(-1, 1);
      y = rnd.Rannyu(-1, 1);
   } while ((x*x + y*y > 1.0) | (y < 0.0));

   double r = sqrt(x*x + y*y);

   double theta;
   theta = acos( x / r) ;
   return theta ;
}


double error(double AV, double AV2, int n){
    if (n==0){
        return 0;
      } else {
        return sqrt((AV2 - AV*AV)/n);
      }
}

BlockResult blocking_method(const vector<double>& data, int nBlocks){
    vector<double> mean;
    vector<double> err;
    
    vector<double> block_av(nBlocks, 0.0);
    vector<double> block_av2(nBlocks, 0.0);
    for (int i = 0; i < nBlocks; i++){
        block_av[i] = data[i];
        block_av2[i] = data[i] * data[i];
    }
    
    for (int i = 0; i < nBlocks; i++){
        double sum = 0.0;
        double sum2 = 0.0;
        for (int j = 0; j <= i; j++){
            sum  += block_av[j];
            sum2 += block_av2[j];
        }
        double avg = sum / (i + 1);
        mean.push_back(avg);
        if(i == 0)
            err.push_back(0.0);
        else
            err.push_back( sqrt((sum2/(i + 1) - avg*avg) / i) );
    }


      BlockResult result;
      result.mean = mean;
      result.error = err;
      return result;
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
