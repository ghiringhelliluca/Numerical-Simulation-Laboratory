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
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include "random.h"

using namespace std;



Random random_gen();
double compute_chi(vector<int> bin, int M, int n);

// main
int main (int argc, char *argv[]){

   Random rnd = random_gen();

   int M = 1e2;
   int n = 1e4;
   int N = 1e2;
   vector<double> chisquare; // vettore chi quadro

   for(int i =0; i<N; i++){
      vector<int> bins(M, 0); //  il vettore contatore
      for(int j = 0; j< n; j ++){
         // genero numero tra 0 e 1, lo motiplico per 100 e prendo il valore intero, quindi avrò nuermo intero tra 0 e 99
         double rand_num = rnd.Rannyu();
         int bin_index = static_cast<int>(rand_num * M); //  indice dell'intervallo, static_cast per convertire in int
         bins[bin_index]++;
      }
      chisquare.push_back(compute_chi(bins, M, n));
   }

   ofstream out("output_chisquare.dat");
   for(size_t i = 0; i <chisquare.size() ; i++){
      //cout << chisquare[i] << endl;
      out << chisquare[i] << endl;
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

double compute_chi(vector<int> bin, int M, int n){

   double E = double(n/M);
   double chi2 = 0.0;

   for (int i = 0; i < M; ++i) {
      double O = bin[i];  // valore osservato
      chi2 += std::pow(O - E, 2) / E;
   }

   return chi2;
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
