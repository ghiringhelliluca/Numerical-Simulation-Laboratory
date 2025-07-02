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


Random random_gen();

// main
int main (int argc, char *argv[]){

   Random rnd = random_gen();

   int n = 1e4;
   vector<int> N = {1, 2 , 10 ,100};
   double lambda = 1.;
   double  mu = 0.;
   double gamma = 1.;
   
   for(int k = 0; k < N.size(); k++) {

      string filename = "output_N" + to_string(N[k]) + ".dat";
      ofstream out(filename);

      for(int j = 0; j < n; j++) {
         double sum_unif = 0.;
         double sum_exp = 0.;
         double sum_lorentz = 0.;

         for(int i = 0; i < N[k]; i++) {
            sum_unif += rnd.Rannyu();
            sum_exp += -1./lambda * log(1. - rnd.Rannyu());
            sum_lorentz += gamma * tan(M_PI * (rnd.Rannyu() - 0.5)) + mu;
         }

         out << sum_unif / N[k] << ' ' << sum_exp / N[k] << ' ' << sum_lorentz / N[k] << endl;
      }

      out.close();
   }

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

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
