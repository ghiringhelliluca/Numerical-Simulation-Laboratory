/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

// RANNYU: genera numeri casuali tramite x_n+1 = (a * x_n + c) mod m
// seed sono quattro numeri tra 0 e 4096, mi servono per genrare x_0
// Cambiare seed ⇒ Cambiare sequenza casuale
// primes sono due coppie di numeri che sevrvono per impostare c
// ci sono coppie di primes in modo che le seqeunze siano tra loro indipendenti
// Se cambi i due numeri primi, il generatore ha c diverso e posso generare sequenze diverse ma indipendenti

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
BlockResult blocking_method(vector<double> r, int N, int L);

// main
int main (int argc, char *argv[]){

   Random rnd = random_gen();

   int M = 100000;
   int N = 100;
   int L = M/N;
   vector<double> r;
   vector<double> sigma;
   vector<double> mean_sigma;
   vector<double> mean_r;
   vector<double> error_r;
   vector<double> error_sigma;


   for(int i=0; i<M; i++){
      r.push_back(rnd.Rannyu());
   }

   for(int i=0; i<M; i++){
      sigma.push_back(pow((r[i]-0.5),2));
   }

   BlockResult result_r = blocking_method(r, N, L);
   BlockResult result_sigma = blocking_method(sigma, N, L);
   for (int i = 0; i < N; i++) {
      mean_r.push_back(result_r.mean[i]);
      error_r.push_back(result_r.error[i]);
      mean_sigma.push_back(result_sigma.mean[i]);
      error_sigma.push_back(result_sigma.error[i]);

   }

   ofstream out("output_r.dat");
   for(int i = 0; i <N; i++){
   //   cout << i * L << " " << mean_r[i] - 1./2 << " " << error_r[i] << endl;
      out << i * L << " " << mean_r[i] - 1./2 << " " << error_r[i] << endl;
   }
   out.close();

   ofstream outs("output_sigma.dat");
   for(int i = 0; i <N; i++){
    //  cout << i * L << " " << mean_sigma[i] - 1./12 << " " << error_sigma[i] << endl;
      outs << i * L << " " << mean_sigma[i] - 1./12 << " " << error_sigma[i] << endl;
   }
   outs.close();

   return 0;
}


// function declaration

Random random_gen(){
   Random rnd;
   int seed[4];

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


double error(double AV, double AV2, int n){
    if (n==0){
        return 0;
      } else {
        return sqrt((AV2 - AV*AV)/n);
      }
}


BlockResult blocking_method(vector<double> r, int N, int L){


    vector<double> av(N,0);
    vector<double> av2(N,0);
      vector<double> mean;
      vector<double> err;


    for(int i = 0; i < N; i++){
      av[i] = 0;
      av2[i] = 0;
      for(int j = 0; j < L; j++){
         int k = j + L*i;
         av[i] += r[k];      
         }
      av[i] /= L;  // media per blocco i-esimo dei r[k] 
      av2[i] = pow(av[i],2); //  quadrato della media per blocco i-esimo dei r[k], la userò per calcolare sigma A
    }

    for(int i = 0; i< N; i++){
      double sum = 0;
      double sum2 = 0;
      for(int j = 0; j < i+1; j++){
         sum += av[j]; // somma delle medie per blocco fino al blocco i-esimo
         sum2 += av2[j]; // somma dei quadrati delle medie per blocco fino al blocco i-esimo
      }
      mean.push_back( sum/(i+1)); // media delle medie per blocco fino al blocco i-esimo, che equivale alla media dei r[k] fino al k-esimo
      err.push_back(error(mean[i], sum2/(i+1), i)); //i perchè devo fare / N-1
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
