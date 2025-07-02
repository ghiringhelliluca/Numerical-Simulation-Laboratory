/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

// Integro con campionamento uniforme e data blocking
// genero M numeri casuali tra 0 e 1 e li divido in N blocchi
// calcolo la media e l'errore con il data blocking, valuto l'integrale come f nei punti generati * (b - a), poi faccio data blocking sui valori medi di I 

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
double f(double x);
double f_importance(double x);

// main
int main (int argc, char *argv[]){

   Random rnd = random_gen();

   int M = 1e4;
   int N = 1e2;
   int L = M/N;
   int b = 1;
   int a = 0.;

   // UNIFORM SAMPLING
   vector<double> f_r;
   vector<double> mean_I;
   vector<double> error_I;

   for(int i =0; i<M; i++){
      f_r.push_back((b - a) * f(rnd.Rannyu()));
   }

   BlockResult result_r = blocking_method(f_r, N, L);
   for (int i = 0; i < N; i++) {
      mean_I.push_back(result_r.mean[i]);
      error_I.push_back(result_r.error[i]);
   }

   ofstream out("uniform_sampling.dat");
   for(int i = 0; i < N ; i++){
      out << i * L << " " << mean_I[i]  << " " << error_I[i] << endl;
   }
   out.close();


   // IMPORTANCE SAMPLING
   vector<double> f_r_importance;
   vector<double> mean_I_importance;
   vector<double> error_I_importance;

   for(int i =0; i< M; i++){
      double x = rnd.Rannyu();
      double x_importance = 1. - sqrt(1. - x); // cumulativa inversa di -2x + 2
      f_r_importance.push_back(f(x_importance) / f_importance(x_importance));
      
   }

   BlockResult result_r_importance = blocking_method(f_r_importance, N, L);
   for (int i = 0; i < N; i++) {
      mean_I_importance.push_back(result_r_importance.mean[i]);
      error_I_importance.push_back(result_r_importance.error[i]);
   }

   ofstream out_importance("importance_sampling.dat");
   for(int i = 0; i < N ; i++){
      out_importance << i * L << " " << mean_I_importance[i]  << " " << error_I_importance[i] << endl;
   }
   out_importance.close();

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


double f(double x) {
    return M_PI / 2.0 * cos(M_PI * x / 2.0);
}

double f_importance(double x) {
    return - 2.0 * x + 2.;
}


double error(double AV, double AV2, int n){
    if (n==0){
        return 0;
      } else {
        return sqrt((AV2 - AV*AV)/n);
      }
}


BlockResult blocking_method(vector<double> r, int N, int L){ //N numero di blocchi, L numero di elementi per blocco

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
   av2[i] = pow(av[i],2); // <A^2> quadrato della media per blocco i-esimo dei r[k], la userò per calcolare sigma A
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
