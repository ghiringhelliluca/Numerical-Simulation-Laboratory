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
BlockResult blocking_method(vector<double> r, int N, int L);
double asset_price(double S0, double r, double sigma, double ti, double tf, Random &rnd);
double N(double x);

// main
int main (int argc, char *argv[]){

   Random rnd = random_gen();

   int M = 100000;
   int N = 100;
   int L = M/N;

   double S0 = 100.;
   double K = 100.;
   double T = 1.;
   double r = 0.1;
   double sigma = 0.25;
   double t=0 ;

   vector<double> C_directly;
   vector<double> P_directly;
   vector<double> C_discretized;
   vector<double> P_discretized;

   double t_new = 0.;
   double t_old = 0.;
   int step = 100;
   double t_step = T/step;

   for(int i = 0; i < M; i++){
      double S = asset_price(S0, r, sigma, t, T, rnd);
      double C = exp(-r * T) * max(0., S - K);
      double P = exp(-r * T) * max(0., K - S);
      C_directly.push_back(C);
      P_directly.push_back(P);

      t_new = 0.;
      double S_discretized = S0;

      for(int j = 0; j < step; j++){
         t_old = t_new;
         t_new += t_step;

         S_discretized = asset_price(S_discretized, r, sigma, t_old, t_new, rnd);
      }

      C = exp(-r * T) * max(0., S_discretized - K);
      P = exp(-r * T) * max(0., K - S_discretized);
      C_discretized.push_back(C);
      P_discretized.push_back(P);
   }

   BlockResult result_C_directly = blocking_method(C_directly, N, L);
   BlockResult result_P_directly = blocking_method(P_directly, N, L);
   BlockResult result_C_discretized = blocking_method(C_discretized, N, L);
   BlockResult result_P_discretized = blocking_method(P_discretized, N, L);


   ofstream outd("output_C_directly.dat");
   ofstream out("output_C_discretized.dat");
   for (int i = 0; i < N; i++) {
      outd << i * L << " " << result_C_directly.mean[i] << " " << result_C_directly.error[i] << endl;
      out << i * L << " " << result_C_discretized.mean[i] << " " << result_C_discretized.error[i] << endl;
   }
   outd.close();
   out.close();

   ofstream outP("output_P_discretized.dat");
   ofstream outdP("output_P_directly.dat");
   for (int i = 0; i < N; i++) {
      outP << i * L << " " << result_P_discretized.mean[i] << " " << result_P_discretized.error[i] << endl;
      outdP << i * L << " " << result_P_directly.mean[i] << " " << result_P_directly.error[i] << endl;
   }
   outP.close();
   outdP.close();
   
   return 0;
}




// function declaration

double N(double x){
   return 0.5 * erfc(-x/sqrt(2)); // erfc(x) = 1 - erf(x)
}

double asset_price(double S0, double r, double sigma, double ti, double tf, Random &rnd){
   return S0 * exp((r - 0.5 * pow(sigma,2)) * (tf - ti) + sigma * rnd.Gauss(0, 1) * sqrt(tf - ti));
}

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
