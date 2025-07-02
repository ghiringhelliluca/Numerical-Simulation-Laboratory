/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

// Random Walk 3D

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
vector<double> single_lattice_rw( double a, int n_step, Random &rnd);
BlockResult blocking_method_radq(vector<double> r, int N, int L);
vector<double> single_continuum_rw( double a, int n_step, Random &rnd);


// main
int main (int argc, char *argv[]){

   Random rnd = random_gen();

   int n_rw = 1e4; // number of random walk 
   int n_step = 1e2; // step for each random walk
   int N = 1e2; // blocks for blocking method
   int L = n_rw/N; // number of elements for each block
   double a = 1.; // step size
  


   // random walk discrete
   // generate n_rw random walks and store the r^2 for each step
   vector<vector<double>> r2_vector;

   for(int i = 0; i < n_rw; i++){
      r2_vector.push_back(single_lattice_rw(a, n_step, rnd));
   }


   
   // // error propagation to compute the error of radq(mean(r^2))   

   // vector<double> mean_r2; 
   // vector<double> error_r2;
   // for(int i =0; i< n_step; i++){
   //    vector<double> r2_step;
   //    for(int j = 0; j < n_rw; j++){
   //       r2_step.push_back(r2_vector[j][i]);
   //    }

   //    BlockResult result_r = blocking_method(r2_step, N, L);
   //    mean_r2.push_back(result_r.mean.back()); // take last value of mean
   //    error_r2.push_back(result_r.error.back());
   // }
   
   // cout << mean_r2.size() << " " << error_r2.size() << endl;
   // ofstream out("random_walk.dat");
   // for (int i = 0; i < n_step; i++) {
   //    out << i << " " << sqrt(mean_r2[i]) << " " << error_r2[i]/(2 * sqrt(mean_r2[i]))  << endl;
   //    cout << i << " " << sqrt(mean_r2[i]) << " " << error_r2[i]/(2 * sqrt(mean_r2[i]))  << endl;

   // }
   // out.close();



   // approssimazione radice media di r^2 media r 
   // Grazie all'approssimazione lineare, la varianza di r block è data dallo stesso fattore derivato dalla propagazione degli errori.
   vector<double> mean_r;
   vector<double> error_r;

   for(int i =0; i< n_step; i++){
      vector<double> r2_step;
      for(int j = 0; j < n_rw; j++){
         r2_step.push_back(r2_vector[j][i]);
      }

      BlockResult result_r = blocking_method_radq(r2_step, N, L);
      mean_r.push_back(result_r.mean.back()); // take last value of mean
      error_r.push_back(result_r.error.back());
   }
   
   ofstream out("random_walk_lattice.dat");
   for (int i = 0; i < n_step; i++) {
      out << i << " " << mean_r[i] << " " << error_r[i]  << endl;
   }
   out.close();




   // random walk continuum

   vector<vector<double>> r2_vector_c;

   for(int i = 0; i < n_rw; i++){
      r2_vector_c.push_back(single_continuum_rw(a, n_step, rnd));
   }

   vector<double> mean_r_c;
   vector<double> error_r_c;

   for(int i =0; i< n_step; i++){
      vector<double> r2_step_c;
      for(int j = 0; j < n_rw; j++){
         r2_step_c.push_back(r2_vector_c[j][i]);
      }

      BlockResult result_r_c = blocking_method_radq(r2_step_c, N, L);
      mean_r_c.push_back(result_r_c.mean.back()); // take last value of mean
      error_r_c.push_back(result_r_c.error.back());
   }
   
   ofstream outs("random_walk_continuum.dat");
   for (int i = 0; i < n_step; i++) {
      outs << i << " " << mean_r_c[i] << " " << error_r_c[i]  << endl;
   }
   outs.close();

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


vector<double> single_lattice_rw( double a, int n_step, Random &rnd) {

   vector<double> position(3, 0.0);
   vector<double> r2(n_step, 0.0);

   for (int i = 0; i < n_step; i++) {
      int direction = (int)(rnd.Rannyu() * 6); // Random number between 0 and 5
      if (direction == 0) position[0] += a;  // +x
      if (direction == 1) position[0] -= a;  // -x
      if (direction == 2) position[1] += a;  // +y
      if (direction == 3) position[1] -= a;  // -y
      if (direction == 4) position[2] += a;  // +z
      if (direction == 5) position[2] -= a;  // -z
      r2[i] = position[0]*position[0] + position[1]*position[1] + position[2]*position[2];
   }

   return r2;
}


vector<double> single_continuum_rw( double a, int n_step, Random &rnd) {

   vector<double> position(3, 0.0);
   vector<double> r2(n_step, 0.0);

   for (int i = 0; i < n_step; i++) {

      double phi = 2* M_PI * rnd.Rannyu();
      double theta = acos(1 - 2 * rnd.Rannyu());

      position[0] += a * sin(theta) * cos(phi);
      position[1] += a * sin(theta) * sin(phi);
      position[2] += a * cos(theta);

      r2[i] = position[0]*position[0] + position[1]*position[1] + position[2]*position[2];
   }

   return r2;
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
      av2[i] = pow(av[i],2); //A_i ^2 quadrato della media per blocco i-esimo dei r[k], la userò per calcolare sigma <A^2>
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



BlockResult blocking_method_radq(vector<double> r, int N, int L){ //N numero di blocchi, L numero di elementi per blocco

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
      av[i] = sqrt( av[i] / L);  // radice della media per blocco i-esimo 
      av2[i] = pow(av[i],2); // A_i ^2 quadrato della media per blocco i-esimo dei r[k], la userò per calcolare sigma <A^2>
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
