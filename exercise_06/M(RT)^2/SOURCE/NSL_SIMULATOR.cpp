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
#include "system.h"

using namespace std;

// code to simulate the system at different temperatures
// and measure the properties of the system
// the system is equilibrated at T = 2.0 and then the properties are measured

int main (int argc, char *argv[]){

  System SYS;
  SYS.initialize();  
  SYS.initialize_properties_output();
  SYS.block_reset(0);  

  for(int i = 0; i < 16; i++){
    double T = 2.0 - i*0.1;
    SYS.set_temperature(T); 

    // Equilibrate system only for T = 2.0
    if(T == 2.0){
      for(int j = 0; j < SYS.get_eqsteps(); j++){
        SYS.step();
      }
    }

    //  simulation
    for(int j = 0; j < SYS.get_nbl(); j++){
      for(int k = 0; k < SYS.get_nsteps(); k++){
        SYS.step();
        SYS.measure();
      }
      SYS.averages_final(j+1, T);
      SYS.block_reset(j+1);
    }

    SYS.finalize(T); 
  }

  return 0;
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
