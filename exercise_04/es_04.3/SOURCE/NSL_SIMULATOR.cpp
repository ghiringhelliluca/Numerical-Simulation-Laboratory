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

int main (int argc, char *argv[]){

  int nconf = 1;
  System SYS;
  SYS.initialize();  //read input.dat file
  SYS.initialize_properties(); // properties to be measured and create output files in OUTPUT directory
  SYS.block_reset(0);  // reset block accumulators to zero

  // begin simulation
  for(int i=0; i < SYS.get_nbl(); i++){ // loop over blocks
    for(int j=0; j < SYS.get_nsteps(); j++){ // loop over steps in a block
      SYS.step(); // make a step, i.e. move a particle with Verlet
      SYS.measure(); // measure properties
    }
    SYS.averages(i+1); // compute averages and errors
    SYS.block_reset(i+1);
  }
  SYS.finalize();

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
