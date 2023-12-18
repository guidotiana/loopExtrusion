#include "extrusion.h"
#include "interface_lmp.h"
#include <sstream>
#include <iostream>
#include <string>

#ifndef HPARAMETERS
#define HPARAMETERS
#include "parameters.h"
#endif

int main(int argc, char **argv)
{   
    //Defining variables
    double time=0, tau;
    bool ok;
    int iStep=0;
    double tau_0=0; //minimum time between dynamics runs 
    string data_line = "write_data last.data"; 
    
    //Reading Gillespie parameters
    Parameters parm(argc, argv);

    //Initializing extrusion algorithm
    Extrusion e( parm );

    //Reading CTCF sites
    e.ReadCTCF(parm.ctcf_file);

    //Reading state
    e.ReadState(parm.state_file, true); 

    //Initializing lammps and opening interface
    Interface_lmp inter_lmp(argc, argv, parm.screen); 
    
    //Setting integration timestep
    inter_lmp.set_timestep(parm.timestep);

    //Loading initial extruders in lammps
    for (int i=0; i<e.n_extr_bound; i++)
       {
         inter_lmp.load_bond(2, e.extrList[i][0]+1, e.extrList[i][1]+1);
       }

    //Main Gillespie loop    
    do
    {  
       while (tau_0 <= parm.tau_min)
       {
          // Gillespie event
          ok = e.Event( parm.debug );
          
          if (!ok){
             cout << "Binding probability is zero, no loop extrusion" << endl;
             tau_0 = parm.time_max;
          } 
          else {
             tau_0 += e.tau;
             // Update of links
             inter_lmp.update_bonds(2, e.add_link, e.delete_link, e.add_link_i, e.add_link_j, e.delete_link_i, e.delete_link_j);
          }
       }

       //Minimize energy of new configuration
       inter_lmp.minimize();
       
       //Molecular dynamics with LAMMPS from time to (time + e.tau)
       int time_left = parm.time_max-time;

       if (ceil(tau_0) > time_left){
          inter_lmp.run_dynamics(time_left/parm.timestep);
          time = parm.time_max;
       }
       else {
          inter_lmp.run_dynamics(ceil(tau_0)/parm.timestep);      
          time += tau_0;
       } 

       //Update internal time variables
       iStep ++;
       tau_0 = 0;

       //Print output
       if ( parm.stride_log>0 && !(iStep%parm.stride_log) )
       {
          cout << fixed;
          cout << "Time = " << time << "\t\t" << "# extruders = " << e.n_extr_bound << endl;
          inter_lmp.print_bonds(e.extrList, e.n_extr_bound);   
       }
    } while ( time < parm.time_max );
  
   cout << "Final number of extruders: " << e.n_extr_bound << endl; 
   
   //Writing final configuration
   inter_lmp.write_data(data_line); 
   
   //Closing LAMMPS
   inter_lmp.close_lmp();
    
   cout << "Done!" << endl;
    
   return 0;
}
