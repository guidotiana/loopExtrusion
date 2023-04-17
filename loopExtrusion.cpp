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
    double tau_min=300, tau_0=0; //minimum time between dynamics runs 
    string data_line = "write_data last.data"; 
    
    //Reading Gillespie parameters
    Parameters parm(argc, argv);

    //Initializing extrusion algorithm
    Extrusion e( parm );

    //Initializing lammps and opening interface
    Interface_lmp inter_lmp(argc, argv); 
   
    //Reading CTCF sites
    e.ReadCTCF("ctcf_sites.data");
 
    //Main Gillespie loop    
    do
    {  
       while (tau_0 < tau_min)
       {
          // Gillespie event
          ok = e.Event( parm.debug );
          
          if (!ok){
             cout << "Binding probability is zero, no loop extrusion" << endl;
             tau_0 = parm.time_max;
          } 
          else {
             tau_0 += e.tau;
          }
 
          // Update of links
          inter_lmp.update_bonds(2, e.add_link, e.delete_link, e.add_link_i, e.add_link_j, e.delete_link_i, e.delete_link_j);
       }

       //Minimize energy of new configuration
       inter_lmp.minimize();
       
       //Molecular dynamics with LAMMPS from time to (time + e.tau)
       int time_left = parm.time_max-time;
       time += tau_0;
       if (ceil(tau_0)>parm.time_max){
          inter_lmp.run_dynamics(parm.time_max);
       }
       else if (ceil(tau_0)>time_left){
          inter_lmp.run_dynamics(time_left);
       }
       else inter_lmp.run_dynamics(ceil(tau_0));       

       //Update internal time variables
       iStep ++;
       tau_0 = 0;

       //Print output
       if ( parm.stride_log>0 && !(iStep%parm.stride_log) )
       {
          cout << "Time = " << time << "\t\t" << "# extruders = " << e.n_extr_bound << endl;
          inter_lmp.print_bonds();   
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
