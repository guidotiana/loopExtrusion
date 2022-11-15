#include "extrusion.h"
#ifndef HPARAMETERS
#define HPARAMETERS
#include "parameters.h"
#endif

int main(int argc, char **argv)
{
    double time=0, tau;
    bool ok;
    int iStep=0;

    Parameters parm(argc, argv);

    Extrusion e( parm );

    // time loop of extruder Markov chain
    do
    {
       // Gillespie event
       ok = e.Event( parm.debug );
     
       // Update of links
       // here update of contacts using e.addlink, etc.

       // MD
       // here MD from time to (time + e.tau)
       time += e.tau;
       iStep ++;

       // print output
       if ( parm.stride_log>0 && !(iStep%parm.stride_log) ) cout << time << "\t\t" << e.n_extr_bound << endl;

    } while ( time < parm.time_max );

    return 0;
}
