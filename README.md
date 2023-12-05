# loopExtrusion

C++ wrapper that calls LAMMPS as a library and performs loop extrusion of cohesin along chromatin.

------------------------
------- CONTENTS -------
------------------------

- *loopExtrusion.cpp* is the executable, the main script which initialise LAMMPS, the Gillespie algorithm and launch the simulation.

- *extrusion.cpp/extrusion.h* define the C++ class which drives the Gillespie algorithm simulating the extrusion process

- *parameters.cpp/parameters.h* define the C++ class which reads the parameters of the simulation.

- *interface_lmp.cpp/interface_lmp.h* define the C++ class which calls LAMMPS as a library and update the simulation according to the Gillespie algorithm.

- *test.tar* contains the files to run an example simulation (read the 'RUNNING THE TEST SIMULATION' section below). 


----------------------------
------- INSTALLATION -------
----------------------------

**BUILDING LAMMPS:**

- Download and unpack tarball from LAMMPS site (tested for lammps-29Sep2021).
- Inside the LAMMPS folder, run the following commands:

    mkdir build

    cd build

    cmake3 -D PKG_BROWNIAN=yes -D PKG_MOLECULE=yes -D PKG_OPENMP=yes -D PKG_REPLICA=yes -D PKG_EXTRA-DUMP=yes -D PKG_MISC=yes -D PKG_USER-MISC=yes -D BUILD_SHARED_LIBS=yes -D LAMMPS_INSTALL_RPATH=on ../cmake

    make -j 8

    make install

- Finally, add to your .bashrc file the following line: 

    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/$HOME/lammps-29Sep2021/build

Now you should have a LAMMPS installation that supports being called by external code.

**COMPILING:**

- Download and unpack the 'loopExtrusion' folder.
- Modify *Makefile* with your own directories 
- Run the 'make' command inside the folder to compile.

Now you should have the loopExtrusion executable.

**RUNNING THE TEST SIMULATION:** 

- Unpack the 'test' folder.
- 'polymer.data' is a LAMMPS configuration file for a homo-polymer.
- 'polymer.lam' is a LAMMPS input file for Langevin dynamics and hard-core repulsion between beads.
- 'param.in' is a text file defining the parameters for loop extrusion of the chain.
- Run the following command: 
 
    $PATH/loopExtrusion param.in polymer.lam 
    
----------------------
----- PARAMETERS -----
----------------------

A brief explanation of the parameters you can set for the loop extrusion simulation:

- *time_max* (double): total time of simulation
- *timestep* (double): timestep of integration
- *k_binding* (double): rate of loading of extruders (default=0)
- *k_unbinding* (double): rate of unloading of extruders (default=0)
- *k_step* (double): rate of movement of extruders (default=0)
- *k_cross_ctcf* (double): rate of crossing of a CTCF site (default=0)
- *n_extr_tot* (int): maximum number of extruders available (default=-1, i.e. unlimited extruders available)
- *n_extr_max* (int): maximum number of active extruders on the chain (default=0)
- *seed* (int): seed for the generation of random numbers (default=-1, i.e. the seed is generated)
- *debug*: activate debug mode, which prints real-time information about the extrusion process (default=False)
- *allow_overcome*: allows the extruders to cross themselves (default=False)
- *screen*: output of LAMMPS is printed in the terminal (default=False)
- *stride_log* (int): print output every *stride_log* Gillespie iterations (default=-1, i.e. don't print output)
- *state_file* (str): file with info on active extruders at the start of the simulation
- *ctcf_file* (str): file with positions and type of ctcf sites
