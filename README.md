# loopExtrusion

C++ wrapper that calls LAMMPS as a library and performs loop extrusion of cohesin along chromatin.

----------------------------
------- INSTALLATION -------
----------------------------

BUILDING LAMMPS:

- Download and unpack tarball from LAMMPS site (tested for lammps-29Sep2021).
- Inside the LAMMPS folder, run the following commands:

    mkdir build

    cd build

    cmake3 -D PKG_BROWNIAN=yes -D PKG_MOLECULE=yes -D PKG_OPENMP=yes -D PKG_REPLICA=yes -D PKG_EXTRA-DUMP=yes -D PKG_USER-LE=yes -D PKG_MISC=yes -D PKG_USER-MISC=yes -D BUILD_SHARED_LIBS=yes -D LAMMPS_INSTALL_RPATH=on ../cmake

    make -j 8

    make install

- Finally, add to your .bashrc file the following line: 

    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/$HOME/lammps-29Sep2021/build

Now you should have a LAMMPS installation that supports being called by external code.

COMPILING:

- Download and unpack the 'loopExtrusion' folder.
- Run the 'make' command inside the folder to compile.

Now you should have the loopExtrusion executable.

RUNNING THE TEST SIMULATION: 

- Unpack the 'test' folder.
- 'polymer.data' is a LAMMPS configuration file for a homo-polymer.
- 'polymer.lam' is a LAMMPS input file for Langevin dynamics and hard-core repulsion between beads.
- 'param.in' is a text file defining the parameters for loop extrusion of the chain.
- Run the following command: 
 
    $PATH/loopExtrusion param.in polymer.lam 

