#include <iostream>
#include <lammps.h>
#include <library.h>
#include <sstream>
#include <mpi.h>
#include <cstring>

#ifndef INTERFACE_LMP_H
#define INTERFACE_LMP_H

using namespace std;

class Interface_lmp
{
public:
   
    MPI_Comm comm_lammps;
    LAMMPS_NS::LAMMPS *lmp;
    int lammps;    
    int myProc;

    Interface_lmp(int argc, char **argv);

    void initiate_lmp(int argc, char **argv);
    void load_bond(int bond_type, int new_id1, int new_id2);
    void unload_bond(int bond_type, int old_id1, int old_id2);
    void update_bonds(int bond_type, bool add_link, bool delete_link, int add_link_i, int add_link_j, int delete_link_i, int delete_link_j);
    void minimize();
    void run_dynamics(int steps);
    void print_bonds();
    void write_data(string line);
    void close_lmp();

private:
    
    bool isFirstRun;
};

#endif
