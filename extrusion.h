#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <ctype.h>
#include <ctime>
#include <cmath>
#include <iomanip>

#ifndef HPARAMETERS
#define HPARAMETERS
#include "parameters.h"
#endif

#define LARGE 999999
#define SMALL 1E-15
#define NREACT 4

using namespace std;


class Extrusion
{

   public:

     // input
     double k_binding;
     double k_unbinding;
     double k_step;  
     double k_cross_ctcf;
     bool allow_overcome;
     int n_extr_tot; 		// set to -1 to ignore
     int seed;

     // output
     double tau;
     bool add_link;
     int add_link_i;
     int add_link_j;
     bool delete_link;
     int delete_link_i;
     int delete_link_j;
     int n_extr_bound;		// how many extruders bound
     int (*extrList)[4];	// 0=i, 1=j, 2=time last move i, 3=time last move j
     string exitError;

     // functions
     Extrusion( Parameters );
     bool Event(bool debug);
     bool ReadCTCF(string fileName);
     bool PrintState(string fileName);
     bool ReadState(string fileName, bool debug);
     bool PrintMap( string fileName, bool asList, bool onlyExist );
     void CatchError( bool ok );


   private:


     int iTime;
     int length;
     int *ctcf;
     int nCTCF;
     int *occupiedSites;
     int n_extr_max;
     double propensities[NREACT+1];
     int **map;			// how many extruders between i and j
     string reaction_name[NREACT+1];


     // functions
     int **AlloArrayInt(int n);
     bool RandomBind( bool debug);
     bool RandomUnbind( bool debug);
     bool RandomStepForward(  bool ctcf_cross, bool debug);
     bool AddExtruder(int i, int j, int iTimeI, int iTimeJ);
     bool RemoveExtruder(int i, int j, int iTimeI, int iTimeJ);
     int iRand(int n);
     double DRand( void );
     bool LogicalXOR(bool a, bool b);
     void Randomize(int seed, bool debug);
     bool CalulatePropensities( bool debug);
     bool CheckStepOk( int w, int dir, bool ctcf_cross, bool debug);
     int SelectReaction( void );
     bool ApplyReaction(int r, bool debug);



};
     
