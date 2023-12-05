#include <ctype.h>
#include <cstring>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <ctime>
#include <string>
#include <vector>
#include <sys/types.h>
#include <unistd.h>

using namespace std;

#define VERSION "0.1"


class Parameters
{

      public:

      int length;
      double time_max;
      double timestep;
      int stride_log;
      double k_binding;
      double k_unbinding;
      double k_step;  
      double k_cross_ctcf;
      bool verbose;
      bool debug;
      bool allow_overcome;
      bool screen;
      int n_extr_tot; 		// set to -1 to ignore
      int n_extr_max; 		
      int seed;
      string ctcf_file;
      string state_file;    

      Parameters( int, char ** );
      void Error( string );
      bool ReadFile( string );
      void Welcome( string version );


      private:

      void Split(const string& s, string c, vector<string>& v);
      string BoolToString(bool b);
      




};
