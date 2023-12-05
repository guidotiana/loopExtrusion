#include "parameters.h"

using namespace std;

/////////////////////////////////////////////
// Parameters constructor
/////////////////////////////////////////////
Parameters::Parameters( int argc, char **argv )
{
     string fileName; 

     if ( argc != 3 ) Error("Name of parameters file not provided");
     fileName = argv[1];

     // defaults
     verbose = false;
     length = 0;
     k_binding = 0.;
     k_unbinding = 0.;     
     k_step = 0.;  
     k_cross_ctcf = 0.;
     allow_overcome = false;
     n_extr_tot = -1;
     seed = -1;
     stride_log = -1;
     n_extr_max = 0.;
     screen = false;
     debug = false;

     // read file
     ReadFile( fileName );
}


/////////////////////////////////////////////
// Read parameters from file
/////////////////////////////////////////////
bool Parameters::ReadFile( string fileName )
{
     ifstream fin;
     string line;
     vector<string> word;


     // open file
     fin.open(fileName, ios::in);
     if ( !fin.is_open() ) Error("Cannot open file "+fileName+" for reading parameters");

     // read file
     while( getline( fin, line ) )
     {
        if ( !line.empty() ) 
        { 
           Split(line, " \n\t=", word); 
           if ( word[0] == "time_max" ) time_max = stod( word[1] );
           if ( word[0] == "timestep" ) timestep = stod( word[1] ); 
           if ( word[0] == "stride_log" ) stride_log = stoi( word[1] ); 
           if ( word[0] == "k_binding" ) k_binding = stod( word[1] ); 
           if ( word[0] == "k_unbinding" ) k_unbinding = stod( word[1] ); 
           if ( word[0] == "k_step" ) k_step = stod( word[1] ); 
           if ( word[0] == "k_cross_ctcf" ) k_cross_ctcf = stod( word[1] ); 
           if ( word[0] == "verbose" ) verbose = true; 
           if ( word[0] == "debug" ) debug = true;
           if ( word[0] == "screen" ) screen = true; 
           if ( word[0] == "allow_overcome" ) allow_overcome = true; 
           if ( word[0] == "n_extr_tot" ) n_extr_tot = stoi( word[1] ); 
           if ( word[0] == "seed" ) seed = stoi( word[1] ); 
           if ( word[0] == "length" ) length = stoi( word[1] ); 
           if ( word[0] == "n_extr_max" ) n_extr_max = stoi( word[1] ); 
           if ( word[0] == "ctcf_file" ) ctcf_file = word[1];
           if ( word[0] == "state_file" ) state_file = word[1];
        } 
     }


     // write log
     if  (verbose )
     {
        
        Welcome( VERSION );
        cout << "Parameters read from "+fileName << endl;
        cout << "time_max          = "+to_string(time_max) << endl;
        cout << "timestep          = "+to_string(timestep) << endl;
        cout << "stride_log        = "+to_string(stride_log) << endl;
        cout << "length            = "+to_string(length) << endl;
        cout << "k_binding         = "+to_string(k_binding) << endl;
        cout << "k_unbinding       = "+to_string(k_unbinding) << endl;
        cout << "k_step            = "+to_string(k_step) << endl;
        cout << "k_cross_ctcf      = "+to_string(k_cross_ctcf) << endl;
        cout << "allow_overcome    = "+BoolToString(allow_overcome) << endl;
        cout << "n_extr_tot        = "+to_string(n_extr_tot) << endl;
        cout << "n_extr_max        = "+to_string(n_extr_max) << endl;
        cout << "seed              = "+to_string(seed) << endl;
        cout << "debug             = "+BoolToString(debug) << endl;
        if ( !ctcf_file.empty() ) cout << "ctcf_file         "+ctcf_file << endl;
        if ( !state_file.empty() ) cout << "state_file        = "+state_file << endl;        
        cout << endl;
     }

     // checks
     if (length < 1) Error("The length of the chain mast be larger than 1");
     if (time_max<1E-15) Error("You must define time_max in the parameter file");
     if (timestep<1E-15) Error("You must define timestep in the parameter file");

     // Warnings
     if (time_max <= 2.3/k_binding){
        // 2.3 ~ -log(0.1), in this case there is a probability 0.1 that no extruders will load.  
        cout << "+++++++++++++++++++++" << endl;
        cout << "WARNING: the binding time of extruders seems large for your simulation time. Probable absence of loop-extrusion with these parameters!" << endl;
        cout << "The simulation may also misbehave." << endl;
        cout << "+++++++++++++++++++++" << endl;
        cout << endl;
     }

     fin.close();

     return true;

}

/////////////////////////////////////////////
// Exit with error
/////////////////////////////////////////////
void Parameters::Error( string message )
{
     time_t now = time(0);
     char* dt = ctime(&now);

     cerr << endl;
     cerr << "*** Execution error:" << endl;
     cerr << message << endl;
     cerr << dt << endl;

     exit(1);
}

/////////////////////////////////////////////
// Split line into words
/////////////////////////////////////////////
void Parameters::Split(const string& s, string c, vector<string>& v) 
{
   string::size_type i = 0;
   string::size_type j = 0;
   string w;

   v.clear();
   j = s.find_first_of(c);
   if ( j ==  string::npos ) { v.push_back(s); return; }


   while (j != string::npos) 
   {
      v.push_back(s.substr(i, j-i));
      i = ++j;
      j = s.find_first_of(c, j);

      if (j == string::npos)
      {
         w = s.substr(i, s.length());
         if ( ! w.empty() ) v.push_back(w);
      }
   }
}

/////////////////////////////////////////////
// Welcome screen
/////////////////////////////////////////////
void Parameters::Welcome( string version )
{
   time_t now = time(0);
   char* dt = ctime(&now);
   pid_t pid = getpid();

   cout << "*****************************************************" << endl;
   cout << "*                                                   *" << endl;
   cout << "*             loopExtrusion v. "+version+"                  *" << endl;
   cout << "*                                                   *" << endl;
   cout << "*****************************************************" << endl;
   cout << " Started: " << dt ;
   cout << " pid: " << pid << endl;
   cout << endl;
}

/////////////////////////////////////////////
// Cast bool into string
/////////////////////////////////////////////
string Parameters::BoolToString(bool b)
{
  return b ? "true" : "false";
}
