#include "extrusion.h"



/////////////////////////////////////////////
// Extrusion constructor
/////////////////////////////////////////////
Extrusion::Extrusion(Parameters parm)
{
   if (parm.debug) cerr << "Initializing extrusion..." << endl;

   length = parm.length;

   // allocate memory
   map = AlloArrayInt(parm.length);
   extrList = new int [parm.n_extr_max][4];
   ctcf = new int [parm.length];
   for (int i=0;i<parm.length;i++) ctcf[i] = 0;
   nCTCF = 0;
   occupiedSites = new int [parm.length];
   for (int i=0;i<parm.length;i++) occupiedSites[i] = 0;

   // set private variables
   iTime = 0;
   n_extr_bound = 0;
   exitError = "";

   // set input from parameters 
   k_binding = parm.k_binding;
   k_unbinding = parm.k_unbinding;
   k_step = parm.k_step;
   k_cross_ctcf = parm.k_cross_ctcf;
   allow_overcome = parm.allow_overcome;
   n_extr_tot = parm.n_extr_tot;

   // set output defaults
   bool add_link = false;
   int add_link_i = -1;
   int add_link_j = -1;
   bool delete_link = false;
   int delete_link_i = -1;
   int delete_link_j = -1;


   reaction_name[1] = "Random bind";
   reaction_name[2] = "Random unbind";
   reaction_name[3] = "Random step";
   reaction_name[4] = "Overcome ctcf";

   if (seed != 0) Randomize(seed, parm.debug); 
}



/////////////////////////////////////////////
// Simulate an event of Gillespie algorithm 
/////////////////////////////////////////////
bool Extrusion::Event(bool debug=false)
{
   bool ok;
   int r;

   if ( debug ) cerr << to_string(iTime)+") Starting Gillespie event" << endl;

   // Calculate propensities for the different reactions
   CalulatePropensities(debug);
   if ( propensities[0] < SMALL ) 
   {
      exitError = "All propensities are zero";
      return false;
   }

   // Time of next reaction
   tau = log( 1./ DRand() ) / propensities[0];
   if (debug) cerr << "tau = "+to_string(tau) << endl;

   // Choose which reaction
   r = SelectReaction();
   if (debug) cerr << "Selected reaction is r="+to_string(r)+"  : "+reaction_name[r] << endl;
   if ( r < 0 ) return false;

   // Apply chosen reaction
   ok = ApplyReaction(r, debug);

   return ok;
}


/////////////////////////////////////////////
// Bind an extruder to random sites i and i+1 
/////////////////////////////////////////////
bool Extrusion::RandomBind( bool debug=false )
{
   int i = iRand(length-1);
   if (debug) cerr << to_string(iTime)+") Random bind extruder at sites "+to_string(i)+"-"+to_string(i+1) << endl; 
   return AddExtruder(i, i+1, iTime, iTime);        
}


/////////////////////////////////////////////
// Unbind an extruder chosen at random 
/////////////////////////////////////////////
bool Extrusion::RandomUnbind( bool debug=false )
{
   int w = iRand( n_extr_bound );
   int i = extrList[w][0];
   int j = extrList[w][1];
   if (debug) cerr << to_string(iTime)+") Random unbind extruder from sites "+to_string(i)+"-"+to_string(j)
                       +" (w="+to_string(w)+")" << endl; 
   return RemoveExtruder(i, j);
}


/////////////////////////////////////////////
// Step extruder forward in sites
// the argument ctcf_cross indicates if only CTCF sites or only nonCTCF sites are considered
/////////////////////////////////////////////
bool Extrusion::RandomStepForward(  bool ctcf_cross,  bool debug=false )
{
   int i, j, iTimeI, iTimeJ, w, cnt=0, dir;
   bool go_ahead;

   do
   {
       w = iRand( n_extr_bound );
       i = extrList[w][0];
       j = extrList[w][1];
       iTimeI = extrList[w][2];
       iTimeJ = extrList[w][3];
       dir = iRand(2);  // 0=move i, 1=move j

       if (debug) cerr << " testing extruder step from "+to_string(i)+"-"+to_string(j)+" (w="+to_string(w)+
                           ") direction="+to_string(dir) << endl;
     
       // be sure that i<j 
       if ( j < i) 
       { 
          i = extrList[w][1];
          j = extrList[w][0];
       }       

       // if it has reached the ends then unbinds
       if (i == 0 || j == length-1 )
       {
          RemoveExtruder(i, j);
          if (debug) cerr << to_string(iTime)+") Reaches one of the ends and unbinds" << endl;
          return true;
       }
  
       // check if overlaps with something 
       go_ahead = CheckStepOk( w, dir, ctcf_cross, debug );
       

       // make a step if allowed 
       // from i
       if ( dir == 0 && go_ahead )	
       {
           RemoveExtruder(i, j);
           AddExtruder(i-1, j, iTime, iTimeJ);
           if (debug) cerr << to_string(iTime)+") Accepted move to "+to_string(i-1)+"-"+to_string(j) << endl;
           return true;
       }
       // from j
       else if ( dir == 1 && go_ahead )
       {

           RemoveExtruder(i, j);
           AddExtruder(i, j+1, iTimeI, iTime);
           if (debug) cerr << to_string(iTime)+") Accepted move to "+to_string(i)+"-"+to_string(j+1) << endl;
           return true;
       }

       cnt ++;

   } while ( cnt < LARGE );        

  exitError = "Cannot make extruder step";
  return false;
}



/////////////////////////////////////////////
// Read ctcf from file 
/////////////////////////////////////////////
bool Extrusion::ReadCTCF(string fileName)
{
   int i;
   string line;

   ifstream fin(fileName);
   if ( fin.is_open() )
   {
      while( getline(fin, line) )
      {
          i = stoi(line);
          if (i<0 || i>=length)
          {
             exitError = "CTCF site out of range, i = "+line;
             return false;
          }
          ctcf[ i ] = 1;
          nCTCF ++;
      }
   }
   else
   {
      exitError = "Cannot open CTCF file "+fileName;
      return false;
   }

   fin.close();

   return true;
}    

/////////////////////////////////////////////
// Print state to file 
/////////////////////////////////////////////
bool Extrusion::PrintState(string fileName)
{
   ofstream fout(fileName);
   if ( fout.is_open() )
   {
      fout << length << " " << n_extr_bound << " " << n_extr_max << " " << nCTCF << " " << iTime << endl;
      for (int i=0;i<n_extr_bound;i++)
      {
        for (int j=0;j<4;j++)
          fout << extrList[i][j] << " ";
        fout << endl;
      }
      for (int i=0;i<length;i++)
        if ( ctcf[i] > 0 ) fout << i << endl;
   }
   else
   {
      exitError = "Cannot open file for writing state: " + fileName;
      return false;
   }

   fout.close();

   return true;
}

/////////////////////////////////////////////
// Read state from file 
/////////////////////////////////////////////
bool Extrusion::ReadState(string fileName, bool debug=false)
{
   int k;

   // delete existing arrays
   delete map;
   delete[] extrList;
   delete ctcf;
   delete occupiedSites;

   // read from file
   ifstream fin(fileName);
   if ( fin.is_open() )
   {
       fin >> length;
       fin >> n_extr_bound;
       fin >> n_extr_max;
       fin >> nCTCF;
       fin >> iTime;

       if (debug) cerr << "Reading from file "+fileName+" "+to_string(n_extr_bound)+" extrusors." << endl;

       map = AlloArrayInt(length);			// rebuild arrays
       extrList = new int [n_extr_max][4];
       ctcf = new int [length];
       for (int i=0;i<length;i++) ctcf[i] = 0;	
       occupiedSites = new int [length];
       for (int i=0;i<length;i++) occupiedSites[i] = 0;

       for (int i=0;i<n_extr_bound;i++)		// read extruders
        for (int j=0;j<4;j++) 
          fin >> extrList[i][j];
       for (int i=0;i<nCTCF;i++)			// read ctcf
       {
          fin >> k;
          ctcf[k] = 1;
       }
       
   }
   else
   {
      exitError = "Cannot open state file "+fileName;
      return false;

   }

   // reset all arrays
   for (int i=0;i<length;i++)
     for (int j=0;j<length;j++) map[i][j] = 0;

   // fill the arrays
   for (int k=0;k<n_extr_bound;k++) 
   { 
      int i = extrList[k][0];
      int j = extrList[k][1];
      map[i][j] ++;
      map[j][i] ++;
      occupiedSites[i] ++;
      occupiedSites[j] ++;
   }

   if (debug) cerr << "Read with success." << endl;

   fin.close();

   return true;
}

/////////////////////////////////////////////
// Print extrusor map 
/////////////////////////////////////////////
bool Extrusion::PrintMap( string fileName="", bool asList=false, bool onlyExist=false )
{
   ofstream tmp;

   ostream & fout = ( fileName != "" ) ? tmp.open(fileName, ios::out), tmp : cout; 
      
   if ( !fout )
   {
       exitError = "Cannot open file "+fileName+" to write map";
       return false;
   }

   if ( !asList )
   {
      for (int i=0;i<length;i++)
      {
        for (int j=0;j<length;j++)
          fout << setw(3) << map[i][j] ;
        fout << endl;
      }
   }
   else
      for (int i=0;i<length;i++)
        for (int j=i+1;j<length;j++)
           if ( !onlyExist || map[i][j] > 0 )
           {
             fout << setw(6) << i;
             fout << setw(6) << j;
             fout << setw(3) << map[i][j] << endl;
           }

   if ( fileName == "" ) tmp.close();

   return true;

}
/////////////////////////////////////////////
// Evaluate if error occurs
/////////////////////////////////////////////
void Extrusion::CatchError( bool ok )
{
   if ( !ok )
   {
      time_t now = time(0);
      char* dt = ctime(&now);

      cerr << endl;
      cerr << "*** Execution error:" << endl;
      cerr << exitError << endl;
      cerr << dt << endl;

      exit(1);
   }
} 




/////////////////////////////////////////////
/////////////////////////////////////////////
// Private functions
/////////////////////////////////////////////
/////////////////////////////////////////////





/////////////////////////////////////////////
// Allocate n*n matrix 
/////////////////////////////////////////////
int **Extrusion::AlloArrayInt(int n)
{
   int **x = new int*[ n ];
   for (int i=0; i<n; i++) 
   {
     x[i] = new int[ n ];
     for (int j=0;j<n;j++) x[i][j] = 0;      
   }
   return x;
}


/////////////////////////////////////////////
// Create an extruder at sites i, j 
/////////////////////////////////////////////
bool Extrusion::AddExtruder(int i, int j, int iTimeI, int iTimeJ)
{
  map[i][j] ++; 
  map[j][i] ++; 
  extrList[n_extr_bound][0] = i;
  extrList[n_extr_bound][1] = j;
  extrList[n_extr_bound][2] = iTimeI;
  extrList[n_extr_bound][3] = iTimeJ;
  occupiedSites[i] ++;
  occupiedSites[j] ++;
  n_extr_bound ++;

  if ( n_extr_bound >= n_extr_max ) 
  {
     exitError = "nEntrMax too small.";
     return false;
  }

  // set output
  add_link = true;
  add_link_i = i;
  add_link_j = j;

  return true;
}

/////////////////////////////////////////////
// Destroy an extruder at sites i, j 
/////////////////////////////////////////////
bool Extrusion::RemoveExtruder(int i, int j)
{
  if ( map[i][j]==0) 
  {
     exitError = "Trying to remove extruder that is not there (i="+to_string(i)+", j="+to_string(j)+")";
     return false;
  }

  map[i][j] --;
  map[j][i] --;
  occupiedSites[i] --;
  occupiedSites[j] --;

  for (int n=0;n<n_extr_bound;n++)
     if ( (extrList[n][0] == i && extrList[n][1] == j) || (extrList[n][0] == j && extrList[n][1] == i) )
     {
        for (int k=0;k<4;k++)
           extrList[n][k] = extrList[n_extr_bound-1][k]; 
        break;
     }

  n_extr_bound --;

  // set output
  delete_link = true;
  delete_link_i = i;
  delete_link_j = j;

  return true;
}


/////////////////////////////////////////////
// Random number in [0,n) 
/////////////////////////////////////////////
int Extrusion::iRand(int n)
{
   return rand() % n;
}

/////////////////////////////////////////////
// Random number double in [0,1)
/////////////////////////////////////////////
double Extrusion::DRand( void )
{
   double r = (double) rand() / RAND_MAX;
   return r ;
}


/////////////////////////////////////////////
// logical xor 
/////////////////////////////////////////////
bool Extrusion::LogicalXOR(bool a, bool b)
{
   return a != b;
}

/////////////////////////////////////////////
// randomize random numbers (-1 for seed from time) 
/////////////////////////////////////////////
void Extrusion::Randomize(int seed, bool debug=false)
{
   if (seed<0) seed = time(NULL) % 29799493;
   if ( debug ) cerr << "Random seed = "+to_string(seed) << endl;
   for (int i=0;i<seed;i++) iRand(seed%100);
}


/////////////////////////////////////////////
// calculate propensities for Gillespie algorithm 
/////////////////////////////////////////////
bool Extrusion::CalulatePropensities( bool debug=false )
{
   int n_extr_free, n_steppable=0, n_cross_ctcf=0;

   for (int i=0;i<NREACT+1;i++) propensities[i] = 0.;

   // 1 - random binding
   if ( n_extr_tot > 0 ) n_extr_free = n_extr_tot - n_extr_bound;
   else n_extr_free = 1;
   propensities[1] = k_binding * n_extr_free;

   // 2 - random unbinding
   propensities[2] = k_unbinding * n_extr_bound;

   // 3 - stepping (no ctcf)
   if ( k_step > 0 )
   {
      for (int w=0;w<n_extr_bound;w++)
        for (int dir=0;dir<2;dir++) 
          if ( CheckStepOk(w, dir, false,  false ) ) n_steppable += 1;
      propensities[3] = k_step * n_steppable;
   }

   // 4 - crossing ctcf
   if ( k_cross_ctcf > 0)
   {
      for (int w=0;w<n_extr_bound;w++)
        for (int dir=0;dir<2;dir++) 
          if ( CheckStepOk(w, dir, true,  false ) ) n_cross_ctcf += 1;
      propensities[4] = k_cross_ctcf * n_cross_ctcf;
   }

   for (int i=0;i<NREACT;i++) propensities[0] += propensities[i];

   if (debug)
   {
      cerr << "Propensities:" << endl;
      for (int i=1;i<NREACT;i++) cerr << "a["+to_string(i)+"]="+to_string(propensities[i])+"  - "+reaction_name[i] << endl;
      cerr << "=> a[0]="+to_string(propensities[0]) << endl;
   }

   return true;
}


/////////////////////////////////////////////
// check if suggested step clashes with another extrusor and if is on ctcf 
/////////////////////////////////////////////
bool Extrusion::CheckStepOk( int w, int dir, bool ctcf_cross, bool debug=false )
{
   int i = extrList[w][0];
   int j = extrList[w][1];
   int iTimeI = extrList[w][2];
   int iTimeJ = extrList[w][3];

   if (debug) cerr << " testing extruder step from "+to_string(i)+"-"+to_string(j)+" (w="+to_string(w)+
                           ") direction="+to_string(dir) << endl;
     
   
   // check if it is overcoming another extrusor
   if ( !allow_overcome )
   {
      if ( dir == 0 )
      {
        for (int k;k<n_extr_bound;k++) 	// if there is an extrusor in i from more time, skip.
           if ( (extrList[k][0] == i && extrList[k][2] < iTimeI ) || 
                    (extrList[k][1] == i && extrList[k][3] < iTimeI ) ) 
           {
              if (debug) cerr << "  step is stopped by overlap with w="+to_string(k)+" ("+
                                     to_string(extrList[k][0])+"-"+to_string(extrList[k][1])+")" << endl;
              return false;
           }
      }
      else if ( dir == 1 )
      {
        for (int k;k<n_extr_bound;k++) 	// if there is an extrusor in j from more time, skip.
           if ( (extrList[k][0] == j && extrList[k][2] < iTimeI ) || 
                    (extrList[k][1] == j && extrList[k][3] < iTimeI ) ) 
           {
              if (debug) cerr << "  step is stopped by overlap with w="+to_string(k)+" ("+
                                  to_string(extrList[k][0])+"-"+to_string(extrList[k][1])+")" << endl;
              return false;
           }
      }
   }

   // check if meeting ctcf condition of the function argument
   // of i
       if ( dir == 0 && LogicalXOR( ctcf[i-1] == 0, ctcf_cross )  )	
           return true;
       // of j
       else if ( dir == 1 && LogicalXOR( ctcf[j+1] == 0,  ctcf_cross )  )
           return true;

   return false;

}



/////////////////////////////////////////////
// Select Gillespie reaction 
/////////////////////////////////////////////
int Extrusion::SelectReaction( void )
{
   double r,aSum=0.;

   r = propensities[0] * DRand();

   for (int i=1;i<NREACT;i++)
   {
      aSum += propensities[i];
      if ( r < aSum ) return i;
   }

  exitError = "All propensities are zero."; 
  return -1;
}

/////////////////////////////////////////////
// Apply Gillespie reaction 
/////////////////////////////////////////////
bool Extrusion::ApplyReaction(int r, bool debug=false)
{
   bool ok;

   switch (r)
   {
      case 1:	// random binding
      {
        ok = RandomBind( debug );
        break;
      }

      case 2:	// random unbinding
      {
        ok = RandomUnbind( debug );
        break;
      }

      case 3:	// step far from ctcf
      {
        ok = RandomStepForward( false, debug );
        break;
      }
      case 4:	// overcome ctcf
      {
        ok = RandomStepForward( true, debug );
        break;
      }
   }
   
   if ( ok ) iTime ++; 
   return ok;
}




