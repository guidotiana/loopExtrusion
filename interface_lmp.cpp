#include "interface_lmp.h"
#include "update.h"

Interface_lmp::Interface_lmp(int argc, char **argv)
{
   cout << "Opening interface with LAMMPS..." << endl;
   cout << "" << endl;
   initiate_lmp(argc, argv);
}

void Interface_lmp::initiate_lmp(int argc, char **argv)
{
   cout << "Initializing LAMMPS..." << endl;
   cout << "" << endl;

   MPI_Init(&argc,&argv);
  
  /*
  if (argc != 3) {
    printf("Syntax: simpleCC P in.lammps\n");
    exit(1);
  }*/

  int me,nprocs;
  MPI_Comm_rank(MPI_COMM_WORLD,&me);
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);

  //int lammps;
  lammps = MPI_UNDEFINED;

  // open LAMMPS input script
  
  FILE *fp;
  if (me == 0) {
     fp = fopen(argv[2],"r");
     if (fp == NULL) {
        printf("ERROR: Could not open LAMMPS input script\n");
     }
  }
  
  // run the input script thru LAMMPS one line at a time until end-of-file
  // driver proc 0 reads a line, Bcasts it to all procs
  // (could just send it to proc 0 of comm_lammps and let it Bcast)
  // all LAMMPS procs call input->one() on the line
  
  //LAMMPS_NS::LAMMPS *lmp = NULL;
  lmp = new LAMMPS_NS::LAMMPS(0,NULL,MPI_COMM_WORLD);
  
  //Turn off screen output
  lmp->screen = NULL;  
  //Turn off log output
  lmp->logfile = NULL;  

  isFirstRun = true;

  //lammps_command(lmp, "screen none");
  int n;
  char line[1024];
  while (1) {
      if (fgets(line,1024,fp) == NULL) n = 0;
      else n = strlen(line) + 1;
      if (n == 0) fclose(fp);
      if (n == 0) break;
      lammps_command(lmp,line);
  }
 
}

void Interface_lmp::load_bond(int bond_type, int new_id1, int new_id2)
{
   //create bond
   ostringstream line;
   line << "create_bonds single/bond " << bond_type << " " << new_id1 << " " << new_id2;
   string MyString  = line.str();
   lammps_command(lmp, MyString.c_str());
   line.str("");
   line.clear();
}

void Interface_lmp::unload_bond(int bond_type, int old_id1, int old_id2)
{
   stringstream line;
   line << "group to_remove id " << old_id1 << " " << old_id2;
   string MyString = line.str();

   //create group with atoms whose bond must be removed 
   lammps_command(lmp, MyString.c_str());
   line.str("");
   line.clear();
   
   //delete bond
   line << "delete_bonds to_remove bond " << bond_type << " remove";
   MyString = line.str();
   lammps_command(lmp, MyString.c_str());
   line.str("");
   line.clear();
   
   //delete group 
   lammps_command(lmp,"group to_remove delete");
}

void Interface_lmp::update_bonds(int bond_type, bool add_link, bool delete_link, int add_link_i, int add_link_j, int delete_link_i, int delete_link_j)
{
   if (add_link) load_bond(bond_type, add_link_i+1, add_link_j+1);
   if (delete_link) unload_bond(bond_type, delete_link_i+1, delete_link_j+1);
} 

void Interface_lmp::run_dynamics(int steps)
{  
   lmp->update->restrict_output = 0;

   if (isFirstRun = true) {
      stringstream line;
      line << "run " << steps;
      string MyString = line.str();

      //launch lammps dynamics
      lammps_command(lmp, MyString.c_str());
      line.str("");
      line.clear();
   
      isFirstRun = false;
   }
   else {
      stringstream line;
      line << "run " << steps << " pre no post no";
      string MyString = line.str();

      //launch lammps dynamics
      lammps_command(lmp, MyString.c_str());
      line.str("");
      line.clear(); 
   }

}

void Interface_lmp::write_data(string str)
{
   stringstream line;
   line << str;
   string MyString = line.str();

   //write data file
   lammps_command(lmp, MyString.c_str());
   line.str("");
   line.clear();  
}

void Interface_lmp::minimize()
{  
   //don't dump/output minimization data
   lmp->update->restrict_output = 1;

   //getting dynamics time
   int * ntimestep_ptr;
   int ntimestep;
   ntimestep_ptr = (int *) lammps_extract_global(lmp, "ntimestep"); 
   ntimestep = *ntimestep_ptr; 
  
   //minimize
   stringstream line; 
   line << "minimize 1e-5 1e-5 1000 1000";
   string MyString = line.str();
   lammps_command(lmp, MyString.c_str());
   line.str("");
   line.clear();
 
   //don't count minimization steps as dynamics steps
   line << "reset_timestep " << ntimestep;
   MyString = line.str();
   lammps_command(lmp, MyString.c_str());
   line.str("");
   line.clear(); 
} 

void Interface_lmp::close_lmp()
{
  //closing lammps and MPI
  delete lmp;

  MPI_Finalize();
  
}
