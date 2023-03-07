CPP = mpicxx
CFLAGS = -I. -std=c++0x -I/home/edoardo/prog/lammps-29Sep2021/src -g 
LFLAGS = -lm -L/home/edoardo/prog/lammps-29Sep2021/build -llammps 
DEPS = extrusion.h parameters.h interface_lmp.h
OBJ = loopExtrusion.o extrusion.o parameters.o interface_lmp.o

%.o:  %.cpp $(DEPS)
	$(CPP) -c -o $@ $< $(CFLAGS)

loopExtrusion: $(OBJ)
	$(CPP) -o $@ $(OBJ) $(LFLAGS)

clean:
	rm -f *.o loopExtrusion 
