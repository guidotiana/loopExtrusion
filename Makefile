CPP = c++
CFLAGS = -I.
LFLAGS = -lm
DEPS = extrusion.h parameters.h
OBJ = loopExtrusion.o extrusion.o parameters.o

%.o:  %.cpp $(DEPS)
	$(CPP) -c -o $@ $< $(CFLAGS)

loopExtrusion: $(OBJ)
	$(CPP) -o $@ $(OBJ) $(LFLAGS)

clean:
	rm -f *.o loopExtrusion 
