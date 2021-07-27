
OBJS = GTMatrix.o GeneGlobDriver.o  TDIMatrix.o GeneTDIGlobalDriver_Main.o 
DEPS = GTMatrix.h  GeneGlobDriver.h TDIMatrix.h

%.o: %.cpp $(DEPS)
	g++  -O3 -c $< -o $@


TCI_GD : $(OBJS)
	g++  -O3  $(OBJS) -o $@ -lm

clean:
	rm -f *.o TCI_GD