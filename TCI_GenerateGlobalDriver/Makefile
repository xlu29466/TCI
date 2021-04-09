
OBJS = GTMatrix.o PanCanGTMatrix.o GeneGlobDriver.o GenePanCanGlobDriver.o TDIMatrix.o GeneTDIGlobalDriver_Main.o 
DEPS = GTMatrix.h PanCanGTMatrix.h GeneGlobDriver.h GenePanCanGlobDriver.h TDIMatrix.h

%.o: %.cpp $(DEPS)
	g++  -O3 -c $< -o $@


TCI_GD : $(OBJS)
	g++  -O3  $(OBJS) -o $@ -lm

clean:
	rm -f *.o TCI_GD