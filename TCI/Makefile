
OBJS = GTMatrix.o PanCanGTMatrix.o PanCanTDIC.o TDIC.o TDIMatrix.o TDIC_Main.o 
DEPS = GTMatrix.h PanCanGTMatrix.h PanCanTDIC.h TDIC.h TDIMatrix.h

%.o: %.cpp $(DEPS)
	g++  -O3 -c $< -o $@


TCI : $(OBJS)
	g++  -O3  $(OBJS) -o $@ -lm

clean:
	rm -f *.o TCI