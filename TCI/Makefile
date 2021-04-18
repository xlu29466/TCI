
OBJS = GTMatrix.o  TDIC.o TDIMatrix.o TDIC_Main.o 
DEPS = GTMatrix.h  TDIC.h TDIMatrix.h

%.o: %.cpp $(DEPS)
	g++  -O3 -c $< -o $@


TCI : $(OBJS)
	g++  -O3  $(OBJS) -o $@ -lm

clean:
	rm -f *.o TCI