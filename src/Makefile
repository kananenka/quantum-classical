COMP=icpc
FFLAGS=-L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl -lfftw3 
LIBS=-L/usr/local/lib
INCLUDES=#-I${MKLROOT}/include
OBJECTS=main.o read_write.o move.o inter.o util.o energy.o pes.o qs.o\
        cstr.o

run: $(OBJECTS)
	$(COMP) -o run $(INCLUDES) $(LIBS) $(OBJECTS) $(FFLAGS)
	cp run ../bin

main.o: main.cpp 
	$(COMP) -c main.cpp 

read_write.o: read_write.cpp
	$(COMP) -c read_write.cpp

move.o: move.cpp
	$(COMP) -c move.cpp

inter.o: inter.cpp
	$(COMP) -c inter.cpp

util.o: util.cpp
	$(COMP) -c util.cpp

energy.o: energy.cpp
	$(COMP) -c energy.cpp

pes.o: pes.cpp
	$(COMP) -c pes.cpp

qs.o: qs.cpp
	$(COMP) -c qs.cpp

cstr.o: cstr.cpp
	$(COMP) -c cstr.cpp

clean:
	rm *.o run

