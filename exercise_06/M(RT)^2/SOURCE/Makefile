CC = g++
CFLAGS = -O3 --std=c++11
AFLAGS = -larmadillo
FILENAME = NSL_SIMULATOR
EXE = simulator
all: $(EXE).exe
$(EXE).exe : system.o $(FILENAME).o particle.o random.o
	$(CC) system.o particle.o random.o $(FILENAME).o -o $(EXE).exe $(AFLAGS)
$(FILENAME).o : $(FILENAME).cpp system.h particle.h
	$(CC) -c $(FILENAME).cpp -o $(FILENAME).o $(CFLAGS)
system.o : system.cpp system.h random.h particle.h
	$(CC) -c system.cpp -o system.o $(CFLAGS)
particle.o : particle.cpp particle.h
	$(CC) -c particle.cpp -o particle.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o $(EXE).exe
remove:
	rm ../OUTPUT/*.*
	rm ../OUTPUT/CONFIG/*.*
