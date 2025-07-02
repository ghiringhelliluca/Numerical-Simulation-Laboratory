CC = g++
CFLAGS = -O3 --std=c++11
AFLAGS = -larmadillo
all: equilibrator.exe
equilibrator.exe : system.o NSL_EQUILIBRATOR.o particle.o random.o
	$(CC) system.o particle.o random.o NSL_EQUILIBRATOR.o -o equilibrator.exe $(AFLAGS)
NSL_EQUILIBRATOR.o : NSL_EQUILIBRATOR.cpp system.h particle.h
	$(CC) -c NSL_EQUILIBRATOR.cpp -o NSL_EQUILIBRATOR.o $(CFLAGS)
system.o : system.cpp system.h random.h particle.h
	$(CC) -c system.cpp -o system.o $(CFLAGS)
particle.o : particle.cpp particle.h
	$(CC) -c particle.cpp -o particle.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o equilibrator.exe
remove:
	rm ../OUTPUT/*.*
	rm ../OUTPUT/CONFIG/*.*
