CC = g++
CFLAGS = -O3 --std=c++11
AFLAGS = -larmadillo
all: main.exe
main.exe : population.o main.o individual.o random.o
	$(CC) population.o individual.o random.o main.o -o main.exe $(AFLAGS)
main.o : main.cpp population.h individual.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)
population.o : population.cpp population.h random.h individual.h
	$(CC) -c population.cpp -o population.o $(CFLAGS)
individual.o : individual.cpp individual.h
	$(CC) -c individual.cpp -o individual.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o main.exe
remove:
	rm ../OUTPUT/*.*
	rm ../OUTPUT/CONFIG/*.*
