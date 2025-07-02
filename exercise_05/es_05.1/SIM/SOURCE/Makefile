CC = g++
CFLAGS = -O3 --std=c++11
all: main.exe
main.exe : hydrogen_atom.o main.o random.o
	$(CC) hydrogen_atom.o random.o main.o -o main.exe 
main.o : main.cpp hydrogen_atom.h 
	$(CC) -c main.cpp -o main.o $(CFLAGS)
hydrogen_atom.o : hydrogen_atom.cpp hydrogen_atom.h random.h 
	$(CC) -c hydrogen_atom.cpp -o hydrogen_atom.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o main.exe
remove:
	rm ../OUTPUT/*.*
