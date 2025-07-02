CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o function.o
	$(CC) random.o function.o main.o -o main.exe

main.o : main.cpp function.h random.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)

function.o : function.cpp function.h
	$(CC) -c function.cpp -o function.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

clean :
	rm *.o main.exe
