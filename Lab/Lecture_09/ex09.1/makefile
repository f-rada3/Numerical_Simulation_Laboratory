CC = g++
CFLAGS = -O3 --std=c++14 -I/opt/homebrew/include
AFLAGS = -L/opt/homebrew/lib -larmadillo
all: main.exe

main.exe: main.o random.o population.o individual.o
	$(CC) random.o population.o individual.o main.o -o main.exe $(LIBS)
main.o: main.cpp population.h random.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
individual.o : individual.cpp individual.h random.h
	$(CC) -c individual.cpp -o individual.o $(CFLAGS)
population.o : population.cpp population.h individual.h random.h
	$(CC) -c population.cpp -o population.o $(CFLAGS)
clean :
	rm -f *.o *.exe seed.out