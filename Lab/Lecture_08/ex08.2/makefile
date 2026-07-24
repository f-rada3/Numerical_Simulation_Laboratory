CC = g++
CFLAGS = -Wall -std=c++11

all: main.exe PostProcessing.exe

main.exe: main.o random.o
	$(CC) $(CFLAGS) -o main.exe main.o random.o

PostProcessing.exe: PostProcessing.o random.o
	$(CC) $(CFLAGS) -o PostProcessing.exe PostProcessing.o random.o

main.o: main.cpp random.h
	$(CC) $(CFLAGS) -c main.cpp

PostProcessing.o: PostProcessing.cpp random.h
	$(CC) $(CFLAGS) -c PostProcessing.cpp

random.o: random.cpp random.h
	$(CC) $(CFLAGS) -c random.cpp

clean:
	rm -f *.o *.exe seed.out
