CC = g++
CFLAGS = -Wall -std=c++11

main: main.o random.o
	$(CC) $(CFLAGS) -o main main.o random.o
main.o: main.cpp
	$(CC) $(CFLAGS) -c main.cpp
random.o: random.cpp random.h
	$(CC) $(CFLAGS) -c random.cpp
clean:
	rm -f main main.o