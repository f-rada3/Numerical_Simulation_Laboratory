# =================================================================
# Makefile per il TSP con MPI + Armadillo
# Assicurarsi che OpenMPI sia installato (brew install open-mpi)
# =================================================================

# wrapper MPI per compilare e linkare
CC       = mpicxx

# flags per compilazione
CXXFLAGS = -g -O0 -std=c++14 -I/opt/homebrew/include

# flags per il link (Armadiilo)
LDFLAGS  = -L/opt/homebrew/lib -larmadillo

# target ed oggetti
TARGET   = main.exe
OBJS     = main.o random.o population.o individual.o

.PHONY: all clean run

all: $(TARGET)

# link finale
$(TARGET): $(OBJS)
	$(CC) $(OBJS) -o $(TARGET) $(LDFLAGS)

# regola generica per .cpp → .o
%.o: %.cpp
	$(CC) -c $< -o $@ $(CXXFLAGS)

clean:
	rm -f $(OBJS) $(TARGET)

# regola per lanciare con 8 processi
run: $(TARGET)
	mpirun -np 8 ./$(TARGET)