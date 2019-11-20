CC = mpicc
CXX = mpic++
CFLAGS = -Wall -Wextra -g
NP = 2
FILE = ../data/test_sm.mtx

all:	main

main:	mmio.o matrix.o graph.o main.o 
	$(CXX) $(CFLAGS) -o $@ $^ -lrt


mmio.o:	mmio.c
	$(CC) -c mmio.c

%.o:	%.cc
	$(CXX) $(CFLAGS) -c $<


run: clean all
	clear
	# ./$(BIN)/$(EXECUTABLE)
	mpirun -np $(NP) ./main $(FILE)

debug: clean all
	clear
	mpirun -np $(NP) xterm -hold -e gdb -ex run --args ./main $(FILE)
	
run_again:
	clear
	# ./$(BIN)/$(EXECUTABLE)
	mpirun -np $(NP) main $(FILE)

clean:
	rm -f *o
	rm -f lu
