CC = mpicc
CXX = mpic++
CFLAGS = -Wall -Wextra -O3 -g
NP = 1
FILE = data/test_sm.mtx

all:	main

main:	main.o mmio.o matrix.o sparse.o
	$(CXX) $(CFLAGS) -o $@ $^ -lrt


mmio.o:	mmio.c
	$(CC) -c mmio.c

%.o:	%.cc
	$(CXX) $(CFLAGS) -c $<


run: clean all
	clear
	# ./$(BIN)/$(EXECUTABLE)
	mpirun -np $(NP) main $(FILE)

run_again:
	clear
	# ./$(BIN)/$(EXECUTABLE)
	mpirun -np $(NP) main $(FILE)

clean:
	rm -f *o
	rm -f lu
