CC=/home/jtmeng/software/mpich-3.1.3/install/bin/mpic++

all	:swap stats
swap    : graph.o kmerGraph.o sequence.o mympi.o
	$(CC) -O3 graph.o kmerGraph.o sequence.o mympi.o  -o swap
stats   :stats.cpp
	$(CC) stats.cpp -o stats	

graph.o : graph.cpp kmerGraph.h mympi.h sequence.h
	$(CC) -Wno-deprecated -c -O3 graph.cpp -o graph.o
kmerGraph.o : kmerGraph.cpp kmerGraph.h mympi.h sequence.h
	$(CC) -Wno-deprecated -c  -O3 kmerGraph.cpp -o kmerGraph.o
sequence.o : sequence.cpp sequence.h mympi.h
	$(CC) -Wno-deprecated -c  -O3 sequence.cpp -o sequence.o
mympi.o : mympi.cpp mympi.h
	$(CC) -Wno-deprecated -c  -O3 mympi.cpp -o mympi.o

clean :
	rm -f stats swap graph.o kmerGraph.o sequence.o mympi.o

#	$(CC) -c -Wno-deprecated -O2 graph.cpp -o graph.o
