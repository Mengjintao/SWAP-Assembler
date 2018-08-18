#ifndef _KMERGRAPH_H_
#define _KMERGRAPH_H_

#include <map>
#include <string>
#include <ext/hash_map>
#include <unistd.h>
#include "sequence.h"
#include "mympi.h"

using namespace std;

class kmerGraph
{
public:
	unsigned long long *kmers;
	unsigned char *arcs;
	unsigned long long size;
	unsigned long long read_pos;
	double commtime, commtimetot;
	double cuttime, storagetime;
	MPI_Datatype commType;

	unsigned long long readRound;
	char *readBuf;
	unsigned long long readStart;
	unsigned long long readEnd;
	unsigned long long readLen;
	
	unsigned long long reverseComplement(unsigned long long kmerDescriptor, parameter *parameters);
    	static unsigned long long stringToLongLong(const char *buf, int start, int end, parameter *parameters);
	static string longLongToString(unsigned long long a, parameter *parameters);
	unsigned long long getProcsID(unsigned long long kmerID, int hashLength, MPIEnviroment *MPIcontrol);
	int arcPos(unsigned long long &A, unsigned long long &B, char directA, char directB, int hashLength);

	kmerGraph()
	{
		kmers    = NULL;
		arcs     = NULL;
		size     = 0;
		read_pos = 0;
		commtime = 0;
		commtimetot = 0;
		cuttime = 0;
		storagetime = 0;	
		
		readRound = 0;
		readBuf  = NULL;
		readStart = 0;
		readEnd   = 0;

		MPI_Type_contiguous(KMER_COMM_TYPE_LEN, MPI_UNSIGNED_LONG_LONG, &commType);
		MPI_Type_commit(&commType);
	}
	
	~kmerGraph()
	{
		if(kmers!=NULL)		delete(kmers);
		if(arcs!=NULL)		delete(arcs);
		if(readBuf!=NULL)	delete(readBuf);
		kmers = NULL;
		arcs = NULL;
		readBuf = NULL;
	}

public:
	int  constructKmerGraph(parameter *parameters, MPIEnviroment *MPIcontrol);
	void printKmerGraph(parameter *parameters, MPIEnviroment *MPIcontrol);
	void distributeKmerGraph(parameter *parameters, MPIEnviroment *MPIcontrol);
};

#endif
