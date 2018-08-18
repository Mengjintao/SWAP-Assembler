#include "sequence.h"
#include "mympi.h"
#include "kmerGraph.h"
#include <pthread.h>

#include <map>
#include <string>
#include <ext/hash_map>
#include <unistd.h>

#define EndTunnel       2
#define DataTunnel      1
#define ServiceTunnel   100000
#define ComputingTunnel 0 
//using namespace __gnu_cxx;

#define LDEBUG 0 

using namespace std;
using namespace __gnu_cxx;

namespace __gnu_cxx
{
        template<> struct hash< std::string > {
                size_t operator()( const std::string& x ) const {
                        return hash< const char* >()( x.c_str() );
                }
        };
        template<> struct hash< unsigned long long > {
                size_t operator()( const unsigned long long & x ) const {
                        return (x >> 32L) ^ hash< int >()( x & 0xFFFFFFFF );
                }
        };
}

class arc
{
public:
	unsigned char multiplicity[8];	
};

class node
{
    public:
	string arcs[8];				//0-3 arc, 4-7 twinarc
	unsigned char multiplicity[8];	
	unsigned long long nodeID;
	pthread_mutex_t lockFlag;
	char deleteFlag;	
	
	node()
	{
		nodeID = 0;
		for(int i=0;i<8;i++)	arcs[i].clear();
		for(int i=0;i<8;i++)	multiplicity[i] = 0;
	}
	node(const node& C)
	{
		nodeID = C.nodeID;
		for(int i=0;i<8;i++)	arcs[i] = C.arcs[i];
		for(int i=0;i<8;i++)	multiplicity[i] = C.multiplicity[i];
		deleteFlag = C.deleteFlag;
		lockFlag = C.lockFlag;	
	}
	~node();
		
	void init(unsigned long long ID, arc kmoleculeArc, parameter *parameters);
	void Union(unsigned char arc, parameter *parameters);
	const char *getLeftEdge(unsigned char& curMulti);
	const char *getRightEdge(unsigned char& curMulti);
	
	int  getNodeID(int k, unsigned long long &ID, parameter *parameters);
	int  getReverseEdge(int k, unsigned long long &ID, int &IDindex, parameter *parameters);
	int  getLeftNodeID(unsigned long long& ID,  parameter *parameters);
	int  getRightNodeID(unsigned long long& ID, parameter *parameters);
	int  getSrcNodeEdge(unsigned long long ID, int srcDirection, string srcEdge, int srcMulti, int direction, parameter *parameters);
	void extendEdge(unsigned long long neighbourID, string edge, parameter *parameters);
	void cutoff(int threshold);
};

void node::init(unsigned long long kmerID, arc kmoleculeArc, parameter *parameters)
{
    nodeID = kmerID;
    deleteFlag = 1;

    for(int i=0;i<8;i++)	multiplicity[i] = 0;	

    for(int i=0;i<8;i++)
    {
		arcs[i].clear();
		if(kmoleculeArc.multiplicity[i])
	    	{
			arcs[i]+=parameters->nucleotideArray[i%4];
			multiplicity[i]=kmoleculeArc.multiplicity[i];
		}
    }
    if(pthread_mutex_init( &this->lockFlag,	NULL) != 0)
    {
		printf("pthread_mutex_init error\n");
		exit(0);
    }
}

node::~node()
{
    for(int i=0;i<8;i++)	arcs[i].clear();
    pthread_mutex_destroy(&this->lockFlag);	
}

void node::cutoff(int threshold)
{
	for(int i=0;i<8;i++)
	{
		if(multiplicity[i]<threshold)	
		{	
			arcs[i].clear();
			multiplicity[i] = 0;
		}
	}	
}

void node::Union(unsigned char arc, parameter *parameters)
{
    for(int i=0;i<8;i++)
    {
	if( arc & (1<<i) )
	{
		if(arcs[i].empty())		//add arc, if the arc is not exist.
    			arcs[i] += parameters->nucleotideArray[i%4];	
		multiplicity[i]++;		//add multiplicity.
	}
    }
}

const char *node::getLeftEdge(unsigned char& curMultiplicity)
{
    int count =  0;
    int Index = -1;
    for(int i=4;i<8;i++)
    {
	if(arcs[i].length()!=0) 	{Index=i, count++;}
    }			
    if(count!=1)	return NULL;
//	string ret = string("+") + arcs[Index];
	curMultiplicity = multiplicity[Index]; 
    return arcs[Index].c_str();
}

const char *node::getRightEdge(unsigned char& curMultiplicity)
{
    int count = 0;
    int Index = -1;
    for(int i=0;i<4;i++)
    {
		if(arcs[i].length() != 0) 	{Index=i, count++;}
    }			
    if(count!=1)	return NULL;

	curMultiplicity = multiplicity[Index];
    return arcs[Index].c_str();
}


int node::getNodeID(int k, unsigned long long &ID, parameter *parameters)
{
    int ret;

    string descriptor = kmerGraph::longLongToString(nodeID, parameters);
    
    if(k>=4)	//this is the negative side of kmer	
    {	
	reverse(descriptor.begin(), descriptor.end());
    	for(int i=0;i<descriptor.length();i++)	
	    descriptor[i] = parameters->nucleotideReverse[descriptor[i]];
    }	
   
    descriptor += arcs[k];

    int pos = descriptor.length()-parameters->hashLength;

    string leftDescriptor = descriptor.substr(pos);
    string reverseStr = leftDescriptor;
    reverse(reverseStr.begin(), reverseStr.end());
    for(int i=0;i<reverseStr.length();i++)
	reverseStr[i] = parameters->nucleotideReverse[reverseStr[i]];

    if(leftDescriptor < reverseStr)
    {
	ret = -1;
	leftDescriptor = reverseStr;
    }
    else 	ret = 1;

    ID = kmerGraph::stringToLongLong(leftDescriptor.c_str(), 0, leftDescriptor.length(),parameters);		
    return ret; 
}

int node::getReverseEdge(int k, unsigned long long &ID, int &IDindex, parameter *parameters)
{
    int ret;

    string descriptor = kmerGraph::longLongToString(nodeID, parameters);
    
    if(k>=4)	//this is the negative side of kmer	
    {	
	reverse(descriptor.begin(), descriptor.end());
    	for(int i=0;i<descriptor.length();i++)	
	    descriptor[i] = parameters->nucleotideReverse[descriptor[i]];
    }	
   
    descriptor += arcs[k];

    int pos = descriptor.length()-parameters->hashLength;

    int strIndex = 3 - parameters->nucleotideValue[ descriptor[pos-1] ];

    string leftDescriptor = descriptor.substr(pos);
    string reverseStr = leftDescriptor;
    reverse(reverseStr.begin(), reverseStr.end());
    for(int i=0;i<reverseStr.length();i++)
	reverseStr[i] = parameters->nucleotideReverse[reverseStr[i]];

    if(leftDescriptor < reverseStr)
    {
	ret = -1;
	leftDescriptor = reverseStr;
	IDindex = strIndex;
    }
    else
    {
	ret = 1;
	IDindex = 4+strIndex; 	
    }

    ID = kmerGraph::stringToLongLong(leftDescriptor.c_str(), 0, leftDescriptor.length(),parameters);		
    return ret;
 
}


int node::getLeftNodeID(unsigned long long& ID, parameter *parameters)
{
    int count = 0, index=-1, ret;
    for(int i=4; i<8; i++)	
    {
		if(arcs[i].length()!=0)	{ index=i; count++;}
    }

    if(count!=1)	return count;

    return this->getNodeID(index, ID, parameters); 
}

int node::getRightNodeID(unsigned long long& ID, parameter *parameters)
{
    int count=0, index, ret;
    for(int i=0;i<4;i++)	
	if(arcs[i].length()!=0)	{index=i, count++;}
    if(count!=1)	return 0;
    return getNodeID(index, ID, parameters);
}

int node::getSrcNodeEdge(unsigned long long srcNodeID, int srcDirection, string srcEdge, int srcMulti, int direction, parameter *parameters)
{
	int startIndex, endIndex, index=-1, ret;
	if(direction <0 && direction>2)	
	{
		printf("error in get SrcNodeEdge\n");	
		exit(0);	
	}

	if(direction == 0)		startIndex=4, endIndex=8;	//negative side edge
	if(direction == 1)		startIndex=0, endIndex=4;	//positive side edge
	if(direction == 2)		startIndex=0, endIndex=8;	//two side edge	
	
	for(int i=startIndex; i<endIndex; i++)
	{
		string descriptor    = kmerGraph::longLongToString(this->nodeID, parameters);
      		if(i>3) //negative side k-mers
        	{
                        reverse(descriptor.begin(), descriptor.end());

                        for(int j=0;j<descriptor.length();j++)
                                descriptor[j] = parameters->nucleotideReverse[descriptor[j]];
                }

                descriptor += this->arcs[i];

                int pos = descriptor.length()-parameters->hashLength;

                string retDescriptor = descriptor.substr(pos);

                string reverseStr = retDescriptor;
                reverse(reverseStr.begin(), reverseStr.end());
                for(int j=0;j<reverseStr.length();j++)
                	reverseStr[j] = parameters->nucleotideReverse[reverseStr[j]];

                if(retDescriptor < reverseStr)
                {
                        ret = 0;		//negative side
                        retDescriptor = reverseStr;
                }
                else    ret = 1;		//positive side
                assert(retDescriptor.length() == parameters->hashLength);
                unsigned long long ID = kmerGraph::stringToLongLong(retDescriptor.c_str(),0,retDescriptor.length(),parameters);

		string srcArc = descriptor;
		reverse(srcArc.begin(), srcArc.end());
		for(int j=0;j<srcArc.length();j++)
			srcArc[j] = parameters->nucleotideReverse[srcArc[j]];				

		string srcArc2 = srcArc.substr(parameters->hashLength);

                if(ID == srcNodeID && srcDirection == ret && (srcMulti==0 ||  this->multiplicity[i] == srcMulti) 
		   && (srcEdge.length()==0 || srcEdge == srcArc2 ))
                {
                	index = i;
                        break;
                }
	}
	if(index==-1)	return -1;
	
	return ret+(index<<2);
}

void node::extendEdge(unsigned long long neighbourID, string edge, parameter *parameters)
{
    string localedge;
    string neighbourStr = kmerGraph::longLongToString(neighbourID, parameters);

    for(int i=0;i<8;i++)
    {
	if(arcs[i].length()!=0)
	{
	    edge = arcs[i];
	    string descriptor = kmerGraph::longLongToString(nodeID, parameters);
	    if(i>=4)	reverse(descriptor.begin(), descriptor.end());
	    descriptor += edge;
	    int pos = descriptor.length() - parameters->hashLength;
	    string nextDescriptor = descriptor.substr(pos);

	    if(nextDescriptor == neighbourStr)	
	    {
		arcs[i] += neighbourStr.substr(parameters->hashLength-1) + edge;
	    }			
	}		
    }
}

class distNodeGraph
{
    public:
	hash_map<unsigned long long, arc  > kmolecules;
	hash_map<unsigned long long, node > nodes;
	vector<string> RemovedEdges;
	
	parameter *parameters;
	MPIEnviroment *MPIcontrol;
    	
	MPI_Request recvFuncRequest;	//monitering the Irecv Request
	char        recvFuncFlag;	//Huang up a Irecv Routine at the first time
	unsigned long long recvFuncData[3]; //used to store any command message received by Irecv Routine

	double  totalWorkTime, totalTime;		

	
	distNodeGraph(parameter *parameters, MPIEnviroment *MPIcontrol)
	{
	    this->parameters = parameters;
	    this->MPIcontrol = MPIcontrol;
	    this->RemovedEdges.clear();
	    recvFuncFlag = 0;
	    totalWorkTime= 0;
	    totalTime    = 0;
	}

	~distNodeGraph()
	{
	    kmolecules.clear();
	    nodes.clear();
	    RemovedEdges.clear();
	}

    public:
	unsigned long long getProcsID(unsigned long long kmerID, int hashLength, MPIEnviroment *MPIcontrol);

	void arcFrequency(parameter *parameters);
	void printNode(unsigned long long nodeID, FILE *fp);
	void printDistNodeGraph(parameter *parameters, MPIEnviroment *MPIcontrol);
	void printLocalNodeGraph(parameter *parameters, MPIEnviroment *MPIcontrol);
	void readDistNodeGraph(parameter *parameters, MPIEnviroment *MPIcontrol);
	
	void recvBufDistNodeGraph(char *recvbuf, parameter *parameters, MPIEnviroment *MPIcontrol);
	char* gatherDistNodeGraph(parameter *parameters, MPIEnviroment *MPIcontrol);
	bool checkDistNodeGraph();
	bool checkLocalNodeGraph(parameter *parameters);
	void constructDistNodeGraph(kmerGraph *kGraph, parameter *parameters, MPIEnviroment *MPIcontrol);
	void simplifyDistNodeGraph(parameter *parameters, MPIEnviroment *MPIcontrol);
	void cutoffGraph(MPIEnviroment *MPIcontrol, int threshold);
	void buildKmoleculeGraph(MPIEnviroment *MPIcontrol, parameter *parameters);
	void printContigs(parameter *parameters, MPIEnviroment *MPIcontrol);
	int  recvProc(parameter *parameters, MPIEnviroment *MPIcontrol, int n, MPI_Request *data_reqs);
	void tipsRemoval(parameter *parameters, MPIEnviroment *MPIcontrol);
	void bubbleRemoval(parameter *parameters, MPIEnviroment *MPIcontrol);

	//associated functions 
	void graphStatistics(parameter *parameters, MPIEnviroment *MPIcontrol);
	
	void masterSimplifyNodeGraph(parameter *parameters, MPIEnviroment *MPIcontrol); 
	void masterTipsRemoval(parameter  *parameters, MPIEnviroment *MPIcontrol);
	void masterLoopBubbleRemoval(parameter *parameters, MPIEnviroment *MPIcontrol);
	void masterMultipleEdgeBubbleRemoval(parameter *parameters, MPIEnviroment *MPIcontrol);
	
	void masterGraphStatistic(parameter *parameters, MPIEnviroment *MPIcontrol);
	void masterLowCoverageEdgeRemoval(parameter *parameters, MPIEnviroment *MPIcontrol);
	
	void masterRemoveCrossNode(parameter *parameters, MPIEnviroment *MPIcontrol);
	void masterRemoveCrossEdge(parameter *parameters, MPIEnviroment *MPIcontrol);
	void masterPrintContigs(parameter *parameters, MPIEnviroment *MPIcontrol);

	void stringContigs(parameter *parameters, MPIEnviroment *MPIcontrol);
	void printJungGraph(parameter *parameters, MPIEnviroment *MPIcontrol);
        void printKmerGraph(parameter *parameters, MPIEnviroment *MPIcontrol);
};


void distNodeGraph::arcFrequency(parameter *parameters)
{
   	hash_map<unsigned long long, arc>::iterator it;

	FILE *FP = fopen(parameters->arcFrequency, "w");
	
	int arcFreq[1000];
	for(int i=0;i<1000;i++)	arcFreq[i]=0;
	int tot = kmolecules.size();
   	for(it=kmolecules.begin(); it!=kmolecules.end(); it++)
   	{
		for(int i=0;i<8;i++)
		{
			if(it->second.multiplicity[i]<250)
				arcFreq[it->second.multiplicity[i]] ++;	
   			else	arcFreq[250]++;
		}
	}

	for(int i=0;i<1000;i++)	fprintf(FP, "%d\n", arcFreq[i]);
   	fclose(FP);	
}

/*
unsigned long long distNodeGraph::getProcsID(unsigned long long kmerID, int hashLength, MPIEnviroment *MPIcontrol)
{
    	unsigned long long ret = 0, tmpID = kmerID;
    	for(int i=0;i<hashLength;i++)
    	{
		ret = ret * 4;
//		ret |= ( (unsigned long long) 3 - (tmpID & (unsigned long long)3) );
		ret += ( 3 - tmpID%4 );
		tmpID = tmpID / 4;
    	}
//    	ret = ret ^ kmerID;
  	//  return ( (ret % (unsigned long long ) MersPrime) % (unsigned long long) MPIcontrol->nprocs);
  	//  return ( (ret ) % (unsigned long long) MPIcontrol->nprocs);
	if(ret<kmerID)	ret = kmerID;

	double tmp = ((sqrt(5.0)-1)/2) * ret;
	double rr  = tmp - floor(tmp);
	return floor(rr * MPIcontrol->nprocs);
}
*/

unsigned long long distNodeGraph::getProcsID(unsigned long long kmerID, int hashLength, MPIEnviroment *MPIcontrol)
{
    	unsigned long long revKmer = 0, tmpID = kmerID;
    	for(int i=0;i<hashLength;i++)
    	{
		revKmer = revKmer * 4;
		revKmer += ( 3 - tmpID%4 );
		tmpID = tmpID / 4;
    	}

	if(revKmer>kmerID) revKmer = kmerID;  

         unsigned int factor = 19;
         unsigned int numBytes = (hashLength + 3) / 4;
			
         unsigned int sum = 0;
         for(int i = 0; i < numBytes; i++){
                 sum = sum * factor + (revKmer & 0xFF);
                 revKmer >>= 8;
         }
         return sum % MPIcontrol->nprocs; 
}

/*
int Alsha::calcDestinationProc(AlshaKmerData kmerData)
{

         unsigned int factor = 19;
         unsigned int numBytes = (AlshaParams::KMERLENGTH + 3) / 4;
         unsigned int sum = 0;
         for(int i = 0; i < numBytes; i++){
                 sum = sum * factor + (kmerData & 0xFF);
                 kmerData >>= 8;
         }
         return sum % AlshaParams::numProcs;
}
*/
void distNodeGraph::printNode(unsigned long long nodeID, FILE *fp)
{
	if(nodes.find(nodeID) == nodes.end())
	{
		fprintf(fp, "|proc:%d nodeID = %llu can not find in this process\n", MPIcontrol->rank, nodeID);
		exit(0);
	}

	node tmp = nodes[nodeID];
        string descriptor = kmerGraph::longLongToString(tmp.nodeID, parameters);

	char line[20000000];

	int len=0;
	for(int i=0;i<8;i++)	len += nodes[nodeID].arcs[i].length();
	
	if(len+1000>10000000)	{printf("error \n"); exit(0);}
	
        sprintf(line, "|proc:%d >%llu@ %s# %s# %s# %s# %s# %s# %s# %s# %d# %d# %d# %d# %d# %d# %d# %d#", MPIcontrol->rank, nodes[nodeID].nodeID, nodes[nodeID].arcs[0].c_str(), nodes[nodeID].arcs[1].c_str(), nodes[nodeID].arcs[2].c_str(), nodes[nodeID].arcs[3].c_str(), nodes[nodeID].arcs[4].c_str(), nodes[nodeID].arcs[5].c_str(), nodes[nodeID].arcs[6].c_str(), nodes[nodeID].arcs[7].c_str(), nodes[nodeID].multiplicity[0], nodes[nodeID].multiplicity[1], nodes[nodeID].multiplicity[2], nodes[nodeID].multiplicity[3], nodes[nodeID].multiplicity[4], nodes[nodeID].multiplicity[5], nodes[nodeID].multiplicity[6], nodes[nodeID].multiplicity[7]);

/*
	if( !(line[1]>='0' &&  line[1]<='9') )	
	{
		printf("|proc:%d nodeID = %llu %llu %s error find in this process\n", MPIcontrol->rank, nodes[nodeID].nodeID, nodeID, line);
		exit(0);
	}
*/
	fprintf(fp, "%s\n", line);
	fflush(stdout);
/*
	printf("|proc:%d ID=%llu, Descriptor=%s %d\n",MPIcontrol->rank, tmp.nodeID, descriptor.c_str(), tmp.lockFlag);
        for(int i=0;i<4;i++)
        {
            if(!tmp.arcs[i].empty())
            printf("|proc:%d -(%s|%d)-",MPIcontrol->rank, tmp.arcs[i].c_str(), tmp.multiplicity[i]);
        }
        printf("\n");
        for(int i=4;i<8;i++)
        {
	        if(!tmp.arcs[i].empty())
	        printf("|proc:%d -(%s|%d)-",MPIcontrol->rank, tmp.arcs[i].c_str(), tmp.multiplicity[i]);
	}
	printf("\n|proc:%d--------------------------------------\n", MPIcontrol->rank);
*/
}

void distNodeGraph::printLocalNodeGraph(parameter *parameters, MPIEnviroment *MPIcontrol)
{
    hash_map<unsigned long long, node>::iterator it;

    FILE *fp = fopen(parameters->kmerGraph, "w");
    for(it=nodes.begin();it!=nodes.end();it++)
    {
	if(it->second.deleteFlag==0)	continue;

	//clear isolated nodes
	int arcFlag=0;
	for(int i=0;i<8;i++)	if(it->second.arcs[i].length() > 0)	arcFlag=1;
	if(arcFlag==0)		continue;	

        fprintf(fp, "%llu %s# %s# %s# %s# %s# %s# %s# %s# %d %d %d %d %d %d %d %d\n",it->second.nodeID, it->second.arcs[0].c_str(), it->second.arcs[1].c_str(), it->second.arcs[2].c_str(), it->second.arcs[3].c_str(), it->second.arcs[4].c_str(), it->second.arcs[5].c_str(), it->second.arcs[6].c_str(), it->second.arcs[7].c_str(), it->second.multiplicity[0], it->second.multiplicity[1], it->second.multiplicity[2], it->second.multiplicity[3], it->second.multiplicity[4], it->second.multiplicity[5], it->second.multiplicity[6], it->second.multiplicity[7]);
    }
    fclose(fp);
}



char* distNodeGraph::gatherDistNodeGraph(parameter *parameters, MPIEnviroment *MPIcontrol)
{
    hash_map<unsigned long long, node>::iterator it;
    vector<unsigned long long> vlen;

    unsigned long long int len = 0;	
    for(it=nodes.begin();it!=nodes.end();it++)
    {
	unsigned long long  length = 0;
	if(it->second.deleteFlag==0)	continue;

	//clear isolated nodes
	int arcFlag=0;
	for(int i=0;i<8;i++)	if(it->second.arcs[i].length() > 0)	arcFlag=1;
	if(arcFlag==0)		continue;	

	for(int i=0;i<8;i++)	length += it->second.arcs[i].length()+1;  //each string is end with '#'

	length += 20 + 8 * 3 + 16 + 1;  //20 nodeID, 24 multiplicity, 16 space, 1 for '\n';
	
	vlen.push_back(length);
	len += length;
    }

    unsigned long long *recvLen=NULL, *recvPos=NULL, Totlen=0;
    MPI_Request *recv_reqs = NULL;
    MPI_Status  *status  = NULL;

    if(MPIcontrol->rank == 0)
    {
	recvLen   = new unsigned long long [MPIcontrol->nprocs];
        assert(recvLen != NULL);
	recvPos   = new unsigned long long [MPIcontrol->nprocs];
        assert(recvPos != NULL);
	recv_reqs = new MPI_Request [MPIcontrol->nprocs];
        assert(recv_reqs != NULL);
	status    = new MPI_Status  [MPIcontrol->nprocs];
    	assert(status != NULL);
    }

    MPI_Gather(&len, 1, MPI_LONG_LONG_INT, recvLen, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);

    if(MPIcontrol->rank == 0)
    {
        recvPos[0] = 0;
        for(int i=1;i<MPIcontrol->nprocs;i++)	recvPos[i] = recvPos[i-1] + recvLen[i-1];
        Totlen = 0;
        for(int i=0;i<MPIcontrol->nprocs;i++)	Totlen += recvLen[i];
       
    	printf("|proc:%d totlen=%llu in this process\n", MPIcontrol->rank, Totlen);
    }
    
    char *buf = new char [len+1];
    assert(buf != NULL);
    char *pbuf = buf;
    int  index = 0;
    for(it=nodes.begin();it!=nodes.end();it++)
    {
	if(it->second.deleteFlag==0)	continue;

	//clear isolated nodes
	int arcFlag=0;
	for(int i=0;i<8;i++)	if(it->second.arcs[i].length() > 0)	arcFlag=1;
	if(arcFlag==0)		continue;	

        sprintf(pbuf, "%llu %s# %s# %s# %s# %s# %s# %s# %s# %d %d %d %d %d %d %d %d",it->second.nodeID, it->second.arcs[0].c_str(), it->second.arcs[1].c_str(), it->second.arcs[2].c_str(), it->second.arcs[3].c_str(), it->second.arcs[4].c_str(), it->second.arcs[5].c_str(), it->second.arcs[6].c_str(), it->second.arcs[7].c_str(), it->second.multiplicity[0], it->second.multiplicity[1], it->second.multiplicity[2], it->second.multiplicity[3], it->second.multiplicity[4], it->second.multiplicity[5], it->second.multiplicity[6], it->second.multiplicity[7]);
    
	for(int i=strlen(pbuf); i<vlen[index]; i++)    pbuf[i] = ' ';
	pbuf[vlen[index]-1] = '\n';
        
        pbuf += vlen[index];
        index ++; 
    }

    char *recvBuf=NULL;
    if(MPIcontrol->rank ==0)	
    {
	recvBuf = new char [Totlen+1];
        assert(recvBuf != NULL);
    }

    if(MPIcontrol->rank ==0)
    {
    	for(int i=0;i<MPIcontrol->nprocs; i++)
	{
		if(recvLen[i]>2147483647)	
		{
    			printf("|proc: %d data size is too large, MPI_recv will fail here on data gathering\n", MPIcontrol->rank);
			exit(0);
		}
//    		printf("|proc:%d Recving %llu data at %llu from Proc:%d (tag %d) \n", MPIcontrol->rank, recvLen[i], recvPos[i], i, i);
        	MPI_Irecv(recvBuf+recvPos[i], recvLen[i], MPI_CHAR, i, i, MPI_COMM_WORLD, &recv_reqs[i]);
   	} 
   }

    if(len>2147483647)	
    {
   	printf("|proc: %d data size is too large, MPI_recv will fail here on data gathering\n", MPIcontrol->rank);
	exit(0);
    }
//    printf("|proc:%d Sending %llu data to Proc:0(tag=%d) start\n", MPIcontrol->rank, len, MPIcontrol->rank);
    MPI_Send(buf, len, MPI_CHAR, 0, MPIcontrol->rank, MPI_COMM_WORLD);

    if(MPIcontrol->rank ==0)  MPI_Waitall(MPIcontrol->nprocs, recv_reqs,status);
//  printf("|proc:%d Sending %llu data to Proc:0 finished n", MPIcontrol->rank, len);

    if(MPIcontrol->rank == 0)  recvBuf[Totlen] = 0;

    delete(buf);
    vlen.clear();

    if(MPIcontrol->rank == 0)
    {
        delete(recvLen);
        delete(recvPos);
        delete(recv_reqs);
        delete(status);
    }
    if(MPIcontrol->rank ==0)	return recvBuf;
    else			return NULL;
}

void distNodeGraph::printDistNodeGraph(parameter *parameters, MPIEnviroment *MPIcontrol)
{
    hash_map<unsigned long long, node>::iterator it;
    vector<unsigned long long> vlen;

    unsigned long long int len = 0;	
    for(it=nodes.begin();it!=nodes.end();it++)
    {
	unsigned long long  length = 0;
	if(it->second.deleteFlag==0)	continue;

	//clear isolated nodes
	int arcFlag=0;
	for(int i=0;i<8;i++)	if(it->second.arcs[i].length() > 0)	arcFlag=1;
	if(arcFlag==0)		continue;	

	for(int i=0;i<8;i++)	length += it->second.arcs[i].length()+1;  //each string is end with '#'

	length += 20 + 8 * 3 + 16 + 1;  //20 nodeID, 24 multiplicity, 16 space, 1 for '\n';
	
	vlen.push_back(length);
	len += length;
    }

    unsigned long long *recvLen=NULL, *recvPos=NULL, Totlen=0;
    MPI_Request *recv_reqs = NULL;
    MPI_Status  *status  = NULL;

    if(MPIcontrol->rank == 0)
    {
	recvLen   = new unsigned long long [MPIcontrol->nprocs];
        assert(recvLen != NULL);
	recvPos   = new unsigned long long [MPIcontrol->nprocs];
        assert(recvPos != NULL);
	recv_reqs = new MPI_Request [MPIcontrol->nprocs];
        assert(recv_reqs != NULL);
	status    = new MPI_Status  [MPIcontrol->nprocs];
    	assert(status != NULL);
    }

    MPI_Gather(&len, 1, MPI_LONG_LONG_INT, recvLen, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);

    if(MPIcontrol->rank == 0)
    {
        recvPos[0] = 0;
        for(int i=1;i<MPIcontrol->nprocs;i++)	recvPos[i] = recvPos[i-1] + recvLen[i-1];
        Totlen = 0;
        for(int i=0;i<MPIcontrol->nprocs;i++)	Totlen += recvLen[i];
       
    	printf("|proc:%d totlen=%llu in this process\n", MPIcontrol->rank, Totlen);
    }
    
    char *buf = new char [len+1];
    assert(buf != NULL);
    char *pbuf = buf;
    int  index = 0;
    for(it=nodes.begin();it!=nodes.end();it++)
    {
	if(it->second.deleteFlag==0)	continue;

	//clear isolated nodes
	int arcFlag=0;
	for(int i=0;i<8;i++)	if(it->second.arcs[i].length() > 0)	arcFlag=1;
	if(arcFlag==0)		continue;	

        sprintf(pbuf, "%llu %s# %s# %s# %s# %s# %s# %s# %s# %d %d %d %d %d %d %d %d",it->second.nodeID, it->second.arcs[0].c_str(), it->second.arcs[1].c_str(), it->second.arcs[2].c_str(), it->second.arcs[3].c_str(), it->second.arcs[4].c_str(), it->second.arcs[5].c_str(), it->second.arcs[6].c_str(), it->second.arcs[7].c_str(), it->second.multiplicity[0], it->second.multiplicity[1], it->second.multiplicity[2], it->second.multiplicity[3], it->second.multiplicity[4], it->second.multiplicity[5], it->second.multiplicity[6], it->second.multiplicity[7]);
    
	for(int i=strlen(pbuf); i<vlen[index]; i++)    pbuf[i] = ' ';
	pbuf[vlen[index]-1] = '\n';
        
        pbuf += vlen[index];
        index ++; 
    }

    char *recvBuf=NULL;
    if(MPIcontrol->rank ==0)	
    {
	recvBuf = new char [Totlen+1];
        assert(recvBuf != NULL);
    }

    if(MPIcontrol->rank ==0)
    {
    	for(int i=0;i<MPIcontrol->nprocs; i++)
	{
		if(recvLen[i]>2147483647)	
		{
    			printf("|proc: %d data size is too large, MPI_recv will fail here on data gathering\n", MPIcontrol->rank);
			exit(0);
		}
//    		printf("|proc:%d Recving %llu data at %llu from Proc:%d (tag %d) \n", MPIcontrol->rank, recvLen[i], recvPos[i], i, i);
        	MPI_Irecv(recvBuf+recvPos[i], recvLen[i], MPI_CHAR, i, i, MPI_COMM_WORLD, &recv_reqs[i]);
   	} 
   }

    if(len>2147483647)	
    {
   	printf("|proc: %d data size is too large, MPI_recv will fail here on data gathering\n", MPIcontrol->rank);
	exit(0);
    }
//  printf("|proc:%d Sending %llu data to Proc:0(tag=%d) start\n", MPIcontrol->rank, len, MPIcontrol->rank);
    MPI_Send(buf, len, MPI_CHAR, 0, MPIcontrol->rank, MPI_COMM_WORLD);

    if(MPIcontrol->rank ==0)  MPI_Waitall(MPIcontrol->nprocs, recv_reqs,status);
//  printf("|proc:%d Sending %llu data to Proc:0 finished n", MPIcontrol->rank, len);

    if(MPIcontrol->rank == 0)  recvBuf[Totlen] = 0;

    if(MPIcontrol->rank == 0 )
    {
       	FILE *fp = fopen(parameters->graphPath, "w");
       	if(fp==NULL)	
        {
    		printf("|proc: %d Failed to Open File %s\n", MPIcontrol->rank, parameters->graphPath);
		exit(0);
    	}

    	fprintf(fp, "%s\n", recvBuf);
    	fclose(fp);
/*
        for(int i=0; i<MPIcontrol->nprocs; i++)
	{
	    char fileName[100];
	    strcpy(fileName, parameters->graphPath);
	    sprintf(fileName+strlen(fileName), "GatherGraph_%d", i);
	    fp = fopen(fileName, "w");

	    for(int j=recvPos[i]; j<recvPos[i]+recvLen[i];j++)  
		fprintf(fp, "%c", recvBuf[j]);

            fclose(fp); 	
	}	
*/  }

    delete(buf);
    vlen.clear();

    if(MPIcontrol->rank == 0)
    {
        delete(recvLen);
        delete(recvPos);
        delete(recv_reqs);
        delete(status);
        delete(recvBuf);
    }
}

void distNodeGraph::recvBufDistNodeGraph(char *recvbuf, parameter *parameters, MPIEnviroment *MPIcontrol)
{
	if(MPIcontrol->rank != 0)	return;
	
	node tmp;
	tmp.deleteFlag = 1;
    	if(pthread_mutex_init( &tmp.lockFlag,	NULL) != 0)
    	{
		printf("pthread_mutex_init error\n");
		exit(0);
    	}

	unsigned long long tmpID;
	int tmpMul;
	
	char *pbuf = recvbuf;
	char *ebuf;

	int counter = 0;
	while(*pbuf!=0)
	{
		tmpID = 0;
		while(*pbuf>='0' && *pbuf<='9')	tmpID = tmpID*10+ *pbuf++ - '0';
		tmp.nodeID = tmpID;

		while(*pbuf == ' ')  pbuf++; 
		
		for(int i=0;i<8;i++)
		{
			ebuf = pbuf;
			while(*ebuf != '#') ebuf++;
			*ebuf = '\0';
			
			tmp.arcs[i] = string(pbuf);
			pbuf = ebuf+2;
		}
				
		for(int i=0;i<8;i++)
		{
			tmpMul = 0;
			while(*pbuf>='0' && *pbuf<='9')	tmpMul = tmpMul*10+*pbuf++ -'0';
			tmp.multiplicity[i] = tmpMul;	

			while(*pbuf == ' ') pbuf++;
		}	

		nodes[tmp.nodeID] = tmp;
		while(*pbuf == '\n')	pbuf++;	
	}	
	delete(recvbuf);
	printf("|proc:%d %d nodes in distGraph\n", MPIcontrol->rank, this->nodes.size());
}

void distNodeGraph::readDistNodeGraph(parameter *parameters, MPIEnviroment *MPIcontrol)
{
	if(MPIcontrol->rank != 0)	return;

	FILE *fp = fopen(parameters->graphPath, "r");	
    	if(fp==NULL)	
    	{
    		printf("|proc: %d Failed to Open File %s\n", MPIcontrol->rank, parameters->graphPath);
		exit(0);
	}
	
	unsigned long long tmpID;
	int tmpMul;
	char *str = new char [100000000];

	node tmp;
	tmp.deleteFlag = 1;
    	if(pthread_mutex_init( &tmp.lockFlag,	NULL) != 0)
    	{
		printf("pthread_mutex_init error\n");
		exit(0);
    	}

	while(fscanf(fp, "%llu", &tmpID)!=EOF)
	{
		tmp.nodeID = tmpID;
		for(int i=0;i<8;i++)
		{
			fscanf(fp, "%s#", str);
			if(str[strlen(str)-1] == '#')	str[strlen(str)-1] = 0;
			tmp.arcs[i] = string(str);
			
		}
		for(int i=0;i<8;i++)
		{
			fscanf(fp, "%d", &tmpMul);
			tmp.multiplicity[i] = tmpMul;
		}
		nodes[tmp.nodeID] = tmp;
	}

	delete(str);
	printf("|proc:%d %d nodes in distGraph\n", MPIcontrol->rank, this->nodes.size());
}

bool distNodeGraph::checkDistNodeGraph()
{
    hash_map<unsigned long long, node>::iterator it;	
    for(it=nodes.begin(); it!=nodes.end();it++)
    {
   	unsigned long long srcNodeID = it->second.nodeID;
	
	if(it->second.deleteFlag == 0) continue;

	for(int k=0;k<8;k++)
	{
		if(nodes[srcNodeID].arcs[k].length()==0)	continue;
		string strdes = kmerGraph::longLongToString(srcNodeID, parameters);
		
		if(k>3)	
		{
			reverse(strdes.begin(), strdes.end());
			for(int j=0;j<strdes.length();j++)
				strdes[j] = parameters->nucleotideReverse[strdes[j]];
		}

		strdes += nodes[srcNodeID].arcs[k];

		int strpos = strdes.length() - parameters->hashLength;
		string newDes = strdes.substr(strpos);

		string reverseStr = newDes;
		reverse(reverseStr.begin(), reverseStr.end());
		for(int j=0;j<reverseStr.length();j++)
			reverseStr[j] = parameters->nucleotideReverse[reverseStr[j]];
		if(newDes<reverseStr)	newDes = reverseStr;
		assert(newDes.length() == parameters->hashLength);
		unsigned long long dstNodeID = kmerGraph::stringToLongLong(newDes.c_str(), 0, newDes.length(), parameters);
			
		int i, index=-1;
    	for(i=0;i<8;i++)
    	{
    		if(nodes[dstNodeID].arcs[i].length()==0)	continue;
			string descriptor = kmerGraph::longLongToString(dstNodeID, parameters);
		
    		if(i>3)	//negative side k-mers
			{
				reverse(descriptor.begin(), descriptor.end());

    			for(int j=0;j<descriptor.length();j++)	
				descriptor[j] = parameters->nucleotideReverse[descriptor[j]];
			}

    			descriptor += nodes[dstNodeID].arcs[i];

			int pos = descriptor.length()-parameters->hashLength;

			string retDescriptor = descriptor.substr(pos);

			string reverseStr = retDescriptor;
			reverse(reverseStr.begin(), reverseStr.end());
			for(int j=0;j<reverseStr.length();j++)
	    		reverseStr[j] = parameters->nucleotideReverse[reverseStr[j]];

			if(retDescriptor < reverseStr)	retDescriptor = reverseStr;
        		
			assert(retDescriptor.length() == parameters->hashLength);
    	   		unsigned long long ID = kmerGraph::stringToLongLong(retDescriptor.c_str(),0,retDescriptor.length(),parameters);
			
			string edge1 = descriptor;
			string edge2 = strdes;		
			reverse(edge2.begin(), edge2.end());
			for(int j=0;j<edge2.length();j++)	edge2[j] = parameters->nucleotideReverse[edge2[j]];
			
			if(ID == srcNodeID && nodes[srcNodeID].multiplicity[k] == nodes[dstNodeID].multiplicity[i] && edge1 == edge2)  break;
		}
		if(i==8)	
		{
			printf("src=%lld dst=%lld Coherence error in distNodeGraph\n", srcNodeID, dstNodeID);
			printNode(srcNodeID, stdout);
			printNode(dstNodeID, stdout);
			assert(0);
			exit(0);
			return false;
		}
	}
    }
    return true;
}

bool distNodeGraph::checkLocalNodeGraph(parameter *parameters)
{
    hash_map<unsigned long long, node>::iterator it;	
    for(it=nodes.begin(); it!=nodes.end();it++)
    {
   	unsigned long long srcNodeID = it->second.nodeID;
	
	if(it->second.deleteFlag == 0) continue;

	for(int k=0;k<8;k++)
	{
	    unsigned long long ID;
	    int IDindex;
            int Ret = it->second.getReverseEdge(k, ID, IDindex, parameters);			
	    assert( nodes[it->first].multiplicity[k]  == nodes[ID].multiplicity[IDindex]  );
	    assert( nodes[it->first].arcs[k].length() == nodes[ID].arcs[IDindex].length() ); 	
	}
    }
    return true;
}

void distNodeGraph::constructDistNodeGraph(kmerGraph *kGraph, parameter *parameters, MPIEnviroment *MPIcontrol)
{

	clock_t t0, t1;
	while(kGraph->constructKmerGraph(parameters, MPIcontrol)!=0)
	{
		kGraph->distributeKmerGraph(parameters, MPIcontrol);	

		t0 = clock();
		for(int i=0;i<kGraph->size; i++)
		{
			if(kmolecules.find(kGraph->kmers[i]) == kmolecules.end())
			{
				arc newArc;
				for(int t=0;t<8;t++)
				{
					if((kGraph->arcs[i])&(1<<t)) 	
						newArc.multiplicity[t] = 1;
					else	newArc.multiplicity[t] = 0;				
						
				}
				kmolecules[kGraph->kmers[i]] = newArc;
			}
			else
			{
				for(int t=0;t<8;t++)
				{
					if((kGraph->arcs[i])&(1<<t) && kmolecules[kGraph->kmers[i]].multiplicity[t]<255) 	
						kmolecules[kGraph->kmers[i]].multiplicity[t]++;
				}	
			}
		}

		delete kGraph->kmers;
		delete kGraph->arcs;
		kGraph->kmers=NULL;
		kGraph->arcs=NULL;
		kGraph->size = 0;
		t1= clock();
		kGraph->storagetime += t1 - t0;	
	}
}

void distNodeGraph::cutoffGraph(MPIEnviroment *MPIcontrol, int threshold)
{
   	hash_map<unsigned long long, arc>::iterator it;

	hash_map<unsigned long long, arc> ret;
	ret.clear();
	
	int tot = kmolecules.size();
   	for(it=kmolecules.begin(); it!=kmolecules.end(); it++)
   	{
		int cnt=0;
		for(int i=0;i<8;i++)
		{
			if(it->second.multiplicity[i] <threshold)	it->second.multiplicity[i]=0;
			if(it->second.multiplicity[i] !=0)	cnt++;
		}
		if(cnt>0)	ret[it->first] = it->second;
   	}
	kmolecules.clear();

	kmolecules = ret;
	printf("proc%d: Deleted %d Nodes(tot: %d)\n", MPIcontrol->rank, tot - ret.size(), tot);
}


void distNodeGraph::buildKmoleculeGraph(MPIEnviroment *MPIcontrol, parameter *parameters)
{
   	hash_map<unsigned long long, arc>::iterator it;
	node newNode;
	int kcount=0;
	char kbuf[100];
	for(it=kmolecules.begin(); it!=kmolecules.end(); it++)
	{
		kcount ++;
/*		if(kcount%100000==0)
		{
			sprintf(kbuf, "[%d-%d] building %d NodeGraph",kcount, kmolecules.size(), nodes.size());
			MPIcontrol->print(kbuf);
		}	
*/		newNode.init(it->first, it->second, parameters);
		nodes[it->first] = newNode;				
	}	
	kmolecules.clear();
	
	printf("proc%d: Build %d NodeGraph\n", MPIcontrol->rank, nodes.size());
}

void distNodeGraph::printContigs(parameter *parameters, MPIEnviroment *MPIcontrol)
{
    FILE *fp = fopen(parameters->contigsPath, "a");
    hash_map<unsigned long long, node>::iterator it;

    int tipsNum = 0;
    int contigNum = 0;	
    for(it=nodes.begin(); it!=nodes.end();it++)
    {
		if(it->second.deleteFlag==0)	continue;
		
		int edgeCount = 0;
		int edgeLength = -1;
		for(int i=0;i<8;i++)
		{

			if(it->second.arcs[i].length()==0)	continue;
			edgeCount++;
			edgeLength = it->second.arcs[i].length();

	    		if(it->second.arcs[i].length()<Contig_Length - parameters->hashLength)	continue;
			
    			string descriptor = kmerGraph::longLongToString(it->second.nodeID, parameters);
    			if(i>3)	//negative side k-mers
			{
				reverse(descriptor.begin(), descriptor.end());

    				for(int j=0;j<descriptor.length();j++)	
					descriptor[j] = parameters->nucleotideReverse[descriptor[j]];
			}

    			descriptor += it->second.arcs[i];

	    		int pos = descriptor.length()-parameters->hashLength;

	        	string retDescriptor = descriptor.substr(pos);

			string reverseStr = retDescriptor;
			reverse(reverseStr.begin(), reverseStr.end());
			for(int j=0;j<reverseStr.length();j++)
		    	reverseStr[j] = parameters->nucleotideReverse[reverseStr[j]];

			int ret;
			if(retDescriptor < reverseStr)
			{
			    	ret = -1;
			    	retDescriptor = reverseStr;
			}
			else   ret = 1;
	        	assert(retDescriptor.length() == parameters->hashLength);

	        	unsigned long long ID = kmerGraph::stringToLongLong(retDescriptor.c_str(),0,retDescriptor.length(),parameters);
		
			if(it->second.nodeID < ID)	continue;
		
			contigNum++;	
			fprintf(fp, ">proc%d: %s\n", MPIcontrol->rank, descriptor.c_str());
    			fflush(fp);
		}
		if(edgeCount==1 && edgeLength < parameters->hashLength) tipsNum++;
    }
    printf(">proc%d: %d tipsNum %d conitgs is found here\n", MPIcontrol->rank, tipsNum, contigNum);
    fclose(fp);
}

int distNodeGraph::recvProc(parameter *parameters, MPIEnviroment *MPIcontrol, int n, MPI_Request *data_reqs)
{
	
	unsigned long long data[3];
    	unsigned long long sdata[3];
    	unsigned long long srcNodeID, dstNodeID, msgType;
    	int length, srcKmerside;
    	unsigned char curMultiplicity;
	
	MPI_Request probeRequest;
	MPI_Status recvStatus, otherStats;
    	MPI_Status *data_status = new MPI_Status [MPIcontrol->nprocs];
    	int flag=0, peepFlag=0;
	int localFlag;

	int ret = 0;


    	clock_t clockStart, clockStop, workClockStart, workClockStop;
	clockStart = clock();

	while(1)
	{
		//if n!=0,  recv or send action was finished, jump out from the while. otherwise, we need to process packets.
		if(n!=0)	
		{
			MPI_Testall(n, data_reqs, &flag, data_status);
			ret = 0;
			if(flag !=0 )	
			{
				break;
			}
		}
/*	
		//we peep the packets before we continue. 
		peepFlag = 0;
		MPI_Iprobe(MPI_ANY_SOURCE, ServiceTunnel, MPI_COMM_WORLD, &peepFlag, &otherStats);

		// no packets exist, and we Test our data_reqs again. 
		if(peepFlag ==0)	{
			continue;
		}
*/
		
		workClockStart = clock();	
		//Hang up a Recv routine for any income packets to the ServiceTunnel. 	
		if(recvFuncFlag == 0)
		{	
			MPI_Irecv(recvFuncData, 3, MPI_LONG_LONG_INT, MPI_ANY_SOURCE, ServiceTunnel, MPI_COMM_WORLD, &recvFuncRequest);
			recvFuncFlag = 1;
		}

		MPI_Test(&recvFuncRequest, &flag, &recvStatus);	
		if(flag == 0)	continue;
		
//		MPI_Recv(data, 3, MPI_LONG_LONG_INT, MPI_ANY_SOURCE, ServiceTunnel, MPI_COMM_WORLD, &otherStats);
		if(LDEBUG)	printf("proc:%d Recv Msg (%llu %llu %llu)\n", MPIcontrol->rank, recvFuncData[0], recvFuncData[1], recvFuncData[2]);      

		//processing the packets acording to the protocol. 	
		srcNodeID = recvFuncData[0];
		dstNodeID = recvFuncData[1];
		msgType   = recvFuncData[2]%10;
		srcKmerside = recvFuncData[2]%100/10;   	//0 is negative side kmer , 1 is positive side kmer 	
		curMultiplicity = (recvFuncData[2]%100000)/100;
		length = recvFuncData[2]/100000;
			
		sdata[0]  = recvFuncData[1];
		sdata[1]  = recvFuncData[0];

		if(msgType != 9 && nodes.find(dstNodeID)==nodes.end())  
		{
	   		printf("proc:%d Error in find the Node (%llu) in the serviceProcess\n", MPIcontrol->rank, dstNodeID);
	   		exit(0);
		}
			
		unsigned long long  SendProc = getProcsID(srcNodeID, parameters->hashLength, MPIcontrol);

		if(msgType == 9)       //Last Msg Was Received
		{
			if(LDEBUG)	printf("proc:%d stoped\n", MPIcontrol->rank);
			ret = 1;  break;	
		}
		else if(msgType == 0)   //Lock Msg Was Received
		{
	   		//this node is on a link, then we lock this node.
	   		if(pthread_mutex_trylock( &nodes[dstNodeID].lockFlag) == 0)
	   		{	
				//Send back a LockSuccess msg to srcNodeID's ComputingThread through ComputingTunnel
				sdata[2] = 1;
				MPI_Send(sdata, 3, MPI_LONG_LONG_INT, SendProc, ComputingTunnel, MPI_COMM_WORLD);
				if(LDEBUG) printf("proc:%d Node(%llu)_%d Sending LockSuccessMsg to -[%llu](%llu)-\n",MPIcontrol->rank, dstNodeID, nodes[dstNodeID].deleteFlag, SendProc, srcNodeID);
	   		}
	   		else 	//This node is not on a link, or has been locked by other process.
	   		{
				//Send back a LockFailed msg to srcNodeID's ComputingThread through ComputingTunnel
				sdata[2] = 2;
				MPI_Send(sdata, 3, MPI_LONG_LONG_INT, SendProc, ComputingTunnel, MPI_COMM_WORLD);
				if(LDEBUG) printf("proc:%d Node(%llu)_%d Sending LockFailedMsg to -[%llu](%llu)-\n",MPIcontrol->rank, dstNodeID,nodes[dstNodeID].deleteFlag, SendProc, srcNodeID);
			}
	   	}
		else if(msgType == 3)
		{
	   		if(pthread_mutex_trylock( &nodes[dstNodeID].lockFlag) == 0)
	   		{
				printf("Unlock msg Received, But The node is not locked\n");
				pthread_mutex_unlock( &nodes[dstNodeID].lockFlag );
				exit(0);
	   		}
	   		pthread_mutex_unlock( &nodes[dstNodeID].lockFlag);
		}
		else if(msgType == 4)
		{
			char *strBuf = NULL;
			string str;
			if(length>0)
			{
	   			strBuf = new char [length + 1];	
	   			assert(strBuf != NULL);
				if(length!=0)
				{
					MPI_Recv(strBuf,length, MPI_CHAR, SendProc, DataTunnel, MPI_COMM_WORLD, &otherStats);
	   			}
				strBuf[length] = 0;
	   			str = string(strBuf);
	   			if(LDEBUG) printf("proc:%d Node(%llu) Recv Edge(%s) from -[%lld](%llu)-\n",MPIcontrol->rank, dstNodeID, strBuf, SendProc, srcNodeID);
			}

		    	//union the Edge here
		    	int index = nodes[dstNodeID].getSrcNodeEdge(srcNodeID, srcKmerside, string(""), 0, 2, parameters);

/*			if(dstNodeID == 3906247590979 && srcNodeID == 3560456582410)	{
				printNode(srcNodeID);
				printNode(dstNodeID);
				printf("srcNodeID=%llu srcKmerside=%d index=%d \n",  srcNodeID, srcKmerside, index);
			}
*/	    		if(index != -1)	
	    		{
				if(length!=0)
				{
					int preMulti = nodes[dstNodeID].multiplicity[(index>>2)] * nodes[dstNodeID].arcs[(index>>2)].length() + curMultiplicity * str.length();
	    				nodes[dstNodeID].arcs[(index>>2)] += str;
					nodes[dstNodeID].multiplicity[(index>>2)] = (unsigned char) (preMulti/nodes[dstNodeID].arcs[(index>>2)].length());
					
	    				//if(LDEBUG) printf("proc:%d Node(%llu) update Edge(%d) to be %s\n",MPIcontrol->rank, dstNodeID, index, nodes[dstNodeID].arcs[index]);
				}
				else	
				{
					nodes[dstNodeID].multiplicity[(index>>2)] = 0;
					nodes[dstNodeID].arcs[(index>>2)].clear(); 	
				}
	    		}
	    		else
	    		{
				printf("%ld %ld Edge was send to the wrong node, please Check\n", srcNodeID, dstNodeID);
				printNode(srcNodeID, stdout);
				printNode(dstNodeID, stdout);
				printf("---------------------\n");
		//		checkDistNodeGraph();
				assert(0);
				exit(0);
	    		}		

	    		delete strBuf;
	   		//Unlock Nodes here
	    		if(pthread_mutex_trylock( &nodes[dstNodeID].lockFlag) == 0)
	    		{
				printf("Unlock msg Received, But The node is not locked\n");
				pthread_mutex_unlock( &nodes[dstNodeID].lockFlag );
				//	exit(0);
	    		}		
	    		else	pthread_mutex_unlock( &nodes[dstNodeID].lockFlag);
		}
	
		//after finish previous packet, hang up a new Recv routine. 	
		MPI_Irecv(recvFuncData, 3, MPI_LONG_LONG_INT, MPI_ANY_SOURCE, ServiceTunnel, MPI_COMM_WORLD, &recvFuncRequest);

		workClockStop  = clock();
		totalWorkTime = totalWorkTime + (workClockStop-workClockStart);	
	}

	delete data_status;

	clockStop = clock();
	totalTime = totalTime + (clockStop-clockStart);

	return ret;
}

void distNodeGraph::simplifyDistNodeGraph(parameter *parameters, MPIEnviroment *MPIcontrol)
{
    	unsigned long long leftNodeID, rightNodeID;
    	int leftRet, rightRet; 
    	unsigned long long leftProc, rightProc, myProc;
    	int leftLock, rightLock, commTag,lockCount=0, curLockCount=0;
    	int leftLockResult, rightLockResult;
    	unsigned long long msgType = 0;
    	unsigned long long data[5], tmpdata[5], rdata[3];
    	hash_map<unsigned long long, node>::iterator it;
	MPI_Request *data_reqs = new MPI_Request [4]; 

    	if(!DEBUG)	printf("proc:%d Sending Thread Start\n",MPIcontrol->rank);

    	int node_num = nodes.size();
    	int kk = 1;
	int CircleNum = 0;
    	for(it=nodes.begin(); it!=nodes.end();it++)
    	{
		if(kk%100000==0)
		    	printf("proc:%d: processed %d (tot: %d)\n", MPIcontrol->rank, kk, node_num);

		kk++;
		
		//avoid all deleted nodes
		if(it->second.deleteFlag==0)  continue;	
		
		if(LDEBUG)	printf("|proc:%d processing %llu(%s)\n", MPIcontrol->rank, it->second.nodeID, kmerGraph::longLongToString(it->second.nodeID, parameters).c_str());
		
		//try to lock this node 
		if(pthread_mutex_trylock(&it->second.lockFlag) != 0)
		{
		    	if(LDEBUG) printf("*proc:%d Node(%llu)(%s) has locked by other process\n", MPIcontrol->rank, it->second.nodeID, kmerGraph::longLongToString(it->second.nodeID, parameters).c_str());
		    	lockCount++; 
	    		continue;	
		}

		leftRet  = it->second.getLeftNodeID(leftNodeID, parameters);
		rightRet = it->second.getRightNodeID(rightNodeID,parameters);

		if(LDEBUG)	printf("|proc:%d %llu leftNode=%llu(%d) rightNode=%llu(%d)\n", MPIcontrol->rank, it->second.nodeID, leftNodeID, leftRet, rightNodeID, rightRet);
		//	this node is not a semi-extended k-molecule 
		if( abs(leftRet) != 1 || abs(rightRet) != 1)
		{
			pthread_mutex_unlock(&it->second.lockFlag);
			continue;
		}

		leftProc    = getProcsID(leftNodeID, parameters->hashLength, MPIcontrol);
		rightProc   = getProcsID(rightNodeID, parameters->hashLength, MPIcontrol);
		myProc      = getProcsID(it->second.nodeID, parameters->hashLength, MPIcontrol);
		//assert(myProc==MPIcontrol->ank);
		//Processing LeftNode First, Send Lockmsg to leftNodeID's serviceThread 
		data[0] = it->second.nodeID;
		data[1] = leftNodeID;
		data[2] = 1*10;	//positive side
		
		MPI_Isend(data, 3, MPI_LONG_LONG_INT, leftProc, ServiceTunnel, MPI_COMM_WORLD, &data_reqs[0]);
		//Recv LockReply from leftNodeID,
		MPI_Irecv(rdata, 3, MPI_LONG_LONG_INT, leftProc, ComputingTunnel, MPI_COMM_WORLD, &data_reqs[1]);
		recvProc(parameters, MPIcontrol, 2, data_reqs);

		leftLockResult = rdata[2];	
		if(LDEBUG)	printf("|proc:%d Node(%llu)(%s) Sending LockMsg to -[%llu](%llu)(%s)-\n",MPIcontrol->rank, it->second.nodeID,kmerGraph::longLongToString(it->second.nodeID, parameters).c_str(), leftProc, leftNodeID, kmerGraph::longLongToString(leftNodeID, parameters).c_str());
		if(LDEBUG)	printf("|proc:%d Node(%llu)(%s) Recv LockReply from -[%llu](%llu)(%s) Ret=%d-\n",MPIcontrol->rank, it->second.nodeID, kmerGraph::longLongToString(it->second.nodeID, parameters).c_str(), leftProc, leftNodeID, kmerGraph::longLongToString(leftNodeID, parameters).c_str(), leftLockResult);
		fflush(stdout);
		
		//Processing RightNode, Send Lockmsg to rightNodeID's serviceThread  
		if(rightNodeID != leftNodeID)	
		{
			data[1] = rightNodeID;
			data[2] = 0*10;	//negative side
			MPI_Isend(data, 3, MPI_LONG_LONG_INT, rightProc, ServiceTunnel, MPI_COMM_WORLD, &data_reqs[0]);

			//Recv LockReply from rightNodeID
			MPI_Irecv(rdata, 3, MPI_LONG_LONG_INT, rightProc, ComputingTunnel, MPI_COMM_WORLD, &data_reqs[1]);
		
			if(LDEBUG) printf("|proc:%d Node(%llu)(%s) Sending LockMsg to -[%llu](%llu)(%s)-\n",MPIcontrol->rank, it->second.nodeID, kmerGraph::longLongToString(it->second.nodeID,parameters).c_str(), rightProc, rightNodeID, kmerGraph::longLongToString(rightNodeID, parameters).c_str());
			if(LDEBUG) printf("|proc:%d Node(%llu)(%s) Recv LockReply from -[%llu](%llu)(%s) Ret=%d-\n",MPIcontrol->rank, it->second.nodeID, kmerGraph::longLongToString(it->second.nodeID, parameters).c_str(), rightProc, rightNodeID, kmerGraph::longLongToString(rightNodeID, parameters).c_str(), rightLockResult);
			fflush(stdout);
			recvProc(parameters, MPIcontrol, 2, data_reqs);
			rightLockResult = rdata[2];
		}
		else	{ rightLockResult = leftLockResult;	CircleNum++; }

		if(leftLockResult==1 && rightLockResult==1)	
		{
	   	 	// Lock Success,  do some thing here.
		    	data[0] = it->second.nodeID;

		    	//Send leftNodeID with rightEdge
			unsigned char curMultiplicity, curMultiplicity2;
		    	char *pstr = const_cast<char *> (it->second.getLeftEdge(curMultiplicity));
	   	 	data[1] = rightNodeID;
	    		data[2] = 4 + 0*10 + 100*curMultiplicity + 100000*(unsigned long long)strlen(pstr);	//negative side kmer

			/*-------------------------------------------------------------------
			data[2] protocol formart: 
		  		Means for each dec position 
		  		0 is for protocol of SWAP.
		  		1 negative or positive side of kmer
		  		3,4,5: multiplicity of the edge 
          		6+:    length of the edge
		  	-------------------------------------------------------------------*/
		
	    		MPI_Isend(data, 3, MPI_LONG_LONG_INT, rightProc, ServiceTunnel, MPI_COMM_WORLD, &data_reqs[0]);
	    		MPI_Isend(pstr, data[2]/100000, MPI_CHAR, rightProc, DataTunnel, MPI_COMM_WORLD, &data_reqs[1]);
	    
	    	if(LDEBUG)       printf("|proc:%d Node(%llu)(%s) Send EdgeLength=%llu to -[%llu](%llu)(%s)-\n",MPIcontrol->rank, it->second.nodeID, kmerGraph::longLongToString(it->second.nodeID, parameters).c_str(), data[2]/100000, rightProc, rightNodeID, kmerGraph::longLongToString(rightNodeID, parameters).c_str());
	    	if(LDEBUG)       printf("|proc:%d Node(%llu)(%s) Send EdgeData(%s) to -[%llu](%llu)(%s)-\n",MPIcontrol->rank, it->second.nodeID, kmerGraph::longLongToString(it->second.nodeID, parameters).c_str(), pstr, rightProc, rightNodeID, kmerGraph::longLongToString(rightNodeID, parameters).c_str());
			fflush(stdout);
	//		recvProc(parameters, MPIcontrol, 2, data_reqs);
		
		//Send rightNodeID with leftEdge, add five s before data	
	    	char *pstr2 = const_cast<char *> (it->second.getRightEdge(curMultiplicity2));
		tmpdata[0] = it->second.nodeID;
	    	tmpdata[1] = leftNodeID;
	    	tmpdata[2] = 4 + 1*10 + 100*curMultiplicity2 + 100000*(unsigned long long)strlen(pstr2);	//positive side kmer
	    	MPI_Isend(tmpdata,  3, MPI_LONG_LONG_INT,     leftProc, ServiceTunnel, MPI_COMM_WORLD, &data_reqs[2]);
	    	MPI_Isend(pstr2, tmpdata[2]/100000, MPI_CHAR, leftProc, DataTunnel,    MPI_COMM_WORLD, &data_reqs[3]);	
		recvProc(parameters, MPIcontrol, 4, data_reqs);
	    /*	
	    	char *pstr2 = const_cast<char *> (it->second.getRightEdge(curMultiplicity));
		sdata[0] = it->second.nodeID;
	    	sdata[1] = leftNodeID;
	    	sdata[2] = 4 + 1*10 + 100*curMultiplicity + 100000*(unsigned long long)strlen(pstr);	//positive side kmer
	    	MPI_Isend(sdata, 3, MPI_LONG_LONG_INT, leftProc, ServiceTunnel, MPI_COMM_WORLD, &data_reqs[0]);
	    	MPI_Isend(pstr2, sdata[2]/100000, MPI_CHAR,    leftProc, DataTunnel, MPI_COMM_WORLD, &data_reqs[1]);	
		recvProc(parameters, MPIcontrol, 2, data_reqs);
	    */	
/*		if(LDEBUG)  printf("|proc:%d Node(%llu)(%s) Send EdgeLength=%llu to -[%llu](%llu)(%s)-\n",MPIcontrol->rank, it->second.nodeID, kmerGraph::longLongToString(it->second.nodeID, parameters).c_str(), data[2]/100000, leftProc, leftNodeID, kmerGraph::longLongToString(leftNodeID, parameters).c_str());
	    	if(LDEBUG) printf("|proc:%d Node(%llu)(%s) Send EdgeData(%s) to -[%llu](%llu)(%s)-\n",MPIcontrol->rank, it->second.nodeID, kmerGraph::longLongToString(it->second.nodeID, parameters).c_str(), pstr, leftProc, leftNodeID, kmerGraph::longLongToString(leftNodeID, parameters).c_str());
*/	    
		    	//Set DeleteFlag for this node
			it->second.deleteFlag = 0; 
		}
		else if(leftLockResult==1 || rightLockResult==1)
		{
	    		data[0] = it->second.nodeID;
	    		data[2] = 3;
	    		// Lock Failed,   do some thing here.
	    		if(leftLockResult == 1)
	    		{	
			//Send unlock msg to leftNodeID's serviceThread
			data[1] = leftNodeID;
			MPI_Isend(data, 3, MPI_LONG_LONG_INT, leftProc, ServiceTunnel, MPI_COMM_WORLD, &data_reqs[0]);
			recvProc(parameters, MPIcontrol, 1, data_reqs);
	
			if(LDEBUG)       printf("|proc:%d Node(%llu)(%s) Send UnLock to -[%llu](%llu)(%s)-\n",MPIcontrol->rank, it->second.nodeID, kmerGraph::longLongToString(it->second.nodeID,parameters).c_str(), leftProc, leftNodeID, kmerGraph::longLongToString(leftNodeID, parameters).c_str());
			}	
	    		else if(rightLockResult == 1)
	    		{
			//Send unlock msg to rightNodeID's serviceThread
			data[1] = rightNodeID;
			MPI_Isend(data, 3, MPI_LONG_LONG_INT, rightProc, ServiceTunnel, MPI_COMM_WORLD, &data_reqs[0]);
			recvProc(parameters, MPIcontrol, 1, data_reqs);

			if(LDEBUG)       printf("|proc:%d Node(%llu)(%s) Send UnLock to -[%llu](%llu)(%s)-\n",MPIcontrol->rank, it->second.nodeID, kmerGraph::longLongToString(it->second.nodeID, parameters).c_str(), rightProc, rightNodeID, kmerGraph::longLongToString(rightNodeID, parameters).c_str());
			}
		}
		pthread_mutex_unlock(&it->second.lockFlag);
    	}
   
	printf("proc%d: -----------------finshed\n", MPIcontrol->rank);
	fflush(stdout);
	
	MPI_Ibarrier(MPI_COMM_WORLD, &data_reqs[0]);
	recvProc(parameters,MPIcontrol, 1, data_reqs);

	delete data_reqs;
	data_reqs=NULL;   

	printf("*proc:%d  (%d) Node has locked by other process\n", MPIcontrol->rank, lockCount);
	printf("|proc:%d CircleNum = %d\n", MPIcontrol->rank, CircleNum);
}


void distNodeGraph::graphStatistics(parameter *parameters, MPIEnviroment *MPIcontrol)
{
	long long int nodeNum=0;	
	long long int edgeNum=0;
	long long int edgeAvgLength=0;
	long long int tipsNum = 0;
	long long int bubbleNum = 0;

	long long int recvNodeNum, recvEdgeNum, recvEdgeAvgLength, recvTipsNum, recvBubbleNum;

    	hash_map<unsigned long long, node>::iterator it;	
    	for(it=nodes.begin(); it!=nodes.end();it++)
    	{
   		unsigned long long srcNodeID = it->second.nodeID;
		if(it->second.deleteFlag == 0) continue;
		nodeNum++;
		
		int localEdgeNum=0;	
		int oneEdgeLength=0;
		int maxIndex=0,  maxMultiplicity = -1;
		int minIndex=0,  minMultiplicity = 1000000;

		for(int k=0;k<8;k++)
		{
			if(nodes[srcNodeID].arcs[k].length()==0)	continue;
			edgeNum ++;
			edgeAvgLength += nodes[srcNodeID].arcs[k].length();
				
				
			localEdgeNum++;
			oneEdgeLength = nodes[srcNodeID].arcs[k].length();


			if(nodes[srcNodeID].multiplicity[k] > maxMultiplicity)	maxIndex=k, maxMultiplicity = nodes[srcNodeID].multiplicity[k];
			if(nodes[srcNodeID].multiplicity[k] < minMultiplicity)	minIndex=k, minMultiplicity = nodes[srcNodeID].multiplicity[k];
		}
		if(localEdgeNum==1 && oneEdgeLength < parameters->hashLength)	tipsNum++;	
		if(minMultiplicity < 0.2 * maxMultiplicity && nodes[srcNodeID].arcs[minIndex].length() < 2 * parameters->hashLength)	bubbleNum ++;		
		
	}

//	printf("rank %d nodeNum = %ld\n",MPIcontrol->rank,  nodeNum);
//	printf("rank %d edgeNum = %ld\n",MPIcontrol->rank,  edgeNum);

	MPI_Reduce(&nodeNum, &recvNodeNum, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD); 		
	MPI_Reduce(&edgeNum, &recvEdgeNum, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD); 		
	MPI_Reduce(&edgeAvgLength, &recvEdgeAvgLength, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD); 		
	MPI_Reduce(&tipsNum, &recvTipsNum, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&bubbleNum, &recvBubbleNum, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	
	if(MPIcontrol->rank == 0)
	{
		printf("recvNodeNum = %ld\n",recvNodeNum);
		printf("recvEdgeNum = %ld\n", recvEdgeNum);
		printf("recvEdgeAvgLength = %ld\n", recvEdgeAvgLength); 
		printf("recvTipsNum = %ld\n", recvTipsNum);
		printf("recvBubbleNum = %ld\n", recvBubbleNum);
	}
}

void distNodeGraph::tipsRemoval(parameter *parameters, MPIEnviroment *MPIcontrol)
{
    unsigned long long leftNodeID, rightNodeID;
    int leftRet, rightRet; 
    unsigned long long leftProc, rightProc, myProc;
    int leftLock, rightLock, commTag,lockCount=0, curLockCount=0;
    int leftLockResult, rightLockResult;
    unsigned long long msgType = 0;
    unsigned long long data[5], rdata[3];
    hash_map<unsigned long long, node>::iterator it;
    MPI_Status stats;

	MPI_Request send_reqs, recv_reqs;

	MPI_Request *data_reqs   = new MPI_Request [MPIcontrol->nprocs];
	MPI_Status  *data_status = new MPI_Status  [MPIcontrol->nprocs];
	char *recvFinishTag      = new char        [MPIcontrol->nprocs];
	
	assert(recvFinishTag!=NULL);
	assert(data_reqs!=NULL);
	assert(data_status!=NULL);

	int flag;
	
    if(!DEBUG)	printf("proc:%d Sending Thread Start\n",MPIcontrol->rank);

    int node_num = nodes.size();
    int kk = 0;
    int CircleNum = 0;
    int tipsNum = 0;
    for(it=nodes.begin(); it!=nodes.end();it++)
    {
//	if(kk%100000==0)
    //		printf("proc:%d: processed %d (tot: %d)\n", MPIcontrol->rank, kk, node_num);

	//avoid all deleted nodes
	if(it->second.deleteFlag==0)  continue;	
	kk++;
	
	if(LDEBUG)	printf("|proc:%d processing %llu(%s)\n", MPIcontrol->rank, it->second.nodeID, kmerGraph::longLongToString(it->second.nodeID, parameters).c_str());
		
	
	//try to lock this node 
	if(pthread_mutex_trylock(&it->second.lockFlag) != 0)
	{
		if(LDEBUG) printf("*proc:%d Node(%llu)(%s) has locked by other process\n", MPIcontrol->rank, it->second.nodeID, kmerGraph::longLongToString(it->second.nodeID, parameters).c_str());
	    	lockCount++; 
    		continue;	
	}
	
	unsigned long long srcNodeID = it->second.nodeID;
	int localEdgeNum=0;	
	int oneEdgeLength=0;
	for(int k=0;k<8;k++)
	{
		if(nodes[srcNodeID].arcs[k].length()==0)	continue;

		localEdgeNum++;
		oneEdgeLength = nodes[srcNodeID].arcs[k].length();
	}
	if(localEdgeNum!=1 || oneEdgeLength > parameters->hashLength-1)	 
	{	
		pthread_mutex_unlock(&it->second.lockFlag);
		continue;
	}
	tipsNum++;	

	leftRet  = it->second.getLeftNodeID(leftNodeID, parameters);
	rightRet = it->second.getRightNodeID(rightNodeID,parameters);

//	printf("%llu %d %d\n", it->second.nodeID, leftRet, rightRet);
	//	this node is not a semi-extended k-molecule
		 
//	unsigned char curMultiplicity=0;
//    	char *pstr = const_cast<char *> (it->second.getLeftEdge(curMultiplicity));
		
	if(abs(leftRet) + abs(rightRet) != 1)
	{
		pthread_mutex_unlock(&it->second.lockFlag);
		continue;
	}

	if(abs(leftRet)==1)
	{
		leftProc    = getProcsID(leftNodeID, parameters->hashLength, MPIcontrol);
		myProc      = getProcsID(it->second.nodeID, parameters->hashLength, MPIcontrol);

		//Processing LeftNode First, Send Lockmsg to leftNodeID's serviceThread 
		data[0] = it->second.nodeID;
		data[1] = leftNodeID;
		data[2] = 1*10;	//positive side
		
		MPI_Isend(data, 3, MPI_LONG_LONG_INT, leftProc, ServiceTunnel, MPI_COMM_WORLD, &data_reqs[0]);
		//Recv LockReply from leftNodeID,
		MPI_Irecv(rdata, 3, MPI_LONG_LONG_INT, leftProc, ComputingTunnel, MPI_COMM_WORLD, &data_reqs[1]);
		recvProc(parameters, MPIcontrol, 2, data_reqs);
	
		leftLockResult = rdata[2];	
		if(LDEBUG)	printf("|proc:%d Node(%llu)(%s) Sending LockMsg to -[%llu](%llu)(%s)-\n",MPIcontrol->rank, it->second.nodeID,kmerGraph::longLongToString(it->second.nodeID, parameters).c_str(), leftProc, leftNodeID, kmerGraph::longLongToString(leftNodeID, parameters).c_str());
		if(LDEBUG)	printf("|proc:%d Node(%llu)(%s) Recv LockReply from -[%llu](%llu)(%s) Ret=%d-\n",MPIcontrol->rank, it->second.nodeID, kmerGraph::longLongToString(it->second.nodeID, parameters).c_str(), leftProc, leftNodeID, kmerGraph::longLongToString(leftNodeID, parameters).c_str(), leftLockResult);
		fflush(stdout);
	}
	else
	{
		rightProc   = getProcsID(rightNodeID, parameters->hashLength, MPIcontrol);
		myProc      = getProcsID(it->second.nodeID, parameters->hashLength, MPIcontrol);
		
		//Processing RightNode, Send Lockmsg to rightNodeID's serviceThread  
		data[0] = it->second.nodeID;
		data[1] = rightNodeID;
		data[2] = 0*10;	//negative side
		MPI_Isend(data, 3, MPI_LONG_LONG_INT, rightProc, ServiceTunnel, MPI_COMM_WORLD, &data_reqs[0]);

		//Recv LockReply from rightNodeID
		MPI_Irecv(rdata, 3, MPI_LONG_LONG_INT, rightProc, ComputingTunnel, MPI_COMM_WORLD, &data_reqs[1]);
		recvProc(parameters, MPIcontrol, 2, data_reqs);

		rightLockResult = rdata[2];
		if(LDEBUG) printf("|proc:%d Node(%llu)(%s) Sending LockMsg to -[%llu](%llu)(%s)-\n",MPIcontrol->rank, it->second.nodeID, kmerGraph::longLongToString(it->second.nodeID,parameters).c_str(), rightProc, rightNodeID, kmerGraph::longLongToString(rightNodeID, parameters).c_str());
		if(LDEBUG) printf("|proc:%d Node(%llu)(%s) Recv LockReply from -[%llu](%llu)(%s) Ret=%d-\n",MPIcontrol->rank, it->second.nodeID, kmerGraph::longLongToString(it->second.nodeID, parameters).c_str(), rightProc, rightNodeID, kmerGraph::longLongToString(rightNodeID, parameters).c_str(), rightLockResult);
		fflush(stdout);
	}

	if(abs(leftRet)==1 && leftLockResult==1)
	{
   		// Lock Success,  do some thing here.
    		data[0] = it->second.nodeID;
    		data[1] = leftNodeID;
    		data[2] = 4 + 1*10 + 100*0 + 100000*0;	//negative side kmer
			
	    	MPI_Isend(data, 3, MPI_LONG_LONG_INT, leftProc, ServiceTunnel, MPI_COMM_WORLD, &data_reqs[0]);
		recvProc(parameters, MPIcontrol, 1, data_reqs);
	    
	    	if(LDEBUG)       printf("|proc:%d Node(%llu)(%s) Send EdgeLength=%llu to -[%llu](%llu)(%s)-\n",MPIcontrol->rank, it->second.nodeID, kmerGraph::longLongToString(it->second.nodeID, parameters).c_str(), data[2]/100000, rightProc, rightNodeID, kmerGraph::longLongToString(rightNodeID, parameters).c_str());
		fflush(stdout);
		
		it->second.deleteFlag = 0; 
	}
	else if(abs(rightRet)==1 && rightLockResult==1)
	{
   		// Lock Success,  do some thing here.
    		data[0] = it->second.nodeID;
    		data[1] = rightNodeID;
    		data[2] = 4 + 0*10 + 100*0 + 100000*0;	//positive side kmer
    		MPI_Isend(data, 3, MPI_LONG_LONG_INT, rightProc, ServiceTunnel, MPI_COMM_WORLD, &data_reqs[0]);
		recvProc(parameters, MPIcontrol, 1, data_reqs);

		if(LDEBUG)		printf("|proc:%d Node(%llu)(%s) Send EdgeLength=%llu to -[%llu](%llu)(%s)-\n",MPIcontrol->rank, it->second.nodeID, kmerGraph::longLongToString(it->second.nodeID, parameters).c_str(), data[2]/100000, leftProc, leftNodeID, kmerGraph::longLongToString(leftNodeID, parameters).c_str());
    		//Set DeleteFlag for this node
		it->second.deleteFlag = 0; 
	}

	pthread_mutex_unlock(&it->second.lockFlag);
    }
   
	printf("proc%d: -----------------finshed\n", MPIcontrol->rank);
	fflush(stdout);

    	char sendFinishTag = 1;
	MPI_Isend(&sendFinishTag, 1, MPI_CHAR, 0, EndTunnel, MPI_COMM_WORLD, &send_reqs);
	
    	if(MPIcontrol->rank == 0)  
	{
		for(int i=0; i<MPIcontrol->nprocs; i++)	recvFinishTag[i] = 0;

		for(int i=0;i<MPIcontrol->nprocs;i++)
		{
			MPI_Irecv(&recvFinishTag[i], 1, MPI_CHAR, i, EndTunnel, MPI_COMM_WORLD, &data_reqs[i]);
			printf("proc:%d recving FinishTag from %d\n", MPIcontrol->rank, i);
			fflush(stdout);
		}
	
		recvProc(parameters, MPIcontrol, MPIcontrol->nprocs, data_reqs);	

		int sum = 0;
		for(int i=0; i<MPIcontrol->nprocs; i++)
	    	sum += recvFinishTag[i];

		printf("---------------------------%d\n", sum);
		fflush(stdout);
		if(sum != MPIcontrol->nprocs)
		{
	    	printf("Gather Tags Error\n");
	    	exit(0);
		}
	}

	recvProc(parameters, MPIcontrol, 1, &send_reqs);

	//	printf("proc:%d sending FinishTag \n", MPIcontrol->rank);
	//	fflush(stdout);

//  MPI_Gather(&sendFinishTag, 1, MPI_CHAR, recvFinishTag, 1, MPI_CHAR, root, MPI_COMM_WORLD);
   
   	unsigned long long endData[5000][3];	
	if(MPIcontrol->rank==0)
   	 {
		for(int i=0;i<MPIcontrol->nprocs;i++)
		{
			endData[i][0] = 0;
	    	endData[i][1] = i;
			endData[i][2] = 9;
	    	MPI_Isend(endData[i], 3, MPI_LONG_LONG_INT, i, ServiceTunnel, MPI_COMM_WORLD, &data_reqs[i]);
//			printf("proc:%d Sending (%llu %llu %llu)\n", MPIcontrol->rank, endData[i][0], endData[i][1], endData[i][2]);
//			printf("proc:%d sending killTag to %d\n", MPIcontrol->rank, i);
//			fflush(stdout);
		}
    	}

	while(recvProc(parameters, MPIcontrol, 0, data_reqs) == 0);	

	if(MPIcontrol->rank==0)
	{	
		MPI_Waitall(MPIcontrol->nprocs, data_reqs, data_status);
	}

	delete data_reqs;
	delete data_status;
	data_reqs = NULL;
	data_status = NULL;
	
	printf("*proc:%d  (%d) Node has locked by other process\n", MPIcontrol->rank, lockCount);
	printf("|proc:%d CircleNum = %d\n", MPIcontrol->rank, CircleNum);
}

void distNodeGraph::bubbleRemoval(parameter *parameters, MPIEnviroment *MPIcontrol)
{
    unsigned long long leftNodeID, rightNodeID;
    int leftRet, rightRet; 
    unsigned long long leftProc, rightProc, myProc;
    int leftLock, rightLock, commTag,lockCount=0, curLockCount=0;
    int leftLockResult, rightLockResult;
    unsigned long long msgType = 0;
    unsigned long long data[5], rdata[3];
    hash_map<unsigned long long, node>::iterator it;
    MPI_Status stats;

	MPI_Request send_reqs, recv_reqs;

	MPI_Request *data_reqs   = new MPI_Request [MPIcontrol->nprocs];
	MPI_Status  *data_status = new MPI_Status  [MPIcontrol->nprocs];
	char *recvFinishTag      = new char        [MPIcontrol->nprocs];
	
	assert(recvFinishTag!=NULL);
	assert(data_reqs!=NULL);
	assert(data_status!=NULL);

	int flag;
	
    if(!DEBUG)	printf("proc:%d Sending Thread Start\n",MPIcontrol->rank);

    int node_num = nodes.size();
    int kk = 0;
    int CircleNum = 0;
    int bubbleNum = 0;
    for(it=nodes.begin(); it!=nodes.end();it++)
    {
//	if(kk%100000==0)
//	printf("proc:%d: processed %d (tot: %d)\n", MPIcontrol->rank, kk, node_num);

	//avoid all deleted nodes
	if(it->second.deleteFlag==0)  continue;	
	kk++;
	
	if(LDEBUG)	printf("|proc:%d processing %llu(%s)\n", MPIcontrol->rank, it->second.nodeID, kmerGraph::longLongToString(it->second.nodeID, parameters).c_str());
	
	//try to lock this node 
	if(pthread_mutex_trylock(&it->second.lockFlag) != 0)
	{
		if(LDEBUG) printf("*proc:%d Node(%llu)(%s) has locked by other process\n", MPIcontrol->rank, it->second.nodeID, kmerGraph::longLongToString(it->second.nodeID, parameters).c_str());
	    	lockCount++; 
    		continue;	
	}
	
	unsigned long long srcNodeID = it->second.nodeID;
	int maxIndex=0,  maxMultiplicity = -1;
	int minIndex=0,  minMultiplicity = 1000000;
	for(int k=0;k<8;k++)
	{
		if(nodes[srcNodeID].arcs[k].length()==0)	continue;

		if(nodes[srcNodeID].multiplicity[k] > maxMultiplicity)	maxIndex=k, maxMultiplicity = nodes[srcNodeID].multiplicity[k];
		if(nodes[srcNodeID].multiplicity[k] < minMultiplicity)	minIndex=k, minMultiplicity = nodes[srcNodeID].multiplicity[k];
	}
	if(minMultiplicity < 0.2 * maxMultiplicity && nodes[srcNodeID].arcs[minIndex].length() < 2 * parameters->hashLength)	
	{	// this node has a bubble
		bubbleNum ++;		

		if(minIndex<4)	
		{
			rightRet = it->second.getNodeID(minIndex, rightNodeID, parameters);
			leftRet = 0;	leftNodeID = 0;			
		}
		else
		{
			leftRet = it->second.getNodeID(minIndex, leftNodeID, parameters);
			rightRet = 0;  rightNodeID = 0;
		}			
	}
	else
	{	// this node has no bubble.
		pthread_mutex_unlock(&it->second.lockFlag);
		continue;
	}	

	//leftRet  = it->second.getLeftNodeID(leftNodeID, parameters);
	//rightRet = it->second.getRightNodeID(rightNodeID,parameters);

//	printf("%llu %d %d\n", it->second.nodeID, leftRet, rightRet);
	//	this node is not a semi-extended k-molecule
		 
//	unsigned char curMultiplicity=0;
//    	char *pstr = const_cast<char *> (it->second.getLeftEdge(curMultiplicity));
		
	if(abs(leftRet) + abs(rightRet) != 1)
	{
		pthread_mutex_unlock(&it->second.lockFlag);
		continue;
	}

	if(abs(leftRet)==1)
	{
		leftProc    = getProcsID(leftNodeID, parameters->hashLength, MPIcontrol);
		myProc      = getProcsID(it->second.nodeID, parameters->hashLength, MPIcontrol);

		//Processing LeftNode First, Send Lockmsg to leftNodeID's serviceThread 
		data[0] = it->second.nodeID;
		data[1] = leftNodeID;
		data[2] = 1*10;	//positive side
		
		MPI_Isend(data, 3, MPI_LONG_LONG_INT, leftProc, ServiceTunnel, MPI_COMM_WORLD, &data_reqs[0]);
		//Recv LockReply from leftNodeID,
		MPI_Irecv(rdata, 3, MPI_LONG_LONG_INT, leftProc, ComputingTunnel, MPI_COMM_WORLD, &data_reqs[1]);
		recvProc(parameters, MPIcontrol, 2, data_reqs);
	
		leftLockResult = rdata[2];	
		if(LDEBUG)	printf("|proc:%d Node(%llu)(%s) Sending LockMsg to -[%llu](%llu)(%s)-\n",MPIcontrol->rank, it->second.nodeID,kmerGraph::longLongToString(it->second.nodeID, parameters).c_str(), leftProc, leftNodeID, kmerGraph::longLongToString(leftNodeID, parameters).c_str());
		if(LDEBUG)	printf("|proc:%d Node(%llu)(%s) Recv LockReply from -[%llu](%llu)(%s) Ret=%d-\n",MPIcontrol->rank, it->second.nodeID, kmerGraph::longLongToString(it->second.nodeID, parameters).c_str(), leftProc, leftNodeID, kmerGraph::longLongToString(leftNodeID, parameters).c_str(), leftLockResult);
		fflush(stdout);
	}
	else
	{
		rightProc   = getProcsID(rightNodeID, parameters->hashLength, MPIcontrol);
		myProc      = getProcsID(it->second.nodeID, parameters->hashLength, MPIcontrol);
		
		//Processing RightNode, Send Lockmsg to rightNodeID's serviceThread  
		data[0] = it->second.nodeID;
		data[1] = rightNodeID;
		data[2] = 0*10;	//negative side
		MPI_Isend(data, 3, MPI_LONG_LONG_INT, rightProc, ServiceTunnel, MPI_COMM_WORLD, &data_reqs[0]);

		//Recv LockReply from rightNodeID
		MPI_Irecv(rdata, 3, MPI_LONG_LONG_INT, rightProc, ComputingTunnel, MPI_COMM_WORLD, &data_reqs[1]);
		recvProc(parameters, MPIcontrol, 2, data_reqs);

		rightLockResult = rdata[2];
		if(LDEBUG) printf("|proc:%d Node(%llu)(%s) Sending LockMsg to -[%llu](%llu)(%s)-\n",MPIcontrol->rank, it->second.nodeID, kmerGraph::longLongToString(it->second.nodeID,parameters).c_str(), rightProc, rightNodeID, kmerGraph::longLongToString(rightNodeID, parameters).c_str());
		if(LDEBUG) printf("|proc:%d Node(%llu)(%s) Recv LockReply from -[%llu](%llu)(%s) Ret=%d-\n",MPIcontrol->rank, it->second.nodeID, kmerGraph::longLongToString(it->second.nodeID, parameters).c_str(), rightProc, rightNodeID, kmerGraph::longLongToString(rightNodeID, parameters).c_str(), rightLockResult);
		fflush(stdout);
	}

	if(abs(leftRet)==1 && leftLockResult==1)
	{
   		// Lock Success,  do some thing here.
    		data[0] = it->second.nodeID;
    		data[1] = leftNodeID;
    		data[2] = 4 + 1*10 + 100*0 + 100000*0;	//negative side kmer
			
	    	MPI_Isend(data, 3, MPI_LONG_LONG_INT, leftProc, ServiceTunnel, MPI_COMM_WORLD, &data_reqs[0]);
		recvProc(parameters, MPIcontrol, 1, data_reqs);
	    
	    	if(LDEBUG)       printf("|proc:%d Node(%llu)(%s) Send EdgeLength=%llu to -[%llu](%llu)(%s)-\n",MPIcontrol->rank, it->second.nodeID, kmerGraph::longLongToString(it->second.nodeID, parameters).c_str(), data[2]/100000, rightProc, rightNodeID, kmerGraph::longLongToString(rightNodeID, parameters).c_str());
		fflush(stdout);
		
		it->second.deleteFlag = 0; 
	}
	else if(abs(rightRet)==1 && rightLockResult==1)
	{
   		// Lock Success,  do some thing here.
    		data[0] = it->second.nodeID;
    		data[1] = rightNodeID;
    		data[2] = 4 + 0*10 + 100*0 + 100000*0;	//positive side kmer
    		MPI_Isend(data, 3, MPI_LONG_LONG_INT, rightProc, ServiceTunnel, MPI_COMM_WORLD, &data_reqs[0]);
		recvProc(parameters, MPIcontrol, 1, data_reqs);

		if(LDEBUG)		printf("|proc:%d Node(%llu)(%s) Send EdgeLength=%llu to -[%llu](%llu)(%s)-\n",MPIcontrol->rank, it->second.nodeID, kmerGraph::longLongToString(it->second.nodeID, parameters).c_str(), data[2]/100000, leftProc, leftNodeID, kmerGraph::longLongToString(leftNodeID, parameters).c_str());
    		//Set DeleteFlag for this node
		it->second.deleteFlag = 0; 
	}

	pthread_mutex_unlock(&it->second.lockFlag);
    }
   
	printf("proc%d: -----------------finshed\n", MPIcontrol->rank);
	fflush(stdout);

    	char sendFinishTag = 1;
	MPI_Isend(&sendFinishTag, 1, MPI_CHAR, 0, EndTunnel, MPI_COMM_WORLD, &send_reqs);
	
    	if(MPIcontrol->rank == 0)  
	{
		for(int i=0; i<MPIcontrol->nprocs; i++)	recvFinishTag[i] = 0;

		for(int i=0;i<MPIcontrol->nprocs;i++)
		{
			MPI_Irecv(&recvFinishTag[i], 1, MPI_CHAR, i, EndTunnel, MPI_COMM_WORLD, &data_reqs[i]);
			printf("proc:%d recving FinishTag from %d\n", MPIcontrol->rank, i);
			fflush(stdout);
		}
	
		recvProc(parameters, MPIcontrol, MPIcontrol->nprocs, data_reqs);	

		int sum = 0;
		for(int i=0; i<MPIcontrol->nprocs; i++)
	    	sum += recvFinishTag[i];

		printf("---------------------------%d\n", sum);
		fflush(stdout);
		if(sum != MPIcontrol->nprocs)
		{
	    	printf("Gather Tags Error\n");
	    	exit(0);
		}
	}

	recvProc(parameters, MPIcontrol, 1, &send_reqs);
   
   	unsigned long long endData[5000][3];	
	if(MPIcontrol->rank==0)
   	 {
		for(int i=0;i<MPIcontrol->nprocs;i++)
		{
			endData[i][0] = 0;
	    	endData[i][1] = i;
			endData[i][2] = 9;
	    	MPI_Isend(endData[i], 3, MPI_LONG_LONG_INT, i, ServiceTunnel, MPI_COMM_WORLD, &data_reqs[i]);
//			printf("proc:%d Sending (%llu %llu %llu)\n", MPIcontrol->rank, endData[i][0], endData[i][1], endData[i][2]);
//			printf("proc:%d sending killTag to %d\n", MPIcontrol->rank, i);
//			fflush(stdout);
		}
    	}

	while(recvProc(parameters, MPIcontrol, 0, data_reqs) == 0);	

	if(MPIcontrol->rank==0)
	{	
		MPI_Waitall(MPIcontrol->nprocs, data_reqs, data_status);
	}

	delete data_reqs;
	delete data_status;
	data_reqs = NULL;
	data_status = NULL;
	
	printf("*proc:%d  (%d) Node has locked by other process\n", MPIcontrol->rank, lockCount);
	printf("|proc:%d CircleNum = %d\n", MPIcontrol->rank, CircleNum);
}



void distNodeGraph::masterSimplifyNodeGraph(parameter *parameters, MPIEnviroment *MPIcontrol)
{
    	hash_map<unsigned long long, node>::iterator it;

	int num = 0;
    	for(it=nodes.begin();it!=nodes.end();it++)
	{
		if(it->second.deleteFlag == 0)	continue;
		int left=-1, right=-1;
		int leftNum=0, rightNum=0;	
		for(int i=0;i<4;i++)
			if(it->second.arcs[i].length()>0)	right = i, rightNum++;  // srcNode ->  right  	

		for(int i=4;i<8;i++)
			if(it->second.arcs[i].length()>0)	left  = i, leftNum++;  // left    <-  reverse srcNode

		if(rightNum != 1 || leftNum != 1)	continue;
		if(left == -1 || right == -1)		continue;
		
		num ++;
	}	
	printf("|proc:%d Total number of semi-extended nodes is %d\n", MPIcontrol->rank, num);


    	for(it=nodes.begin();it!=nodes.end();it++)
    	{
		if(it->second.deleteFlag == 0)	continue;
		int left=-1, right=-1;
		int leftNum=0, rightNum=0;	
		for(int i=0;i<4;i++)
			if(it->second.arcs[i].length()>0)	right = i, rightNum++;  // srcNode ->  right  	

		for(int i=4;i<8;i++)
			if(it->second.arcs[i].length()>0)	left  = i, leftNum++;  // left    <-  reverse srcNode

		if(rightNum != 1 || leftNum != 1)	continue;
		if(left == -1 || right == -1)		continue;
		

		int leftDirection, rightDirection;
		unsigned long long leftNode, rightNode;
		rightDirection = it->second.getNodeID(right, rightNode, parameters);
		leftDirection  = it->second.getNodeID(left,  leftNode, parameters);					

		if( nodes.find(rightNode) == nodes.end() ) 
		{
//			printf("Can not find %llu\n", rightNode );
		}
		if( nodes.find(leftNode) == nodes.end() )
		{
//			printf("can not find %llu\n", leftNode  );
		}

//		printf("%llu %llu %llu\n", it->first, rightNode, leftNode);
		int rightRet = nodes[rightNode].getSrcNodeEdge(it->first, 0, string(""), 0, 2, parameters);
		int leftRet  = nodes[leftNode].getSrcNodeEdge(it->first,  1, string(""), 0, 2,  parameters);	
	
			
	
		if(rightRet==-1 || leftRet ==-1)		
		{	
//		printf("|proc:%d Union (%llu[%d] %llu[%d %d] %llu[%d])\n", MPIcontrol->rank, leftNode,leftRet, it->first, left, right, rightNode, rightRet);
		continue;
		}

		rightRet = rightRet>>2;
		leftRet  = leftRet >>2; 

		if(leftNode==it->first || rightNode==it->first)	continue;

//		printf("|proc:%d Union (%llu[%d] %llu[%d %d] %llu[%d])\n", MPIcontrol->rank, leftNode,leftRet, it->first, left, right, rightNode, rightRet);

		//nodes[rightNode].arcs[rightRet] union with it->second.arcs[left]	
		nodes[rightNode].multiplicity[rightRet] = (nodes[rightNode].multiplicity[rightRet] * nodes[rightNode].arcs[rightRet].length() + it->second.multiplicity[left] * it->second.arcs[left].length()) / (nodes[rightNode].arcs[rightRet].length() + it->second.arcs[left].length()); 
		
		nodes[rightNode].arcs[rightRet] = nodes[rightNode].arcs[rightRet] + it->second.arcs[left];

		nodes[leftNode].multiplicity[leftRet] = (nodes[leftNode].multiplicity[leftRet] * nodes[leftNode].arcs[leftRet].length() + it->second.multiplicity[right]*it->second.arcs[right].length()) / (nodes[leftNode].arcs[leftRet].length() + it->second.arcs[right].length());

		nodes[leftNode].arcs[leftRet] = nodes[leftNode].arcs[leftRet] + it->second.arcs[right];
		
		it->second.deleteFlag = 0;
	}	

    	num = 0;
	for(it=nodes.begin();it!=nodes.end();it++)
    	{
		if(it->second.deleteFlag == 0)		continue;
		int left=-1, right=-1;
		int leftNum=0, rightNum=0;	
		for(int i=0;i<4;i++)
			if(it->second.arcs[i].length()>0)	right = i, rightNum++;  // srcNode ->  right  	

		for(int i=4;i<8;i++)
			if(it->second.arcs[i].length()>0)	left  = i, leftNum++;  // left    <-  reverse srcNode

		if(rightNum != 1 || leftNum != 1)	continue;
		if(left == -1 || right == -1)		continue;
		
		num ++;
	}	
	printf("|proc:%d Total number of semi-extended loop nodes is %d\n", MPIcontrol->rank, num);

}

void distNodeGraph::masterTipsRemoval(parameter  *parameters, MPIEnviroment *MPIcontrol)
{

    	hash_map<unsigned long long, node>::iterator it;

	int num = 0;
    	for(it=nodes.begin();it!=nodes.end();it++)
	{
		if(it->second.deleteFlag == 0)	continue;

		int arcNum = 0, index=-1;
		for(int i=0;i<8;i++)
			if(it->second.arcs[i].length()>0 )	arcNum ++, index = i;	

		if(arcNum==1 && it->second.arcs[index].length()<parameters->hashLength )	num++;		
	}
	
	printf("|proc:%d Total number of tip nodes is %d\n", MPIcontrol->rank, num);

		
    	for(it=nodes.begin();it!=nodes.end();it++)
	{
		if(it->second.deleteFlag == 0)	continue;
		int arcNum = 0, index=-1;
		for(int i=0;i<8;i++)
			if(it->second.arcs[i].length()>0 )	arcNum ++, index = i;
		if(arcNum==1 && it->second.arcs[index].length() < parameters->hashLength) 
		{
			int nextDirection;
			unsigned long long nextNode;		
			nextDirection = it->second.getNodeID(index, nextNode, parameters);

			int Ret = nodes[nextNode].getSrcNodeEdge(it->first, index<4?0:1, string(""), 0, 2, parameters);
			if(Ret==-1)
			{
				it->second.multiplicity[index] = 0;
				it->second.arcs[index].clear();
				continue;
			}
			Ret = Ret>>2;	
		
			if(nextNode == it->first)	continue;

	//		printf("Clean tips %llu -> %llu\n", it->first, nextNode);	

			nodes[nextNode].multiplicity[Ret] = 0;
			nodes[nextNode].arcs[Ret].clear();

			it->second.multiplicity[index] = 0;
			it->second.arcs[index].clear();	

			it->second.deleteFlag = 0;
		}
	}


	num = 0;
    	for(it=nodes.begin();it!=nodes.end();it++)
	{
		if(it->second.deleteFlag == 0)	continue;

		int arcNum = 0, index=-1;
		for(int i=0;i<8;i++)
			if(it->second.arcs[i].length()>0 )	arcNum ++, index = i;	

		if(arcNum==1 && it->second.arcs[index].length()<parameters->hashLength )	num++;		
	}
	
	printf("|proc:%d Total number of tip nodes is %d\n", MPIcontrol->rank, num);


}

void distNodeGraph::masterLowCoverageEdgeRemoval(parameter *parameters, MPIEnviroment *MPIcontrol)
{
    	hash_map<unsigned long long, node>::iterator it;

	unsigned long long TotalCovarage = 0;
	unsigned long long TotalNucleotide = 0;
    	for(it=nodes.begin();it!=nodes.end();it++)
	{
		if(it->second.deleteFlag == 0)	continue;

		for(int i=0;i<8;i++)
		{
			if(it->second.arcs[i].length()==0)	continue;
			TotalNucleotide += it->second.arcs[i].length();
			TotalCovarage += it->second.arcs[i].length() * it->second.multiplicity[i];
		}
	}
	unsigned long long avgCoverage =  TotalCovarage / TotalNucleotide;

	
    	for(it=nodes.begin();it!=nodes.end();it++)
	{
		if(it->second.deleteFlag == 0)	continue;

		for(int i=0;i<8;i++)
		{
			if(it->second.arcs[i].length() < 2 * parameters->hashLength && it->second.multiplicity[i] < 0.25 * avgCoverage )
			{
				it->second.arcs[i].clear();
				it->second.multiplicity[i] = 0;
			}
		}
	}
}

void distNodeGraph::masterLoopBubbleRemoval(parameter *parameters, MPIEnviroment *MPIcontrol)
{
    	hash_map<unsigned long long, node>::iterator it;

    	for(it=nodes.begin();it!=nodes.end();it++)
	{
		if(it->second.deleteFlag == 0)	continue;

		unsigned long long nextNode = 0,  nextNode2 = 0;
		int degree = 0;
		int loopDegree = 0;
		vector<int> loopIndex, loopDirection;
		vector<int> outEdge;
		for(int i=0;i<8;i++)
		{
			if(it->second.arcs[i].length()==0)	continue; // srcNode ->  nextNode 	
			degree ++ ;
			
			int Ret = it->second.getNodeID(i, nextNode, parameters);
			if(nextNode == it->first)	loopDegree ++, loopIndex.push_back(i), loopDirection.push_back(Ret);
			else				outEdge.push_back(i);
		}

//		printf("%llu loopDegree=%d outEdgeNum=%d\n", it->first, loopIndex.size(), outEdge.size());

		//isolated Nodes
		if(outEdge.size()==0)	continue;

		//each kmer size has one loop arc
		if(loopDegree == 2 && loopDirection[0] != loopDirection[1])
		{
			if(degree == 4)	//semi-extend node		
			{
				int Ret = it->second.getNodeID(outEdge[0], nextNode, parameters);
				int nextRevRet = nodes[nextNode].getSrcNodeEdge(it->first, outEdge[0]<4?0:1, string(""), 0, 2, parameters);
				if(nextRevRet == -1)	continue;
				int nextDirection = nextRevRet %2; 
				nextRevRet = nextRevRet >> 2;

				int Ret2 = it->second.getNodeID(outEdge[1], nextNode2, parameters);
				int nextRevRet2 = nodes[nextNode2].getSrcNodeEdge(it->first, outEdge[1]<4?0:1, string(""), 0, 2, parameters);
				if(nextRevRet2 == -1)	continue;
				int nextDirection2 = nextRevRet2 %2;
				nextRevRet2 = nextRevRet2 >> 2;
				if(nextDirection == nextDirection2)	continue; 	//This is a dead sink node 	
			}

			int Ret = it->second.getNodeID(outEdge[0], nextNode, parameters);
			int nextRevRet = nodes[nextNode].getSrcNodeEdge(it->first, outEdge[0]<4?0:1,string(""), 0,  2, parameters);
			if(nextRevRet == -1)	continue;

			int nextD = nextRevRet %2; 
			nextRevRet = nextRevRet >> 2;
			if(degree == 3 || degree == 4)	//tip node 
			{
//				if(nextD==0)		//positive size extending 
			//	int k = it->second.multiplicity[loopIndex[nextD]] / it->second.multiplicity[outEdge[0]];
			//	if(k==0)	k=1;			
				int k=1;
	
				it->second.multiplicity[loopIndex[nextD]] = (it->second.multiplicity[loopIndex[nextD]]*k*it->second.arcs[loopIndex[nextD]].length() + it->second.multiplicity[outEdge[0]]*it->second.arcs[outEdge[0]].length() )/(it->second.arcs[loopIndex[nextD]].length()*k+it->second.arcs[outEdge[0]].length());		
				string localArc = it->second.arcs[loopIndex[nextD]];
				for(int j=0;j<k-1;j++)	it->second.arcs[loopIndex[nextD]] += localArc;
				it->second.arcs[loopIndex[nextD]] += it->second.arcs[outEdge[0]];

				it->second.arcs[outEdge[0]].clear();
				it->second.multiplicity[outEdge[0]] = 0;	


				nodes[nextNode].multiplicity[nextRevRet] = (nodes[nextNode].multiplicity[nextRevRet]*nodes[nextNode].arcs[nextRevRet].length()+it->second.multiplicity[loopIndex[1-nextD]]*k*it->second.arcs[loopIndex[1-nextD]].length())/(nodes[nextNode].arcs[nextRevRet].length() + it->second.arcs[loopIndex[1-nextD]].length()*k);
				for(int j=0; j<k; j++)	nodes[nextNode].arcs[nextRevRet] += it->second.arcs[loopIndex[1-nextD]];
		
				it->second.arcs[loopIndex[1-nextD]].clear();
				it->second.multiplicity[loopIndex[1-nextD]] = 0;

/*
				printf("Node %llu union with %llu \n", it->first, nextNode);
				for(int j=0;j<8;j++)	printf("%s_%d ", it->second.arcs[j].c_str(), it->second.multiplicity[j]);
				printf("\n");		
	
				printf("Node %llu\n", nextNode);
				for(int j=0;j<8;j++)	printf("%s_%d ", nodes[nextNode].arcs[j].c_str(), nodes[nextNode].multiplicity[j]);		
				printf("\n");
*/
			}
		}
		if(loopDegree == 2 && loopDirection[0] == loopDirection[1])
		{
			// The edge path need pair-end information
		}
		if(loopDegree == 1)
		{
			if(degree == 2)	//tip node
			{
				int Ret = it->second.getNodeID(outEdge[0], nextNode, parameters);
			
				int nextRevRet = nodes[nextNode].getSrcNodeEdge(it->first, outEdge[0]<4?0:1,string(""), 0, 2, parameters);
				if(nextRevRet == -1) continue;
				nextRevRet = nextRevRet >> 2;
				
				if(   loopIndex[0]<4 && loopDirection[0]==-1 && outEdge[0]>3  	//Case 1. Extend the negative edge
				   || loopIndex[0]>3 && loopDirection[0]==1  && outEdge[0]<4 )	//Case 2. extend the positive edge
				{
					nodes[nextNode].multiplicity[nextRevRet] = (nodes[nextNode].multiplicity[nextRevRet] * nodes[nextNode].arcs[nextRevRet].length() + it->second.multiplicity[loopIndex[0]] * it->second.arcs[loopIndex[0]].length()/2)/(nodes[nextNode].arcs[nextRevRet].length() + it->second.arcs[loopIndex[0]].length()); 	
					nodes[nextNode].arcs[nextRevRet] += it->second.arcs[loopIndex[0]];

					it->second.multiplicity[loopIndex[0]] = ( it->second.multiplicity[loopIndex[0]] * it->second.arcs[loopIndex[0]].length()/2 + it->second.multiplicity[outEdge[0]] * it->second.arcs[outEdge[0]].length() ) / (it->second.arcs[loopIndex[0]].length() + it->second.arcs[outEdge[0]].length()) ; 
					it->second.arcs[loopIndex[0]] += it->second.arcs[outEdge[0]];

					it->second.arcs[outEdge[0]].clear();
					it->second.multiplicity[outEdge[0]] = 0;
				}
			}
			
			else if(degree == 3)	//semi-extend node
			{
				int Ret = it->second.getNodeID(outEdge[0], nextNode, parameters);
				int nextRevRet = nodes[nextNode].getSrcNodeEdge(it->first, outEdge[0]<4?0:1,string(""), 0,  2, parameters);
				if(nextRevRet==-1)	continue;
				int nextDirection = nextRevRet %2;
				nextRevRet = nextRevRet >> 2;
				

				int Ret2 = it->second.getNodeID(outEdge[1], nextNode2, parameters);
				int nextRevRet2 = nodes[nextNode2].getSrcNodeEdge(it->first, outEdge[1]<4?0:1,string(""), 0,  2, parameters);
				if(nextRevRet2==-1)	continue;
				int nextDirection2 = nextRevRet2 %2; 
				nextRevRet2 = nextRevRet2 >> 2;

				if( (nodes[nextNode2].multiplicity[outEdge[1]]*0.5 > it->second.multiplicity[loopIndex[0]] || nodes[nextNode].multiplicity[outEdge[0]]*0.5 > it->second.multiplicity[loopIndex[0]]) && nextDirection2  != nextDirection && it->second.arcs[loopIndex[0]].length()<5)
				{
					//delete this bubble edge, and we can continue extend the two other edges later
					it->second.arcs[loopIndex[0]].clear();
					it->second.multiplicity[loopIndex[0]] = 0;							
				}
				if(nextDirection2 == nextDirection && nextDirection == (loopIndex[0]>3?0:1) )
				{
					nodes[nextNode].multiplicity[nextRevRet] = (nodes[nextNode].multiplicity[nextRevRet] * nodes[nextNode].arcs[nextRevRet].length() + it->second.multiplicity[loopIndex[0]] * it->second.arcs[loopIndex[0]].length())/(nodes[nextNode].arcs[nextRevRet].length() + it->second.arcs[loopIndex[0]].length()); 	
					nodes[nextNode].arcs[nextRevRet] += it->second.arcs[loopIndex[0]];

					it->second.multiplicity[loopIndex[0]] = ( it->second.multiplicity[loopIndex[0]] * it->second.arcs[loopIndex[0]].length() + it->second.multiplicity[outEdge[0]] * it->second.arcs[outEdge[0]].length() ) / (it->second.arcs[loopIndex[0]].length() + it->second.arcs[outEdge[0]].length()) ; 
					it->second.arcs[loopIndex[0]] += it->second.arcs[outEdge[0]];

					it->second.arcs[outEdge[0]].clear();
					it->second.multiplicity[outEdge[0]] = 0;
				}
			}
		}	

		nextNode = 0;	
		int nextNodeRet = -1,  nextRevRet = -1;
		string arcStr;
		int arcMulti;
		int index  = -1;
		for(int i=0;i<8;i++)
		{
			if(it->second.arcs[i].length()==0)	continue; // srcNode ->  nextNode 	
			int Ret = it->second.getNodeID(i, nextNode, parameters);
			if(it->first == 46912496118442)		printf("%llu[%d] -> %llu %d\n", it->first, i, nextNode, Ret);	
			if(nextNode == it->first)
			{
				index       = i; 
				nextNodeRet = Ret;
				arcStr      = it->second.arcs[i];
				arcMulti    = it->second.multiplicity[i];
				
				int tmpDirection = (Ret==-1?1:0);  //negative size in -> positive size out
				nextRevRet  = it->second.getSrcNodeEdge(it->first, (index<4?0:1), string(""), 0,  tmpDirection, parameters);
				if(nextRevRet==-1)	continue;
				nextRevRet  = nextRevRet>>2;
//				printf("Node %llu: Length of looparc is%d(%d) (%d_%d, %d) degree %d\n", it->first, it->second.arcs[index].length(), it->second.multiplicity[index], index,Ret, nextRevRet, degree);	
			}
		}

	}
}


void distNodeGraph::masterMultipleEdgeBubbleRemoval(parameter *parameters, MPIEnviroment *MPIcontrol)
{
    	hash_map<unsigned long long, node>::iterator it;

	int totNodeNum = 0;
	int multipleEdgeNodeNum=0, resolvableNodeNum=0, diffNodeNum=0, sameNodeNum=0;
	int diffSolvNum = 0, diffSolvNum2 = 0,  diffSolvRemNum = 0;
	int positiveNum =0, negativeNum = 0, otherNum=0;

	unsigned long long nextNodei = 0, nextNodej = 0;
	int Directioni, Directionj;
	int retIndexi, retIndexj;
	int retDireci, retDirecj;
	vector<int> multiSrcIndex, multiEndIndex;
	vector<int> multiSrcDirec, multiEndDirec;
	vector<unsigned long long> multiSrc, multiEnd;
	vector<int> outSrcEdgeIndex, outEndEdgeIndex;
	vector<int> outSrcArcMul, outEndArcMul;
	vector<int> outSrcArcLen, outEndArcLen;

    	for(it=nodes.begin();it!=nodes.end();it++)
	{
		if(it->second.deleteFlag == 0)	continue;
		totNodeNum ++;

		//multiplicity edge node number
		for(int i=0;i<8;i++)
		{
			if(it->second.arcs[i].length() == 0)	continue;	
	
			Directioni = it->second.getNodeID(i, nextNodei, parameters);	
			if(Directioni==-1)	Directioni = 0;
			else			Directioni = 1;	
	
			for(int j=i+1;j<8;j++)
			{
				if(it->second.arcs[j].length()==0)	continue;

				Directionj = it->second.getNodeID(j,nextNodej,parameters);	
				if(Directionj == -1)	Directionj = 0;		//negative side
				else			Directionj = 1;		//positive side

				if(nextNodei == nextNodej) 
				{
					multipleEdgeNodeNum ++;
					int kdegree = 0, degree = 0;;
					for(int k=0;k<8;k++)
						if(it->second.arcs[k].length()>0)	degree++;						

/*					for(int k=0;k<8;k++)
						if(nodes[nextNodei].arcs[k].length() > 0)	kdegree ++;
	
					if(kdegree !=3 || degree != 3)	 
					{	
						printf("degree = %d kdegree = %d\n", degree, kdegree);		
						continue;	//Local resolveable Node will be left
					}	
*/
					if(degree != 3)		continue;
					resolvableNodeNum ++;
					
					//get outEdgeNodei and outEdgeNodej
					unsigned long long outEdgeNodei=0, outEdgeNodej=0, tmpNode;
					int outEdgeIndexi=-1, outEdgeIndexj=-1;
					int outEdgeReti, outEdgeRetj;						

					for(int k=0;k<8;k++)
					{
						if(it->second.arcs[k].length()==0)	continue;
						int Ret = it->second.getNodeID(k, tmpNode, parameters);
						if(tmpNode != nextNodei)	outEdgeNodei=tmpNode, outEdgeReti=(Ret==-1?0:1), outEdgeIndexi=k;
					}

					for(int k=0;k<8;k++)
					{
						if(nodes[nextNodei].arcs[k].length()==0)	continue;
						int Ret = nodes[nextNodei].getNodeID(k, tmpNode, parameters);
						if(tmpNode != it->first)	outEdgeNodej=tmpNode, outEdgeRetj=Ret=-1?0:1, outEdgeIndexj=k;
					}

			//		printf("%llu nextNodei&j=%llu, outEdgeNodei=%llu(%d), outEdgeNodej=%llu(%d)\n", it->first, nextNodei, outEdgeNodei, outEdgeIndexi, outEdgeNodej, outEdgeIndexj);
					

					//This function has error, fuck !
					//Get Index of reverse Edge from nextNodei or nextNodej
					retIndexi = nodes[nextNodei].getSrcNodeEdge(it->first, (i<4?0:1), it->second.arcs[i], it->second.multiplicity[i], 1-Directioni, parameters);
					retIndexj = nodes[nextNodej].getSrcNodeEdge(it->first, (j<4?0:1), it->second.arcs[j], it->second.multiplicity[j], 1-Directionj, parameters);
					if(retIndexi==-1 || retIndexj == -1)	continue;

					retDireci = retIndexi%2,  retDirecj = retIndexj%2;
					retIndexi = retIndexi>>2, retIndexj = retIndexj>>2;	

					int RetIndex=nodes[outEdgeNodei].getSrcNodeEdge(it->first,outEdgeIndexi<4?0:1,string(""), 0, 2,parameters);
					if(RetIndex==-1)	continue;

					int RetDirec=RetIndex%2;
					RetIndex = RetIndex >>2; 

					if(Directioni != Directionj )
					{
						diffNodeNum ++;
						if( !((i<4&&j>3) || (i>3&&j<4)) )  
						{
							diffSolvRemNum ++; 	
							if( (i<4 && RetDirec==1) || (i>3 && RetDirec==0)) 
							{ 	
		/*						int selfDirect = outEdgeIndexj<4?1:0;
								if(Directioni == selfDirect) 
								{
									it->second.arcs[j].clear();
									it->second.multiplicity[j] = 0;
									nodes[nextNodej].arcs[retIndexj].clear();
									nodes[nextNodej].multiplicity[retIndexj]=0;
								}
								else
								{
									it->second.arcs[i].clear();
									it->second.multiplicity[i]=0;
									nodes[nextNodei].arcs[retIndexi].clear();
									nodes[nextNodei].multiplicity[retIndexi]=0;
								}
								continue;	
		*/					}
							continue; 
						}

						if( (outEdgeIndexi<4&&outEdgeIndexj<4) || (outEdgeIndexi>3&&outEdgeIndexj>3) )  
						{
							diffSolvNum ++;	
//							if(it->second.multiplicity[i] <= it->second.multiplicity[j] && outEdgeIndexi<4 )
							if(outEdgeIndexi<4)
							{
								nodes[outEdgeNodei].arcs[RetIndex] = nodes[outEdgeNodei].arcs[RetIndex] + it->second.arcs[j] + nodes[nextNodei].arcs[retIndexi];
								it->second.arcs[outEdgeIndexi] = it->second.arcs[i] + nodes[nextNodej].arcs[retIndexj] + it->second.arcs[outEdgeIndexi];
								it->second.arcs[i].clear();
								it->second.multiplicity[i]=0;
								nodes[nextNodei].arcs[retIndexi].clear();
								nodes[nextNodei].multiplicity[retIndexi]=0;

							}
							else 
							{
								it->second.arcs[outEdgeIndexi] = it->second.arcs[j]+nodes[nextNodei].arcs[retIndexi] + it->second.arcs[outEdgeIndexi];
								nodes[outEdgeNodei].arcs[RetIndex] = nodes[outEdgeNodei].arcs[RetIndex] + it->second.arcs[i] + nodes[nextNodej].arcs[retIndexj];		
								it->second.arcs[j].clear();
								it->second.multiplicity[j] = 0;
								nodes[nextNodej].arcs[retIndexj].clear();
								nodes[nextNodej].multiplicity[retIndexj] = 0;

							}
						}
						else
						{
							diffSolvNum2++;
							if( outEdgeIndexi<4) 
							{
								
								nodes[outEdgeNodei].arcs[RetIndex] = nodes[outEdgeNodei].arcs[RetIndex] + it->second.arcs[j] + nodes[nextNodei].arcs[retIndexi];
								it->second.arcs[outEdgeIndexi] = it->second.arcs[i] + nodes[nextNodej].arcs[retIndexj] + it->second.arcs[outEdgeIndexi];
								it->second.arcs[i].clear();
								it->second.multiplicity[i] = 0;
								nodes[nextNodei].arcs[retIndexi].clear();
								nodes[nextNodei].multiplicity[retIndexi] = 0;

							}
							else 
							{

								it->second.arcs[outEdgeIndexi] = it->second.arcs[j]+nodes[nextNodei].arcs[retIndexi] + it->second.arcs[outEdgeIndexi];
								nodes[outEdgeNodei].arcs[RetIndex] = nodes[outEdgeNodei].arcs[RetIndex] + it->second.arcs[i] + nodes[nextNodej].arcs[retIndexj];		
								it->second.arcs[j].clear();
								it->second.multiplicity[j] = 0;
								nodes[nextNodej].arcs[retIndexj].clear();
								nodes[nextNodej].multiplicity[retIndexj] = 0;
							}
						}

					}
					else
					{
						sameNodeNum ++;
						if( (i<4 && j<4 && outEdgeIndexi>3) || (i>3 && j>3 && outEdgeIndexi<4))
						{
							string str;
							if(it->second.multiplicity[i]> it->second.multiplicity[j]) 
							{
								if(j<4)	{
									str = kmerGraph::longLongToString(it->second.nodeID, parameters);
									str += it->second.arcs[j];
								}
								else	{
									str = kmerGraph::longLongToString(nextNodej, parameters);
									str += nodes[nextNodej].arcs[retIndexj];
								}
													
								it->second.arcs[j].clear();
								it->second.multiplicity[j] = 0;
							//	printf("%llu (%d %d) %llu (%d %d) clean j\n", it->first,i,j,nextNodej, retIndexi,retIndexj);
								nodes[nextNodej].arcs[retIndexj].clear();
								nodes[nextNodej].multiplicity[retIndexj]=0;
							}
							else
							{
								if(i<4)	{
									str = kmerGraph::longLongToString(it->second.nodeID, parameters);
									str += it->second.arcs[i];
								}
								else	{
									str = kmerGraph::longLongToString(nextNodei, parameters);
									str += nodes[nextNodei].arcs[retIndexi];
								}
						
								it->second.arcs[i].clear();
								it->second.multiplicity[i]=0;
							//	printf("%llu (%d %d) %llu (%d %d) clean i\n", it->first,i,j,nextNodej, retIndexi,retIndexj);
								nodes[nextNodei].arcs[retIndexi].clear();
								nodes[nextNodei].multiplicity[retIndexi]=0;
							}
							RemovedEdges.push_back(str);
						}	
						else {				

							multiSrcIndex.push_back(i);
							multiEndIndex.push_back(j);
							multiSrcDirec.push_back(Directioni);
							multiEndDirec.push_back(Directionj);
							multiSrc.push_back(it->first);
							multiEnd.push_back(nextNodei);
/*							
							outSrcEdgeIndex.push_back(outEdgeIndexi);
							outEndEdgeIndex.push_back(outEdgeIndexj);
							outSrcArcMul.push_back(it->second.multiplicity[outEdgeIndexi]);
							outEndArcMul.push_back(nodes[nextNodei].multiplicity[outEdgeIndexj]);
							outSrcArcLen.push_back( it->second.arcs[outEdgeIndexi].length() );
							outEndArcLen.push_back( nodes[nextNodei].arcs[outEdgeIndexj].length() );
*/						}
					} 
				}
			}

		}
	}

	printf("Total node numbler is %d, multipleEdgeNodeNum is %d, resolvableNodeNum is %d, diffDirectionNode is (%d %d)%d, sameDirectionNode is %d \n", totNodeNum, multipleEdgeNodeNum, resolvableNodeNum, diffSolvNum,diffSolvNum2,  diffNodeNum, sameNodeNum);
	printf("diffSolvRemNum = %d, vector size = %d\n", diffSolvRemNum, multiSrcIndex.size());
	for(int i=0;i<multiSrcIndex.size();i++)
	{
//		printf("Bubble (%llu[%d] [%d] -> %llu) ArcLength (%d %d) Multiplicity (%d %d) Direction (%d %d) OutIndex (%d %d) [%d_%d %d_%d]\n", multiSrc[i], multiSrcIndex[i], multiEndIndex[i], multiEnd[i], nodes[multiSrc[i]].arcs[multiSrcIndex[i]].length(), nodes[multiSrc[i]].arcs[multiEndIndex[i]].length(),  nodes[multiSrc[i]].multiplicity[multiSrcIndex[i]], nodes[multiSrc[i]].multiplicity[multiEndIndex[i]], multiSrcDirec[i], multiEndDirec[i], outSrcEdgeIndex[i], outEndEdgeIndex[i], outSrcArcLen[i],  outSrcArcMul[i], outEndArcLen[i], outEndArcMul[i]);
	}
}

void distNodeGraph::masterGraphStatistic(parameter *parameters, MPIEnviroment *MPIcontrol)
{
    	hash_map<unsigned long long, node>::iterator it;

	unsigned long long TotalCovarage = 0;
	unsigned long long TotalNucleotide = 0;

	int totNodeNum = 0;
	int tipNodeNum = 0;
	int loopNodeNum = 0;
	int multipleEdgeNodeNum = 0;
	int resolvableNodeNum = 0;
	int crossNode = 0;
    	
	for(it=nodes.begin();it!=nodes.end();it++)
	{
		if(it->second.deleteFlag == 0)	continue;
		totNodeNum ++;

		for(int i=0;i<8;i++)
		{
			if(it->second.arcs[i].length()==0)	continue;
			TotalNucleotide += it->second.arcs[i].length();
			TotalCovarage += it->second.arcs[i].length() * it->second.multiplicity[i];
		}

		//tips number 
		int arcNum = 0, index=-1;
		for(int i=0;i<8;i++)
			if(it->second.arcs[i].length()>0 )	arcNum ++, index = i;	
		if(arcNum==1 && it->second.arcs[index].length()<parameters->hashLength )	tipNodeNum++;		

		//loop node number 
		unsigned long long nextNode[8];
		for(int i=0;i<8;i++)	nextNode[i] = 0;

		int degree = 0;
		for(int i=0;i<8;i++)
		{
			if(it->second.arcs[i].length()==0)	continue; // srcNode ->  nextNode 	
			degree ++ ;
			it->second.getNodeID(i, nextNode[i], parameters);
			if(nextNode[i] == it->first)	loopNodeNum ++;			
		}

		//multiplicity edge node number
		for(int i=0;i<8;i++)
		{
			for(int j=i+1;j<8;j++)
			{
				if(nextNode[i] == nextNode[j] && nextNode[i] !=0 ) 
				{
					multipleEdgeNodeNum ++;
					int kdegree = 0;
					for(int k=0;k<8;k++)
						if(nodes[nextNode[i]].arcs[k].length() > 0)	kdegree ++; 
					if(kdegree == 3 && degree == 3) resolvableNodeNum ++;
				}
			}

		}

		//Cross node number
		int positiveIn=0, positiveOut=0;
		int negativeIn=0, negativeOut=0; 
		for(int i=0;i<8;i++)
		{
			if(it->second.arcs[i].length()==0)	continue;
			if(i<4)	positiveOut ++,  negativeIn++;
			else	negativeOut ++,  positiveIn++;
		}	
		if(positiveIn == positiveOut ) crossNode++;
	}

	unsigned long long avgCoverage =  TotalCovarage / TotalNucleotide;
	printf("AvgCoverage = %d\n", avgCoverage);
	
	printf("Total node numbler is %d, tipNodeNum is %d, loopNodeNum is %d, multipleEdgeNodeNum is %d, resolvableNodeNum is %d CrossNode %d\n", totNodeNum, tipNodeNum,loopNodeNum, multipleEdgeNodeNum, resolvableNodeNum, crossNode);
}

void distNodeGraph::masterRemoveCrossNode(parameter *parameters, MPIEnviroment *MPIcontrol)
{
    	hash_map<unsigned long long, node>::iterator it;
	int crossNode[9];
	for(int i=0;i<9;i++)	crossNode[i]=0;

	for(it=nodes.begin();it!=nodes.end();it++)
	{
		if(it->second.deleteFlag == 0)	continue;

		//Cross node number
		int positiveIn=0, positiveOut=0;
		int negativeIn=0, negativeOut=0; 
		for(int i=0;i<8;i++)
		{
			if(it->second.arcs[i].length()==0)	continue;
			if(i<4)	positiveOut ++,  negativeIn++;
			else	negativeOut ++,  positiveIn++;
		}	
		if(positiveIn>=2 && positiveOut>=2) 
		{
			crossNode[positiveIn+positiveOut]++;
		
			for(int i=0;i<8;i++)	
			{
				if(it->second.arcs[i].length()==0)	continue;
			
			//	printf("%llu[%d] (arc:%d mul:%d)\n",it->first, i, it->second.arcs[i].length(), it->second.multiplicity[i]); 	
			}	
		}
		
		if(positiveIn==2 && positiveOut==2)
		{
			int positiveIndex[2], negativeIndex[2];	
			int tk=0;
			for(int t=0;t<4;t++)
			{
				if(it->second.arcs[t].length()>0)	positiveIndex[tk++] = t;
			}		

			tk=0;
			for(int t=4;t<8;t++)
			{
				if(it->second.arcs[t].length()>0)	negativeIndex[tk++] = t;
			}
		
			if(it->second.multiplicity[positiveIndex[0]]>it->second.multiplicity[positiveIndex[1]])  
			{
				int tmp = positiveIndex[0];
				positiveIndex[0] = positiveIndex[1];
				positiveIndex[1] = tmp;
			}				
			
			
			if(it->second.multiplicity[negativeIndex[0]]>it->second.multiplicity[negativeIndex[1]])  
			{
				int tmp = negativeIndex[0];
				negativeIndex[0] = negativeIndex[1];
				negativeIndex[1] = tmp; 
			}
		
			unsigned long long ID;
			int IDindex;

			for(int i=0;i<2;i++)
			{
				it->second.getReverseEdge(negativeIndex[i], ID, IDindex, parameters);
				nodes[ID].arcs[IDindex] += it->second.arcs[positiveIndex[i]];
				nodes[ID].multiplicity[IDindex] = (nodes[ID].multiplicity[IDindex] + it->second.multiplicity[positiveIndex[i]])/nodes[ID].arcs[IDindex].length();
				it->second.getReverseEdge(positiveIndex[i], ID, IDindex,parameters);
				nodes[ID].arcs[IDindex] += it->second.arcs[negativeIndex[i]];
				nodes[ID].multiplicity[IDindex] = (nodes[ID].multiplicity[IDindex] + it->second.multiplicity[negativeIndex[i]])/nodes[ID].arcs[IDindex].length();
			}
			it->second.deleteFlag = 0;
		}
	}
//	for(int i=0;i<9;i++)
//		printf("CrossNode[%d] = %d\t", i, crossNode[i]);
	printf("\n");
}
	
void distNodeGraph::masterRemoveCrossEdge(parameter *parameters, MPIEnviroment *MPIcontrol)
{
    	hash_map<unsigned long long, node>::iterator it;
	int crossEdge=0;
	
	
	for(it=nodes.begin();it!=nodes.end();it++)
	{
	/*	
	printf("Node %llu\n", 4373594623408022034);
	for(int i=0;i<8;i++)
	{
		if(nodes[4373594623408022034].arcs[i].length()==0)	continue;
		printf("[%d] (%d %d) \t", i, nodes[4373594623408022034].arcs[i].length(), nodes[4373594623408022034]. multiplicity[i]);
	}
*/
		if(it->second.deleteFlag == 0)	continue;
		
		int positiveOut=0, negativeOut=0;
		for(int i=0;i<8;i++)
		{
			if(it->second.arcs[i].length()==0)	continue;
			if(i<4)	positiveOut ++;
			else	negativeOut ++;
		}

		if( positiveOut+negativeOut != 3 )		continue;
		if(positiveOut==0 || negativeOut==0)		continue;		

		unsigned long long nextNodei;
		int indexi=-1, Retnext, Ret;
		unsigned long long ID, ID2;
		int IDindex, IDindex2;
		//Cross edge number
		for(int i=0;i<8;i++)
		{
			if(positiveOut==1 && i>3)	continue;
			if(positiveOut==2 && i<4)	continue;		
			
			if(it->second.arcs[i].length()==0)	continue;
/*
			it->second.getReverseEdge(i, nextNodei, Retnext, parameters);
			if(nodes[nextNodei].arcs[Retnext].length()==0)	continue;

			nodes[nextNodei].getReverseEdge(Retnext, ID, IDindex, parameters);
			if(ID==it->first && IDindex==i)			
				indexi = i;	
*/
			Ret = it->second.getNodeID(i,nextNodei,parameters);	
			Ret = (Ret==-1?0:1);	

			Retnext = nodes[nextNodei].getSrcNodeEdge(it->first, i<4?0:1, string(""), 0, 2, parameters);
			if(Retnext==-1)	continue;
			Retnext = Retnext>>2;	
			indexi = i;
		}	
		if(indexi==-1)	continue;

		if(it->first == nextNodei)	continue;
/*
		printf("Node %llu [%d] (%d %d)\t", it->first, indexi, it->second.arcs[indexi].length(), it->second.multiplicity[indexi] );
		for(int i=0;i<8;i++)
		{
			if(it->second.arcs[i].length()==0|| i==indexi)	continue;
			printf("[%d] (%d %d) \t", i, it->second.arcs[i].length(), it->second.multiplicity[i]);
		}
		printf("\nNode %llu [%d] (%d %d)\t", nextNodei, Retnext, nodes[nextNodei].arcs[Retnext].length(), nodes[nextNodei].multiplicity[Retnext]);
		for(int i=0;i<8;i++)
		{
			if(nodes[nextNodei].arcs[i].length()==0|| i==Retnext)	continue;
			printf("[%d] (%d %d) \t", i, nodes[nextNodei].arcs[i].length(), nodes[nextNodei].multiplicity[i]);
		}
		printf("\n\n");
*/		
		int degreeOut=0, degreetotal=0;
		for(int i=0;i<8;i++)
		{
			if(nodes[nextNodei].arcs[i].length()>0)	degreetotal++;
			
			if(Retnext<4 && i<4)	continue;
			if(Retnext>3 && i>3)	continue;	
			if(nodes[nextNodei].arcs[i].length()>0)	degreeOut++;
		}
		if(degreetotal!=3)	continue;		
		if(degreeOut!=2)	continue;

//		if(it->second.arcs[indexi].length() > 90-parameters->hashLength) continue;

	
		crossEdge ++;

		int leftIndex[2], rightIndex[2];		
		int tk=0;

		for(int i=0;i<8;i++)
		{
			if(positiveOut==1 && i<4)	continue;
			if(positiveOut==2 && i>3)	continue;
		
			if(it->second.arcs[i].length()==0)	continue;

			Ret = it->second.getNodeID(i, ID, parameters);		//get the ID	
			IDindex = nodes[ID].getSrcNodeEdge(it->first, i<4?0:1, string(""), 0, 2, parameters);
			if(IDindex==-1)	continue;
		//	it->second.getReverseEdge(i, ID, IDindex, parameters);
		//	nodes[ID].getReverseEdge(IDindex, ID2, IDindex2, parameters);
		//	if(ID2==it->first && IDindex2==i )	
			leftIndex[tk++] = i;			
		}	
		if(tk!=2)	continue;	
		assert(tk==2);
	
		tk = 0;
		for(int i=0;i<8;i++)
		{
			if(Retnext<4 && i<4)	continue;
			if(Retnext>3 && i>3)	continue;
			if(nodes[nextNodei].arcs[i].length()==0)	continue;


			Ret = nodes[nextNodei].getNodeID(i, ID, parameters);		//get the ID	
			IDindex = nodes[ID].getSrcNodeEdge(nextNodei, i<4?0:1, string(""), 0, 2, parameters);
			if(IDindex==-1)	continue;
		//	nodes[nextNodei].getReverseEdge(i, ID, IDindex,parameters);
		//	nodes[ID].getReverseEdge(IDindex, ID2, IDindex, parameters);
		//	if(ID2==nextNodei && IDindex2==i)	
			rightIndex[tk++] = i;	
		}
		if(tk!=2)	continue;
		assert(tk==2);

		if(it->second.multiplicity[leftIndex[0]] > it->second.multiplicity[leftIndex[1]] )
		{
			int tmp = leftIndex[0];
			leftIndex[0] = leftIndex[1];
			leftIndex[1] = tmp;
		}

		if(nodes[nextNodei].multiplicity[rightIndex[0]] > nodes[nextNodei].multiplicity[rightIndex[1]] )
		{
			int tmp = rightIndex[0];
			rightIndex[0] = rightIndex[1];
			rightIndex[1] = tmp;
		}	
	

//		it->second.getReverseEdge(indexi, ID, IDindex, parameters);
//		assert(nextNodei==ID && Retnext == IDindex);
//		nextNodei=ID, Retnext = IDindex;		

//		printf("Merging Cross Edge %llu %llu\n", it->first, nextNodei);	
//		printf("%d %d %d %d\n",  leftIndex[0], leftIndex[1], rightIndex[0], rightIndex[1]);

		Ret = it->second.getNodeID(leftIndex[0], ID, parameters);		//get the ID	
		IDindex = nodes[ID].getSrcNodeEdge(it->first, leftIndex[0]<4?0:1, string(""), 0, 2, parameters);
		if(IDindex==-1) continue;
		assert(IDindex!=-1);
		IDindex = IDindex>>2;	
//		it->second.getReverseEdge(leftIndex[0], ID, IDindex, parameters);
//		printf("%llu %d   ", ID, IDindex);
/*
		if(ID==3692893566669934926 || ID==4373594623408022034)
		{	
			printf("Node %llu [%d] (%d %d)\t", ID, IDindex, nodes[ID].arcs[IDindex].length(), nodes[ID].multiplicity[IDindex] );
			for(int i=0;i<8;i++)
			{
				if(nodes[ID].arcs[i].length()==0|| i==IDindex)	continue;
				unsigned long long tmpID; 
				int tmpIDindex;
				nodes[ID].getReverseEdge(i,tmpID, tmpIDindex, parameters);
				printf("[%d-->%llu-%d] (%d %d) \t", i, tmpID, tmpIDindex, nodes[ID].arcs[i].length(), nodes[ID].multiplicity[i]);
			}
		}
*/	
		nodes[ID].arcs[IDindex] +=  (it->second.arcs[indexi] + nodes[nextNodei].arcs[rightIndex[0]]);
		nodes[ID].multiplicity[IDindex] = ( (int) nodes[ID].multiplicity[IDindex] + it->second.multiplicity[indexi] + nodes[nextNodei].multiplicity[rightIndex[0]])/nodes[ID].arcs[IDindex].length();

		
		Ret = nodes[nextNodei].getNodeID(rightIndex[0], ID, parameters);		//get the ID	
		IDindex = nodes[ID].getSrcNodeEdge(nextNodei, rightIndex[0]<4?0:1, string(""), 0, 2, parameters);
		if(IDindex==-1)	continue;
		assert(IDindex!=-1);
		IDindex = IDindex>>2;	
//		nodes[nextNodei].getReverseEdge(rightIndex[0], ID, IDindex, parameters);
//		printf("%llu %d   ", ID, IDindex);	
		
/*		if(ID==3692893566669934926 || ID==4373594623408022034)
		{	
			printf("Node %llu [%d] (%d %d)\t", ID, IDindex, nodes[ID].arcs[IDindex].length(), nodes[ID].multiplicity[IDindex] );
			for(int i=0;i<8;i++)
			{
				if(nodes[ID].arcs[i].length()==0|| i==IDindex)	continue;
				unsigned long long tmpID;
				int tmpIDindex;
				nodes[ID].getReverseEdge(i,tmpID, tmpIDindex, parameters);
				printf("[%d-->%llu-%d] (%d %d) \t", i, tmpID, tmpIDindex, nodes[ID].arcs[i].length(), nodes[ID].multiplicity[i]);
			}
		}
*/
		nodes[ID].arcs[IDindex] +=  (nodes[nextNodei].arcs[Retnext] + it->second.arcs[leftIndex[0]]);
		nodes[ID].multiplicity[IDindex] = ( (int) nodes[ID].multiplicity[IDindex] + nodes[nextNodei].multiplicity[Retnext] + it->second.multiplicity[leftIndex[0]])/nodes[ID].arcs[IDindex].length();	

		Ret = it->second.getNodeID(leftIndex[1], ID, parameters);			//get the ID	
		IDindex = nodes[ID].getSrcNodeEdge(it->first, leftIndex[1]<4?0:1, string(""), 0, 2, parameters);
		if(IDindex==-1)	continue;
		assert(IDindex!=-1);
		IDindex = IDindex>>2;	
//		it->second.getReverseEdge(leftIndex[1], ID, IDindex, parameters);
//		printf("%llu %d   ", ID, IDindex);	
		
/*		if(ID==3692893566669934926 || ID==4373594623408022034)
		{	
			printf("Node %llu [%d] (%d %d)\t", ID, IDindex, nodes[ID].arcs[IDindex].length(), nodes[ID].multiplicity[IDindex] );
			for(int i=0;i<8;i++)
			{
				if(nodes[ID].arcs[i].length()==0|| i==IDindex)	continue;
				unsigned long long tmpID;
				int  tmpIDindex;
				nodes[ID].getReverseEdge(i,tmpID, tmpIDindex, parameters);
				printf("[%d-->%llu-%d] (%d %d) \t", i, tmpID, tmpIDindex, nodes[ID].arcs[i].length(), nodes[ID].multiplicity[i]);
			}
		}
*/
		nodes[ID].arcs[IDindex] +=  (it->second.arcs[indexi] + nodes[nextNodei].arcs[rightIndex[1]]);
		nodes[ID].multiplicity[IDindex] = ( (int) nodes[ID].multiplicity[IDindex] + it->second.multiplicity[indexi] + nodes[nextNodei].multiplicity[rightIndex[1]])/nodes[ID].arcs[IDindex].length();

		Ret = nodes[nextNodei].getNodeID(rightIndex[1], ID, parameters);		//get the ID	
		IDindex = nodes[ID].getSrcNodeEdge(nextNodei, rightIndex[1]<4?0:1, string(""), 0, 2, parameters);
		if(IDindex==-1)	continue;
		assert(IDindex!=-1);
		IDindex = IDindex>>2;	
//		nodes[nextNodei].getReverseEdge(rightIndex[1], ID, IDindex, parameters);
//		printf("%llu %d   ", ID, IDindex);	
		
/*		if(ID==3692893566669934926 || ID==4373594623408022034)
		{	
			printf("Node %llu [%d] (%d %d)\t", ID, IDindex, nodes[ID].arcs[IDindex].length(), nodes[ID].multiplicity[IDindex] );
			for(int i=0;i<8;i++)
			{
				if(nodes[ID].arcs[i].length()==0|| i==IDindex)	continue;
				unsigned long long tmpID;
				int tmpIDindex;
				nodes[ID].getReverseEdge(i,tmpID, tmpIDindex, parameters);
				printf("[%d-->%llu-%d] (%d %d) \t", i, tmpID, tmpIDindex, nodes[ID].arcs[i].length(), nodes[ID].multiplicity[i]);
			}
		}
*/		nodes[ID].arcs[IDindex] +=  (nodes[nextNodei].arcs[Retnext]+it->second.arcs[leftIndex[1]]);
		nodes[ID].multiplicity[IDindex] = ((int) nodes[ID].multiplicity[IDindex] + nodes[nextNodei].multiplicity[Retnext] + it->second.multiplicity[leftIndex[1]])/nodes[ID].arcs[IDindex].length();	
			
		it->second.deleteFlag = 0;	
		nodes[nextNodei].deleteFlag = 0;

//		printf("end\n");
/*
		printf("Node %llu [%d] (%d %d)\t", it->first, indexi, it->second.arcs[indexi].length(), it->second.multiplicity[indexi] );
		for(int i=0;i<8;i++)
		{
			if(it->second.arcs[i].length()==0|| i==indexi)	continue;
			printf("[%d] (%d %d) \t", i, it->second.arcs[i].length(), it->second.multiplicity[i]);
		}
		printf("\nNode %llu [%d] (%d %d)\t", nextNodei, Retnext, nodes[nextNodei].arcs[Retnext].length(), nodes[nextNodei].multiplicity[Retnext]);
		for(int i=0;i<8;i++)
		{
			if(nodes[nextNodei].arcs[i].length()==0|| i==Retnext)	continue;
			printf("[%d] (%d %d) \t", i, nodes[nextNodei].arcs[i].length(), nodes[nextNodei].multiplicity[i]);
		}
		printf("\n\n");
*/
	}
	printf("crossEdge is %d\n", crossEdge);

}


void distNodeGraph::stringContigs(parameter *parameters, MPIEnviroment *MPIcontrol)
{
    hash_map<unsigned long long, node>::iterator it;

    stringstream contigs;
    int tipsNum = 0;
    int contigNum = 0;	
    for(it=nodes.begin(); it!=nodes.end();it++)
    {
		if(it->second.deleteFlag==0)	continue;
		
		int edgeCount = 0;
		int edgeLength = -1;
		for(int i=0;i<8;i++)
		{
			if(it->second.arcs[i].length()==0)	continue;
			edgeCount++;
			edgeLength = it->second.arcs[i].length();

	    		if(it->second.arcs[i].length()<Contig_Length - parameters->hashLength)	continue;
			
    			string descriptor = kmerGraph::longLongToString(it->second.nodeID, parameters);
    			if(i>3)	//negative side k-mers
			{
				reverse(descriptor.begin(), descriptor.end());

    				for(int j=0;j<descriptor.length();j++)	
					descriptor[j] = parameters->nucleotideReverse[descriptor[j]];
			}

    			descriptor += it->second.arcs[i];

	    		int pos = descriptor.length()-parameters->hashLength;

	        	string retDescriptor = descriptor.substr(pos);

			string reverseStr = retDescriptor;
			reverse(reverseStr.begin(), reverseStr.end());
			for(int j=0;j<reverseStr.length();j++)
			    	reverseStr[j] = parameters->nucleotideReverse[reverseStr[j]];

			int ret;
			if(retDescriptor < reverseStr)
			{
			    	ret = -1;
			    	retDescriptor = reverseStr;
			}
			else   ret = 1;
	        	assert(retDescriptor.length() == parameters->hashLength);

	        	unsigned long long ID = kmerGraph::stringToLongLong(retDescriptor.c_str(),0,retDescriptor.length(),parameters);
		
			if(it->second.nodeID < ID)	continue;

			contigNum++;	
			
			contigs<<">contig_NodeID"<<it->first<<"_Index"<<i<<"_Len"<<descriptor.size()<<"_Mul"<<(int) (it->second.multiplicity[i])<<"_"<<contigNum<<"\n"<<descriptor<<"\n";	
				
		}
		if(edgeCount==1 && edgeLength < parameters->hashLength) tipsNum++;
    }

    int removedContigNum = 0;
    for(int i=0;i<RemovedEdges.size();i++)
    {
	if(RemovedEdges[i].length()>Contig_Length) 
		contigs<<">contig_RemovedNo"<<removedContigNum<<"\n"<<RemovedEdges[i]<<"\n";
		removedContigNum++;	
    }

    MPIcontrol->File_write(parameters->contigsPath, contigs.str().c_str(), contigs.str().length());	
    printf(">proc%d: tipsNum %d  removedContigNum %d, Contigs %d is found here\n", MPIcontrol->rank, tipsNum, removedContigNum, contigNum);
}

void distNodeGraph::masterPrintContigs(parameter *parameters, MPIEnviroment *MPIcontrol)
{
    if(MPIcontrol->rank != 0 ) return;
	
    FILE *fp = fopen(parameters->masterContigPath, "w");
    if(fp==NULL)	
    {
	printf("Failed in Open masterContig.txt\n");
	exit(0); 
    }

    hash_map<unsigned long long, node>::iterator it;

    int tipsNum = 0;
    int contigNum = 0;	
    for(it=nodes.begin(); it!=nodes.end();it++)
    {
		if(it->second.deleteFlag==0)	continue;
		
		int edgeCount = 0;
		int edgeLength = -1;
		for(int i=0;i<8;i++)
		{
			if(it->second.arcs[i].length()==0)	continue;
			edgeCount++;
			edgeLength = it->second.arcs[i].length();

	    		if(it->second.arcs[i].length()<Contig_Length - parameters->hashLength)	continue;
			
    			string descriptor = kmerGraph::longLongToString(it->second.nodeID, parameters);
    			if(i>3)	//negative side k-mers
			{
				reverse(descriptor.begin(), descriptor.end());

    				for(int j=0;j<descriptor.length();j++)	
					descriptor[j] = parameters->nucleotideReverse[descriptor[j]];
			}

    			descriptor += it->second.arcs[i];

	    		int pos = descriptor.length()-parameters->hashLength;

	        	string retDescriptor = descriptor.substr(pos);

			string reverseStr = retDescriptor;
			reverse(reverseStr.begin(), reverseStr.end());
			for(int j=0;j<reverseStr.length();j++)
			    	reverseStr[j] = parameters->nucleotideReverse[reverseStr[j]];

			int ret;
			if(retDescriptor < reverseStr)
			{
			    	ret = -1;
			    	retDescriptor = reverseStr;
			}
			else   ret = 1;
	        	assert(retDescriptor.length() == parameters->hashLength);

	        	unsigned long long ID = kmerGraph::stringToLongLong(retDescriptor.c_str(),0,retDescriptor.length(),parameters);
		
			if(it->second.nodeID < ID)	continue;
		
			contigNum++;	
			
			fprintf(fp, ">contig_NodeID%llu_Index%d_Len%d_Mul%d_%d\n%s\n", it->first, i, descriptor.size(),it->second.multiplicity[i], contigNum, descriptor.c_str());
		}
		if(edgeCount==1 && edgeLength < parameters->hashLength) tipsNum++;
    }

    int removedContigNum = 0;
    for(int i=0;i<RemovedEdges.size();i++)
    {
	if(RemovedEdges[i].length()>Contig_Length) 
		fprintf(fp,">contig_RemovedNo%d\n%s\n", removedContigNum++, RemovedEdges[i].c_str());
    }

    printf(">proc%d: tipsNum %d  removedContigNum %d, Contigs %d is found here\n", MPIcontrol->rank, tipsNum, removedContigNum, contigNum);
    fclose(fp);
}

void distNodeGraph::printJungGraph(parameter *parameters, MPIEnviroment *MPIcontrol)
{
	FILE *fp1 = fopen(parameters->JungGraph_arc, "w");
	FILE *fp2 = fopen(parameters->JungGraph_mul, "w");
       	if(fp1==NULL || fp2==NULL)	
        {
    		printf("|proc: %d Failed to Open File %s\n", MPIcontrol->rank, parameters->graphPath);
		exit(0);
    	}
	
	if(MPIcontrol->rank != 0 ) return;

    	hash_map<unsigned long long, node>::iterator it;
    	for(it=nodes.begin();it!=nodes.end();it++)
    	{
		if(it->second.deleteFlag == 0)	continue;
			
		for(int i=0;i<8;i++)	
		{
			if(it->second.arcs[i].length()>0)
			{
				unsigned long long nextID;
				it->second.getNodeID(i, nextID, parameters);					
				fprintf(fp1, "%llu\t%llu\t%llu\n", it->first, nextID, it->second.arcs[i].length());//, it->second.multiplicity[i]);
				fprintf(fp2, "%llu\t%llu\t%llu\n", it->first, nextID, it->second.multiplicity[i]);
			}
		}
	}
	fclose(fp1);
	fclose(fp2);
}

void distNodeGraph::printKmerGraph(parameter *parameters, MPIEnviroment *MPIcontrol)
{
	FILE *fp = fopen(parameters->kmerGraph, "w");

       	if(fp==NULL)	
        {
    		printf("|proc: %d Failed to Open File %s\n", MPIcontrol->rank, parameters->graphPath);
		exit(0);
    	}
	
	if(MPIcontrol->rank != 0 ) return;

    	hash_map<unsigned long long, node>::iterator it;
    	for(it=nodes.begin();it!=nodes.end();it++)
    	{
		if(it->second.deleteFlag == 0)	continue;
		
		string descriptor = kmerGraph::longLongToString(it->first, parameters);

		string reverseStr = descriptor;
    		reverse(reverseStr.begin(), reverseStr.end());
    		for(int i=0;i<reverseStr.length();i++)
        		reverseStr[i] = parameters->nucleotideReverse[reverseStr[i]];

		unsigned long long ID = kmerGraph::stringToLongLong(reverseStr.c_str(), 0, reverseStr.length(), parameters);

		unsigned long long Next[4];
		int kc=0;
		for(int i=0;i<4;i++)
		{
			if(it->second.arcs[i].length()>0)
			{
				unsigned long long nextID;
				int Ret = it->second.getNodeID(i, nextID, parameters);
	
				if(Ret==-1)
				{				
					string nReverseStr = kmerGraph::longLongToString(nextID, parameters);

    					reverse(nReverseStr.begin(), nReverseStr.end());
    					for(int j=0;j<nReverseStr.length();j++)
        					nReverseStr[i] = parameters->nucleotideReverse[nReverseStr[i]];

					nextID = kmerGraph::stringToLongLong(nReverseStr.c_str(), 0, nReverseStr.length(), parameters);
				}
				Next[kc++] = nextID;	
			}		
		}

                if(kc>0)
                {
                        fprintf(fp, "%llu", it->first);
                        for(int i=0;i<kc;i++)
                                fprintf(fp, " %llu",  Next[i]);
                        fprintf(fp, "\n");
                }

                kc = 0;
                for(int i=4;i<8;i++)
                {
                        if(it->second.arcs[i].length()>0)
                        {
                                unsigned long long nextID;
                                int Ret = it->second.getNodeID(i, nextID, parameters);

                                if(Ret==-1)
                                {
                                        string nReverseStr = kmerGraph::longLongToString(nextID, parameters);

                                        reverse(nReverseStr.begin(), nReverseStr.end());
                                        for(int j=0;j<nReverseStr.length();j++)
                                                nReverseStr[i] = parameters->nucleotideReverse[nReverseStr[i]];

                                        nextID = kmerGraph::stringToLongLong(nReverseStr.c_str(), 0, nReverseStr.length(), parameters);
                                }
                                Next[kc++] = nextID;
                        }
                }

                if(kc>0)
                {
                        fprintf(fp, "%llu", ID);
                        for(int i=0;i<kc;i++)
                                fprintf(fp, " %llu",  Next[i]);
                        fprintf(fp, "\n");
                }
	}
	fclose(fp);
}


int main (int argc, char *argv[])
{
    clock_t t0, t1, t2, t3, t4, t5, t6, t7;
    t0 = clock();
    FILE *FP;

    char message[300];
    int i, size, namelen, lgsize;
    MPIEnviroment MPIcontrol;
    MPIcontrol.init(argc, argv);
    
    //Get Reads
    parameter parameters;
    parameters.getParameters(argc, argv, &MPIcontrol);

    if(MPIcontrol.rank==0)	FP = fopen(parameters.LogPath, "w");

//    sequence *sequences = new sequence;
//    assert(sequences!=NULL);
//    sequences->getSequences(&parameters, &MPIcontrol);

    MPI_Barrier(MPI_COMM_WORLD);
    if(MPIcontrol.rank==0)
    {	
	t1 = clock();
	sprintf(message, "Time spend in Geting Reads %.5f\n", (t1-t0)/(float)CLOCKS_PER_SEC);
	MPIcontrol.print(message, FP);
//	fprintf(FP, "Time spend in Geting Reads %.5f\n", (t1-t0)/(float)CLOCKS_PER_SEC);
    }
    
    //Distributed K-mer to processors
    kmerGraph *mygraph = new kmerGraph;

    //Construct De Bruijn Graph
    distNodeGraph *nodeGraph = new distNodeGraph(&parameters, &MPIcontrol);
    assert(nodeGraph!=NULL);
    nodeGraph->constructDistNodeGraph(mygraph, &parameters, &MPIcontrol);
  
    double communicationTime;
    MPI_Reduce(&mygraph->commtime, &communicationTime, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    communicationTime = communicationTime / MPIcontrol.nprocs;


    double communicationTimetot;
    MPI_Reduce(&mygraph->commtimetot, &communicationTimetot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    communicationTimetot = communicationTimetot / MPIcontrol.nprocs;

    double cutTimetot;
    MPI_Reduce(&mygraph->cuttime, &cutTimetot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    cutTimetot = cutTimetot / MPIcontrol.nprocs;

/*
    double storageTimetot;
    MPI_Reduce(&mygraph->storagetime, &storageTimetot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    storageTimetot = storageTimetot / MPIcontrol.nprocs;
*/
/*
    double elapsedTime;
    MPI_Reduce(&MPIcontrol.elapsedTime, &elapsedTime, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    elapsedTime = elapsedTime / MPIcontrol.nprocs;

    double locateTime;
    MPI_Reduce(&MPIcontrol.locateTime, &locateTime, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    locateTime = locateTime / MPIcontrol.nprocs;

    double readTime;
    MPI_Reduce(&MPIcontrol.locateTime, &readTime, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    readTime = readTime / MPIcontrol.nprocs;
*/
    if(MPIcontrol.rank==0)
    {
// 	fprintf(FP, "Time spend in I/O part in MPIcontrol %.5f\n", MPIcontrol.elapsedTime/(float)CLOCKS_PER_SEC);
//	fprintf(FP, "Time spend in locate function %.5f\n", MPIcontrol.locateTime/(float)CLOCKS_PER_SEC);
//	fprintf(FP, "Time spend in read function %.5f\n", MPIcontrol.readTime/(float)CLOCKS_PER_SEC);
 
	sprintf(message,"Time spend in I/O part in MPIcontrol %.5f\n", MPIcontrol.elapsedTime/(float)CLOCKS_PER_SEC);
 	MPIcontrol.print(message, FP);	
	sprintf(message, "Time spend in locate function %.5f\n", MPIcontrol.locateTime/(float)CLOCKS_PER_SEC);
	MPIcontrol.print(message,FP);
	sprintf(message,"Time spend in read function %.5f\n", MPIcontrol.readTime/(float)CLOCKS_PER_SEC);
	MPIcontrol.print(message,FP);
    }

    delete mygraph;
//    delete sequences;

    MPI_Barrier(MPI_COMM_WORLD);
    if(MPIcontrol.rank==0)
    {
	t2 = clock();
	sprintf(message, "Time spend in Constructing Graph %.5f (cutTimetot=%.5f, communication time=%.5f, communicationtimetot=%.5f)\n", (t2-t1)/(float)CLOCKS_PER_SEC, cutTimetot/(float)CLOCKS_PER_SEC, communicationTime/(float)CLOCKS_PER_SEC, communicationTimetot/(float)CLOCKS_PER_SEC);
	MPIcontrol.print(message, FP);
//	fprintf(FP, "Time spend in Constructing Graph %.5f (communication time=%.5f)\n", (t2-t1)/(float)CLOCKS_PER_SEC, communicationTime/(float)CLOCKS_PER_SEC);
    }

    sprintf(message, "Construct nodeGraph finished");
    MPIcontrol.print(message);

    MPI_Barrier(MPI_COMM_WORLD);

    //nodeGraph->arcFrequency();	
    nodeGraph->cutoffGraph(&MPIcontrol, parameters.cutoffThreshold);

    sprintf(message, "Cutoff Graph finished");
    MPIcontrol.print(message);
    
    MPI_Barrier(MPI_COMM_WORLD);
    if(MPIcontrol.rank==0)
    {
	t3 = clock();
	sprintf(message, "Time spend in cutoffGraph %.5f\n", (t3-t2)/(float)CLOCKS_PER_SEC);
	MPIcontrol.print(message, FP);
//	fprintf(FP, "Time spend in cutoffGraph %.5f\n", (t3-t2)/(float)CLOCKS_PER_SEC);
    }

    nodeGraph->buildKmoleculeGraph(&MPIcontrol, &parameters);

//    nodeGraph->checkDistNodeGraph(); 
     nodeGraph->graphStatistics(&parameters, &MPIcontrol);
//    nodeGraph->printDistNodeGraph(&parameters, &MPIcontrol);


    MPI_Barrier(MPI_COMM_WORLD);
    if(MPIcontrol.rank==0)
    {
	t4 = clock();
	sprintf(message, "Time spend in Constructing DistGraph %.5f\n", (t4-t3)/(float)CLOCKS_PER_SEC);
	MPIcontrol.print(message, FP);

//	fprintf(FP, "Time spend in Constructing DistGraph %.5f\n", (t4-t3)/(float)CLOCKS_PER_SEC);
    }

    //edge merging 
    if(parameters.kmerGraphFlag==0)   nodeGraph->simplifyDistNodeGraph(&parameters, &MPIcontrol);

    double worktotal, timetotal;
    MPI_Reduce(&nodeGraph->totalWorkTime, &worktotal, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&nodeGraph->totalTime, &timetotal, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);	

    MPI_Barrier(MPI_COMM_WORLD);
    if(MPIcontrol.rank==0)
    {
	t5 = clock();
	sprintf(message, "Time spend in Simplify DistGraph %.5f (workTime=%.5f(%.5f), totalTime=%.5f(%.5f))\n", (t5-t4)/(float)CLOCKS_PER_SEC, nodeGraph->totalWorkTime/(float)CLOCKS_PER_SEC, worktotal/(float)CLOCKS_PER_SEC, nodeGraph->totalTime/(float)CLOCKS_PER_SEC, timetotal/(float)CLOCKS_PER_SEC);
	MPIcontrol.print(message,FP);
//	fprintf(FP, "Time spend in Simplify DistGraph %.5f (workTime=%.5f(%.5f), totalTime=%.5f(%.5f))\n", (t5-t4)/(float)CLOCKS_PER_SEC, nodeGraph->totalWorkTime/(float)CLOCKS_PER_SEC, worktotal/(float)CLOCKS_PER_SEC, nodeGraph->totalTime/(float)CLOCKS_PER_SEC, timetotal/(float)CLOCKS_PER_SEC);
    }

    if(parameters.performanceFlag)
    {
    	nodeGraph->stringContigs(&parameters, &MPIcontrol);
    	nodeGraph->nodes.clear();

    	MPI_Barrier(MPI_COMM_WORLD);
    	if(MPIcontrol.rank==0)
    	{
		t1 = clock();
		fprintf(FP, "Time spend in printDistGraph %.5f\n", (t1-t5)/(float)CLOCKS_PER_SEC);
    	}

    	MPIcontrol.finalize();
    	if(MPIcontrol.rank==0)
    	{
		t5 = clock();
		fprintf(FP, "Time spend in All step %.5f\n", (t5)/(float)CLOCKS_PER_SEC);
		fclose(FP);
    	}
    	return (0);
    } 

//  nodeGraph->printDistNodeGraph(&parameters, &MPIcontrol);
    char *rBuf = nodeGraph->gatherDistNodeGraph(&parameters, &MPIcontrol);
    nodeGraph->nodes.clear();

    if(parameters.distGraphFlag && parameters.kmerGraphFlag==0)	
    {
	FILE *fp = fopen(parameters.graphPath,"w");
	fprintf(fp, "%s\n", rBuf);
	fclose(fp);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if(MPIcontrol.rank==0)
    {
	t1 = clock();
	fprintf(FP, "Time spend in printDistGraph %.5f\n", (t1-t5)/(float)CLOCKS_PER_SEC);
    }


    if(MPIcontrol.rank==0)
    { 
	nodeGraph->recvBufDistNodeGraph(rBuf, &parameters, &MPIcontrol);
	if(parameters.kmerGraphFlag)  
	{
		nodeGraph->printKmerGraph(&parameters, &MPIcontrol);
		nodeGraph->nodes.clear();
		MPIcontrol.finalize();
		return 0;
        }

//    	nodeGraph->checkDistNodeGraph();
//	nodeGraph->checkLocalNodeGraph(&parameters);
	
	t6 = clock();
	fprintf(FP, "Time spend in reading DistGraph  %.5f\n", (t6-t1)/(float)CLOCKS_PER_SEC);

	nodeGraph->masterSimplifyNodeGraph(&parameters, &MPIcontrol); 
    	printf("proc:%d masterSimplifyNodeGraph  end\n", MPIcontrol.rank);

        nodeGraph->printContigs(&parameters, &MPIcontrol);
	
        for(int i=0;i<4;i++)
	{	
		nodeGraph->masterTipsRemoval(&parameters, &MPIcontrol);
    		printf("proc:%d masterTipsRemoval end\n", MPIcontrol.rank);
		nodeGraph->masterSimplifyNodeGraph(&parameters, &MPIcontrol); 
    		printf("proc:%d masterSimplifyNodeGraph  end\n", MPIcontrol.rank);
	
		nodeGraph->masterLowCoverageEdgeRemoval(&parameters, &MPIcontrol);
    		printf("proc:%d masterLowCoverageEdgeRemoval  end\n", MPIcontrol.rank);

		nodeGraph->masterSimplifyNodeGraph(&parameters, &MPIcontrol); 
    		printf("proc:%d masterSimplifyNodeGraph  end\n", MPIcontrol.rank);

                nodeGraph->masterLoopBubbleRemoval(&parameters, &MPIcontrol);
                printf("proc:%d LoopBubbleRemoval end\n", MPIcontrol.rank);
                nodeGraph->masterMultipleEdgeBubbleRemoval(&parameters, &MPIcontrol);
                printf("proc:%d MultipleEdgeRemoval end\n", MPIcontrol.rank);
                nodeGraph->masterSimplifyNodeGraph(&parameters, &MPIcontrol);
                printf("proc:%d masterSimplifyNodeGraph end\n", MPIcontrol.rank);
                nodeGraph->masterGraphStatistic(&parameters, &MPIcontrol);
        }

        nodeGraph->masterRemoveCrossEdge(&parameters, &MPIcontrol);
        printf("proc:%d masterRemoveCrossEdge end\n", MPIcontrol.rank);

        t7 = clock();
        fprintf(FP, "Time spend in Contig Extension %.5f\n", (t7-t1)/(float)CLOCKS_PER_SEC);

        nodeGraph->masterPrintContigs(&parameters, &MPIcontrol);
        if(parameters.JungGraphFlag)    nodeGraph->printJungGraph(&parameters, &MPIcontrol);
        nodeGraph->nodes.clear();
    }

    if(MPIcontrol.rank==0)
    {
        t5 = clock();
//	sprintf(message, "Time spend in All step %.5f\n", (t5)/(float)CLOCKS_PER_SEC);
//	MPIcontrol.print(message, FP);
        fprintf(FP, "Time spend in All step %.5f\n", (t5)/(float)CLOCKS_PER_SEC);
        fclose(FP);
    }
    
    MPIcontrol.finalize();
    return (0);
} 


