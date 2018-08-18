#include "mympi.h"
#include "sequence.h"

void parameter::getopt(int argc, char **argv)
{
    this->fastaPath[0] = 0;
    this->outputPath[0] = 0;
    this->hashLength=23;
    this->cutoffThreshold = 5;
	
    this->kmerGraphFlag=0;
    this->JungGraphFlag=0;
    this->distGraphFlag=0;
    this->performanceFlag=0;

    char comStr[100];
    int  ret;
    for(int i=1;i<argc;i++)
    {
	sscanf(argv[i], "%s", comStr);
 	if(comStr[0] != '-')  
	{
	    printf("Usage: %s -k <kmerlength> -c <cutoff Threshod> -i <input Fasta data> -o <output data directory>\n",argv[0]);
//	    MPI_Abort(MPI_COMM_WORLD,1);		
	    exit(0);
	}	
	switch(comStr[1])
	{
	    case 'h':
	    	printf("Usage: %s -k <kmerlength> -c <cutoff Threshod> -i <input Fasta data> -o <output data directory>\n",argv[0]);
		exit(0);
	    case 'v':
		printf("SWAP-Assembler version 0.2\n");
		exit(0);
	    case 'k':
		if(i+1 == argc) 
		{
		    printf("Error: -k needs a value for the kmer size\n");
	    	    printf("Usage: %s -k <kmerlength> -c <cutoff Threshod> -i <input Fasta data> -o <output data directory>\n",argv[0]);
		    exit(0);
		}
		ret = sscanf(argv[i+1], "%d", &this->hashLength);
		i++;
		if(ret!=1 || this->hashLength<19 || this->hashLength>31 || this->hashLength%2!=1)	
		{
		    printf("Error: -k needs a value (odd number between 18 and 32) for the kmer size\n");
		    printf("Usage: %s -k <kmerlength> -c <cutoff Threshod> -i <input Fasta data> -o <output data directory>\n",argv[0]);
	            exit(0);
		}
		break;	
	    case 'c':
		if(i+1 == argc) 
		{
		    printf("Error: -c needs a value for the cutoff threshold\n");
	    	    printf("Usage: %s -k <kmerlength> -c <cutoff Threshod> -i <input Fasta data> -o <output data directory>\n",argv[0]);
		    exit(0);
		}
		ret = sscanf(argv[i+1], "%d", &this->cutoffThreshold);
		i++;
		if(ret!=1)	
		{
		    printf("Error: -c needs a value for the cutoff threshold (it is suggested to be 3~10% of the coverage of dataset)\n");
		    printf("Usage: %s -k <kmerlength> -c <cutoff Threshod> -i <input Fasta data> -o <output data directory>\n",argv[0]);
	            exit(0);
		}
		break;
            case 'i':
		if(i+1 == argc) 
		{
		    printf("Error: -i needs a path for the file in fasta format\n");
	    	    printf("Usage: %s -k <kmerlength> -c <cutoff Threshod> -i <input Fasta data> -o <output data directory>\n",argv[0]);
		    exit(0);
		}
		ret = sscanf(argv[i+1], "%s", this->fastaPath);
		i++;
		if(ret==0)
		{
		    printf("Error: -i needs a path for the file in fasta format\n");
	    	    printf("Usage: %s -k <kmerlength> -c <cutoff Threshod> -i <input Fasta data> -o <output data directory>\n",argv[0]);
		    exit(0);
		}
		break;
            case 'o':
		if(i+1 == argc) 
		{
		    printf("Error: -o needs a directory name used for the output data\n");
	    	    printf("Usage: %s -k <kmerlength> -c <cutoff Threshod> -i <input Fasta data> -o <output data directory>\n",argv[0]);
		    exit(0);
		}
		ret = sscanf(argv[i+1], "%s", this->outputPath);
		i++;
		if(ret==0)
		{
		    printf("Error: -o needs a directory path for the output data\n");
	    	    printf("Usage: %s -k <kmerlength> -c <cutoff Threshod> -i <input Fasta data> -o <output data directory>\n",argv[0]);
		    exit(0);
		}
		break;
            case 's':
		this->kmerGraphFlag = 1;
		break;
	    case 'j':
		this->JungGraphFlag = 1;
		break;
	    case 'd':
		this->distGraphFlag = 1;
		break;
	    case 'p':
		this->performanceFlag = 1;
		break;
	}	
    }
    if(this->fastaPath[0]==0 || this->outputPath[0]==0)
    {
        printf("Usage: %s -k <kmerlength> -c <cutoff Threshod> -i <input Fasta data> -o <output data directory>\n",argv[0]);
        exit(0);
    }

}

void parameter::getParameters(int argc, char **argv, MPIEnviroment *MPIcontrol)
{
    
    getopt(argc, argv);

    if(MPIcontrol->rank==0)
    {    
    	int stat = mkdir(this->outputPath, S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);
    	if(stat==-1)	
    	{
		printf("Error: Creating directory %s\n", this->outputPath);
		exit(0);
    	} 
       

    }


    strcpy(this->contigsPath, this->outputPath);
    strcpy(this->LogPath, this->outputPath);
    strcpy(this->graphPath, this->outputPath);
    strcpy(this->masterContigPath,this->outputPath);
    strcpy(this->JungGraph_arc, this->outputPath);
    strcpy(this->JungGraph_mul, this->outputPath);
    strcpy(this->arcFrequency, this->outputPath);
    strcpy(this->kmerGraph, this->outputPath);

    strcat(this->contigsPath, "/noCEcontigs.fasta");
    strcat(this->LogPath, "/logtime.txt"); 
    strcat(this->graphPath, "/contigGraph.txt");
    strcat(this->masterContigPath,"/CEContig.fasta");
    strcat(this->JungGraph_arc, "/JungGraph_arc.txt");
    strcat(this->JungGraph_mul, "/JungGraph_mul.txt");
    strcat(this->arcFrequency, "/arcFrequency.txt");    
    strcat(this->kmerGraph, "/kmerGraph.txt");

    if(MPIcontrol->rank==0)
    {    printf("Runing command: %s -k %d -c %d -i %s -o %s\n", argv[0], this->hashLength, this->cutoffThreshold, this->fastaPath, this->outputPath);
         printf("contigsPath=%s, logPath=%s, graphPath=%s\n", this->contigsPath, this->LogPath, this->graphPath);
    }
//    strcpy(this->LogPath, "logtime.txt");
//    strcpy(this->graphPath, "distGraph.txt");

    char ch0='A', ch1='C', ch2='G', ch3='T';
    this->nucleotideValue[(int)ch0] = 0;
    this->nucleotideValue[(int)ch1] = 1;
    this->nucleotideValue[(int)ch2] = 2;
    this->nucleotideValue[(int)ch3] = 3;

    this->nucleotideArray[0] = 'A';
    this->nucleotideArray[1] = 'C';
    this->nucleotideArray[2] = 'G';
    this->nucleotideArray[3] = 'T';

    this->nucleotideReverse[(int)ch0] = 'T';
    this->nucleotideReverse[(int)ch1] = 'G';
    this->nucleotideReverse[(int)ch2] = 'C';
    this->nucleotideReverse[(int)ch3] = 'A';

    this->MASK = ~(3ull<<(2*hashLength));
    this->filterThreshold = Filter_Threshold;
 //   this->arcThreshold    = 5;
}

sequence::~sequence()
{
    for(int i=0;i<this->readCount;i++)
	delete this->reads[i];	
    delete this->reads;
}

void sequence::getSequences(parameter *p, MPIEnviroment *MPIcontrol)
{
    char buf[BUF_SIZE*2];
    char tmpBuf[BUF_SIZE*2];
    char tmp[BUF_SIZE];

    MPIcontrol->File_open(p->fastaPath);
    MPIcontrol->File_locate();

    this->readCount=0;
    while(MPIcontrol->File_read(buf,BUF_SIZE) != 0)	
    {	
	for(int i=0;i<BUF_SIZE;i++)	
	    if(buf[i]=='>') this->readCount++;
    }

    MPIcontrol->File_locate();       

    this->reads = new char* [this->readCount];
    assert(this->reads!=NULL);

    unsigned long long totNumReads = this->readCount;
    this->readCount = 0;

    unsigned long long readStart=0;
    unsigned long long len;
    while((len=MPIcontrol->File_read(buf+readStart, BUF_SIZE)) != 0)
    {
	if(DEBUG)	MPIcontrol->print(buf);	
	buf[readStart+len] = 0;
	unsigned long long i=0, local;
	while(i<readStart+len)
	{
	    local = i;
	    if(buf[i]!='>')	{
		printf("Error at getSequences %100s\n", buf+i);
		exit(0);
	    }	

	    while(buf[i]!='\n'&& i<readStart+len) i++;
	    if(i==readStart+len)	break;			

	    int start_pos=i+1;
	    while(buf[i]!='>' && i<readStart+len) i++;
	    if(i==readStart+len)	break;			

	    //between start_pos and i, there is one read
	    this->reads[this->readCount] = new char[i-start_pos];
	    assert(this->reads[this->readCount]!=NULL);

	    int index = 0;
	    for(int j=start_pos; j<i;j++)
	    {
		if(buf[j]=='A'||buf[j]=='C'||buf[j]=='G'||buf[j]=='T'||
			buf[j]=='a'||buf[j]=='c'||buf[j]=='g'||buf[j]=='t')
		this->reads[this->readCount][index++]=buf[j];	
	    }
	    this->reads[this->readCount][index]=0;
/*	    if(index != 36)  {
		sprintf(tmp,"Get %d reads: %s", readCount, reads[readCount]);		
		MPIcontrol->print(tmp);
		exit(0);
	    }
*/
	    this->readCount++;

	    if(this->readCount%100000==0)
	    {
		sprintf(tmp,"Get %d reads: %llu", readCount, totNumReads);		
		MPIcontrol->print(tmp);
	    }
	}
	strncpy(tmpBuf, buf+local, readStart+len - local);
	strncpy(buf, tmpBuf, readStart+len - local);			
	readStart = readStart+len -local;
    }

    //collect the last read
    if(buf[0]=='>')
    {
	int i=0, local;
	while(buf[i]!='\n' && i<readStart)	i++;
	i++;
	this->reads[this->readCount] = new char[readStart-i];
	assert(this->reads[this->readCount] != NULL);	
	int index = 0;
	for(int j=i;j<readStart;j++)
	{
	    if(buf[j]=='A'||buf[j]=='C'||buf[j]=='G'||buf[j]=='T'||
		    buf[j]=='a'||buf[j]=='c'||buf[j]=='g'||buf[j]=='t')
		this->reads[this->readCount][index++] = buf[j];
	}
	this->reads[this->readCount][index] = 0;
	if(DEBUG) {	
	    sprintf(tmp, "Get %d reads: %s", readCount, reads[readCount]);
	    MPIcontrol->print(tmp);
	}
	this->readCount++;		
    }
    MPIcontrol->File_close();

    //if(DEBUG) 
    {
	sprintf(tmp, "Number of Reads %d", readCount);
	MPIcontrol->print(tmp);
    }

    MPI_Reduce(&readCount, &totReads, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);		

    if(MPIcontrol->rank==0) 
    {
	sprintf(tmp, "Total Number of Reads %d", totReads);
	MPIcontrol->print(tmp);
    }
}

/*
   int main (int argc, char *argv[])
   {
   int i, size, namelen, lgsize;
   MPIEnviroment MPIcontrol;
   MPIcontrol.init(argc, argv);

   parameter parameters;
   sequence sequences;
   parameters.getParameters(argc, argv);	
   sequences.getSequences(&parameters, &MPIcontrol);
   MPIcontrol.finalize();
   return (0);
   }
*/
