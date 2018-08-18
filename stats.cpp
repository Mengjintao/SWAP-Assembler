#include <string>
#include <algorithm>
#include <vector>
#include <iostream>
#include <stdio.h>
#include <fstream>

using namespace std;
char str[100000000];
char tmp[100000000];

bool Greater(const string& a, const string& b)
{
	return a.length() > b.length();
}


int main(int argc, char *argv[])
{
        char contigsFile[300];
        long long int refSeqLength;

        sscanf(argv[1],"%s", contigsFile);
        sscanf(argv[2],"%lld", &refSeqLength);

	vector<string> vs;
	ifstream ifs;
	ifs.open (contigsFile, std::ifstream::in);

	string tmp;
	tmp.clear();
	while(ifs.getline(str, 100000000))
	{
		if(str[0]=='>')	
		{
			if(tmp.length() < 200)	continue;
			vs.push_back(tmp);
			tmp.clear();
			continue;
		}
		tmp += string(str);	
	}
	if(tmp.length()>=200)	vs.push_back(tmp);
	
	ifs.close();
	printf("Input %d contigs\n", vs.size());
	
	sort( vs.begin(), vs.end(), Greater );

	printf("Sort end\n");	

	long long int sum = 0, tot = 0;
	long long int N50_abs = 0, N50_self=0;

	for(int i=0;i<vs.size();i++)
	{
		tot += vs[i].length();
	}

	sum = 0;
	for(int i=0;i<vs.size();i++)
	{
		sum += vs[i].length();
		
		if(sum>tot/2)
		{
			N50_self = vs[i].length(); 
			break;	
		}
	}	

	sum = 0;
	for(int i=0;i<vs.size();i++)
	{
		sum += vs[i].length();
		
		if(sum>refSeqLength/2)
		{
			N50_abs = vs[i].length(); 
			break;	
		}
	}

	printf("ContigNum=%d, tot = %lld, N50_self=%lld, N50_abs=%lld\n", vs.size(), tot, N50_self, N50_abs);
	printf("Max = %d, mean = %d\n", vs[0].length(), vs[vs.size()/2].length());
}	
