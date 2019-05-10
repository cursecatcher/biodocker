/***************************************************************************
 *   Copyright (C) 2013 by Marco Beccuti                                   *
 *   beccuti@di.unito.it                                                   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#ifndef __CIO_H__
	#define __CIO_H__
	#include "classIO.h"
#endif

#ifndef __MybitC_CPP__
	#define __MybitC_CPP__
    #include "classMybitset.cpp"
#endif

#ifndef _hCh_CPP__
	#define __hCh_CPP__
	#include "classHashChaining.cpp"
#endif

extern  int* NotFReads;
namespace IOFn
{
  
 
  
  /**************************************************************/
  /* NAME :  Class IOF*/
  /* DESCRIPTION : Constructor*/
  /**************************************************************/
  IOF::IOF(std::string ReadFname1,std::string fileExt1,std::string ReadFname2,std::string fileExt2,const std::string& OFname,const int& window,const int& files1,const int& files2,const unsigned int& numker){

    this->ReadFname1=ReadFname1;
    this->fileExt1=fileExt1;
    this->ReadFname2=ReadFname2;
    this->fileExt2=fileExt2;
    this->OFname=OFname;
    this->window=window;
    this->files1=files1;
    this->files2=files2;
    this->numker=numker;
    count=(unsigned long long *)malloc(files1*sizeof(unsigned long long));
    memset(count,0, files1*sizeof(unsigned long long));
  };
  
  /**************************************************************/
/* NAME :  Class IOF AllocHash() */
/* DESCRIPTION : hash table allocator */
/**************************************************************/  
void IOF::AllocHash(const unsigned int& hash,const unsigned int&list){
  //initialized HASH TABLE
  if (hash>0)
    sh.UpdateHash(hash,list);
  else
  {
    std::cerr<<"\n\nERROR: hash size is"<<hash<<"\n";
    exit(EXIT_FAILURE);
    
  }  
}
  
/**************************************************************/
/* NAME :  Class IOF  ResetReads() */
/* DESCRIPTION : it resets the read counter*/
/**************************************************************/    
void IOF::ResetReads(void)
{
  reads = 0;
}
 
/**************************************************************/
/* NAME :  Class IOF  Read() */
/* DESCRIPTION : hash table creator by pool files */
/**************************************************************/  
void IOF::Read(int l, int u){
  
  using namespace general;
  clock_t startGlobal,endGlobal;
  char buffer[MAXSIZE];
  ifstream in;
  //reads=0;

  disc=0;
  int wrongreads=0;
  //For each pool files
//  for (int i=0;i<files1;i++)  
  //For all pools in the input interval
  for (int i=l;i<=u;i++)
  {
    startGlobal=clock();
    ostringstream of;
    of<<ReadFname1<<i<<"."<<fileExt1;
    
    in.open(of.str().c_str(),ifstream::in);
    if(!in) 
    {
      cerr << "\n*****Error opening input file "<<of.str().c_str() <<" *****" << endl;
      exit(EXIT_FAILURE);
      
    }
    
    int hh=0;
    //For each read in the pool file
    while (!in.eof())
    {
      buffer[0]='\0';
      //Name
      in.getline(buffer,MAXSIZE);
      //Sequence
      buffer[0]='\0';
      in.getline(buffer,MAXSIZE);
      int num=in.gcount();
      if (buffer[num-1]!='\0')
      {
	buffer[num]='\0';
	num++;
	
      }
      if ((num>=window+1)&&((buffer[0]=='A')||(buffer[0]=='G')||(buffer[0]=='C')||(buffer[0]=='T')||(buffer[0]=='N')))
      {
	hh++;
	//insert window in the hash table
	insert(buffer,num,i);
	if(reads%10000000==0)
	{
	//#pragma omp critical
        //{
	  cout<<"Reads: "<<reads<<" wrong reads: "<<wrongreads<<" windows: "<<count[i]<<" stored windows:  "<<sh.size()<<" rejected  windows: "<<disc<<endl;
	  //sh.print();
	  sh.info();
	//}  
	}
	 in.getline(buffer,MAXSIZE);
	 in.getline(buffer,MAXSIZE);
      }
      else
      {
	if ((buffer[0]!='\0'))
	{
	  //cerr<<"Error input format file "<<of.str().c_str()<<" buffer:"<<buffer<<endl;
	  if ((buffer[0]=='A')||(buffer[0]=='G')||(buffer[0]=='C')||(buffer[0]=='T')||(buffer[0]=='N'))
	  {
	    wrongreads++;
	    in.getline(buffer,MAXSIZE);
	    in.getline(buffer,MAXSIZE);
	    
	  }
	  else
	  {
	    exit(EXIT_FAILURE);
	    
	  }

	}
      }
      
    }
    in.close();
    endGlobal=clock();
    //#pragma omp critical
   //{
    cout<<"\n\n\tReads: "<<reads<<" wrong reads: "<<wrongreads<<" windows: "<<count[i]<<" stored windows: "<<sh.size()<<".\n\tTime to read pool "<<i<<":"<<((double)(endGlobal-startGlobal))/CLOCKS_PER_SEC<<"s."<<endl;
    //}
  }
#if DEBUG   
 sh.print();
 sh.info();
#endif
 
}
 
 
/**************************************************************/
/* NAME :  Class IOF insert() */
/* DESCRIPTION : encoding and inserting windows */
/**************************************************************/				
void IOF::insert(const char buffer[],const int& num,const int& pool){
  
// #pragma omp critical
//      { 
	reads++;
//	}
  class mybitset b,d;
  for(int k=0; k<=num-(window)-1; k++)
  {
    int z=0;
    int j=window*2-1;
    bool findN=false;
    for (int c = 0; c<window; c++) 
    {
      switch (buffer[c+k])
      {
	case 'A':
	  b.set(z,0);
	  z++;
	  b.set(z,0);
	  z++;
	  d.set(j,1);
	  j--;
	  d.set(j,1);
	  j--;
	  break;
	case 'C':
	  b.set(z,0);
	  z++;
	  b.set(z,1);
	  z++;
	  d.set(j,0);
	  j--;
	  d.set(j,1);
	  j--;
	  break;
	case 'G':
	  b.set(z,1);
	  z++;
	  b.set(z,0);
	  z++;
	  d.set(j,1);
	  j--;
	  d.set(j,0);
	  j--;
	  break;
	case 'T':
	  b.set(z,1);
	  z++;
	  b.set(z,1);
	  z++;
	  d.set(j,0);
	  j--;
	  d.set(j,0);
	  j--;
	break;
	default:
	  findN=true;
	  c=window;
      }
      
    }

    if (!findN)
    {
      count[pool]++;
      //cout<<b<<"\n"<<d<<endl<<endl;
      if (d<b)
      {
//#if SEARCHNORMAL
	if (!sh.searchIncFreq(d))
//#else
//	  if (!sh.searchAndmoveSetPool(d,pool))
//#endif
	  {
	    sh.insert(d);
	  }
      }
      else
      {
//#if SEARCHNORMAL
	  if (!sh.searchIncFreq(b))
//#else
//	    if (!sh.searchAndmoveSetPool(b,pool))
//#endif
	    {	
	    sh.insert(b);
	    }
      }
      
    }
    else
//       #pragma omp critical
//        { 
	disc++;
//	}
    
  }
  
}


/**************************************************************/
/* NAME :  Class IOF ReadB()*/
/* DESCRIPTION : hash table reader from binary files*/
/**************************************************************/
bool IOF::ReadB(){
  
  FILE * pFile;
  ostringstream tmp;
  tmp<<ReadFname1<<window;
  pFile = fopen (std::string(tmp.str()+string(".count")).c_str(), "r" );
  if (pFile==NULL)
  {
    cerr<<"\n*****Error opening input file .count *****" << endl;
    cout<<"\nHash Table binary encoding not found!"<<endl;
    return false;
    
  }
  

  unsigned int filetmp=0;
  
  
  if (fread (&filetmp, sizeof(unsigned int),1, pFile )!=1)	
  {
  cerr<<"Reading error"<<endl; 
  exit (EXIT_FAILURE);
  }
  
  if (filetmp!=(unsigned)files1)
  {
    cerr<<"\n*****Error you read a different number of pools ("<<filetmp<<") *****" << endl;
    return false;
    
  }
  
  if (fread (count, sizeof(unsigned long long),files1, pFile )!=(unsigned)files1)
  {
    cerr<<"Reading error"<<endl; 
    exit (EXIT_FAILURE);
    
  }
  
  bool out1=sh.Read(tmp.str());
  
  return ((out1));
  
}


/**************************************************************/
/* NAME :  Class IOF WriteB()*/
/* DESCRIPTION : hash table reader from binary files*/
/**************************************************************/
void IOF::WriteB(){
  
  FILE * pFile;
  //string tmp=ReadFname1+window+".count";
  ostringstream tmp;
  tmp<<ReadFname1<<window;
  pFile = fopen (std::string(tmp.str()+string(".count")).c_str(), "w" );
  if (pFile==NULL)
  {
    cerr<<"\n*****Error opening output file .count *****" << endl;
    exit(EXIT_FAILURE);
    
  }
  
  //write pool number
  fwrite (&files1, sizeof(int), 1, pFile );
  //write number of windows for pools
  fwrite (count, sizeof(unsigned long long), files1, pFile );
  fclose (pFile);
	
  cout<<"Windows frequencies for pool saved in:\n\t"<<tmp.str()<<endl;
  
  //write pool hash table
  sh.Write(tmp.str());
  
};

 
/**************************************************************/
/* NAME :  Class IOF Search() */
/* DESCRIPTION : de-convolution step*/
/**************************************************************/
#if !PAIREDEND
void IOF::Search(int l,int u){
  
  using namespace general;
  HASHCH<mybitset> rd;
  rd.UpdateHash(READSIZE-window,READSIZE-window);
  char buffer[MAXSIZE],name[MAXSIZE];//,name2[MAXSIZE];
#if DEBUGOUTPUT  
  bool vetKmer[MAXSIZE];
#endif  
  ifstream in;
  ofstream out;
  //reads=0;
  int num=0,wrongreads=0;
  clock_t startGlobal,endGlobal;
  
  

  //For all pools in the input interval
  for (int i=l;i<=u;i++)
  {
    startGlobal=clock();
    ostringstream of,of1;


    of<<ReadFname2<<i<<"."<<fileExt2;
    
    in.open(of.str().c_str(),ofstream::in);
    if(!in) 
    {
      cerr << "\n*****Error opening input file "<<of.str().c_str()<<" *****" << endl;
      exit(EXIT_FAILURE);
    }
#if !SIGNATURE    
    of1<<OFname<<i;
    out.open(of1.str().c_str(),ofstream::out);
    if (!out)
    {
      cerr << "\n*****Error opening ouput file "<< of1.str().c_str()<<" *****" << endl;
      exit(EXIT_FAILURE);
    }    
    int dec=0;
    int pdec=0;
#endif 

#if SIMILARITY
    of1.str("");
    of1<<OFname<<"Sim"<<i;
    ofstream outsim(of1.str().c_str(),ofstream::out);
    if (!outsim)
    {
      cerr << "\n*****Error opening ouput file "<< of1.str().c_str()<<" *****" << endl;
      exit(EXIT_FAILURE);
    }
#endif

    
    while (!in.eof())
    {
      buffer[0]='\0';
      in.getline(buffer,MAXSIZE);
      num=in.gcount();
      if (buffer[num-1]!='\0')
      {
	buffer[num]='\0';
	num++;
	
      }
      if ((buffer[0]=='\0')||((num>0)&&((buffer[0]=='@')||(buffer[0]=='>'))))
	{
	  strcpy(name,buffer);
	}
      else
      {
	cerr<<"Error input format file "<<of.str().c_str()<<endl;
	exit(0);
      }
      buffer[0]='\0';
      in.getline(buffer,MAXSIZE);
      num=in.gcount();
      if (buffer[num-1]!='\0')
      {
	buffer[num]='\0';
	num++;
	
      }
      if ((num>=window+1)&&((buffer[0]=='A')||(buffer[0]=='G')||(buffer[0]=='C')||(buffer[0]=='T')||(buffer[0]=='N')))
      {

#if SIGNATURE
       searchEachXSignature(buffer,num);	   
       //+
      in.getline(buffer,MAXSIZE);	  
      //Quality values 
      in.getline(buffer,MAXSIZE);
#else	
       unsigned int countOk=0,countKo=0;
  #if !DEBUGOUTPUT	
	searchAll(rd,buffer,num,out,countOk,countKo);
  #else
	searchEach(vetKmer,buffer,num,out,countOk,countKo);
  #endif	
	
if (countOk >= numker)
	{
	  dec++;
  #if SIMILARITY
	 outsim<<name<<endl<<buffer<<endl;
  #endif	  
	 //+
	 in.getline(buffer,MAXSIZE);
  #if SIMILARITY
	 outsim<<buffer<<endl;
  #endif	
	 //Quality values 
	 in.getline(buffer,MAXSIZE);
  #if SIMILARITY
	 outsim<<buffer<<endl;
  #endif		 
	}
#endif	
else
    {
    in.getline(buffer,MAXSIZE); //+
    in.getline(buffer,MAXSIZE); //Quality     
    }
	if(reads % 10000000 == 0)
	{
	#pragma omp critical
        {//critical section for hash update
	  cout<<"Reads "<<reads<<" wrong reads"<<wrongreads<<endl;
	}  
	//sh.print();
	}		
      }
      else
      {
      if ((buffer[0]!='\0'))
	{
	  //cerr<<"Error input format file "<<of.str().c_str()<<" buffer:"<<buffer<<endl;
	  if ((buffer[0]=='A')||(buffer[0]=='G')||(buffer[0]=='C')||(buffer[0]=='T')||(buffer[0]=='N'))
	  {
	    wrongreads++;
	    in.getline(buffer,MAXSIZE);
	    in.getline(buffer,MAXSIZE);  
	  }
	  else
	  {
	    exit(EXIT_FAILURE);
	    
	  }

	}
      }

    } 
    
    in.close();
#if !SIGNATURE    
    out.close();
#endif
    
#if SIMILARITY
    outsim.close();
#endif    
    endGlobal=clock();
#pragma omp critical
    {//critical section for hash update
    #if !SIGNATURE 
    cout<<"\n\n\tFound reads: "<<dec<<" not found reads: "<<pdec<<" wrong reads"<<wrongreads<<"\n";
    NotFReads[i]=pdec;
    #endif
    cout<<"\n\tTime to search overlappings in "<<i<<": "<<((double)(endGlobal-startGlobal))/CLOCKS_PER_SEC<<"s."<<endl<<endl;
    }
  }
#if SIGNATURE    
    out.open(string(OFname+string(".signature")).c_str(),ofstream::out);
    if (!out)
    {
      cerr << "\n*****Error opening ouput file "<< OFname<<".signature *****" << endl;
      exit(EXIT_FAILURE);
    }
    SignaturePrint(out);
    out.close();
#endif 
}
#endif

/**************************************************************/
/* NAME :  Class IOF Search() */
/* DESCRIPTION : de-convolution step*/
/**************************************************************/
#if PAIREDEND
void IOF::Search(int l, int u) {
	using namespace general;

	HASHCH<mybitset> rd;

	rd.UpdateHash(READSIZE - window, READSIZE - window);

	char buffer1[MAXSIZE], buffer2[MAXSIZE], name1[MAXSIZE], name2[MAXSIZE];

#if DEBUGOUTPUT
	bool vetKmer1[MAXSIZE], vetKmer2[MAXSIZE];
#endif

	ifstream in1;
	ifstream in2;
	ofstream out1;
	ofstream out2;
	int num1 = 0, num2 = 0, wrongreads = 0;
	clock_t startGlobal, endGlobal;

	for (int i = l; i <= u; i++) {
		startGlobal = clock();

		ostringstream infn1, infn2, oufn1, oufn2;

		infn1 << ReadFname2 << ".R1." << i << "." << fileExt2;

		in1.open(infn1.str().c_str(), ofstream::in);

		if (!in1) {
			cerr << "\n***** Error opening input file " << infn1.str().c_str() << " *****" << endl;
			exit(EXIT_FAILURE);
		}

		infn2 << ReadFname2 << ".R2." << i << "." << fileExt2;

		in2.open(infn2.str().c_str(), ofstream::in);

		if (!in2) {
			cerr << "\n***** Error opening input file " << infn2.str().c_str() << " *****" << endl;
			exit(EXIT_FAILURE);
		}

		oufn1 << OFname << i << ".R1." << fileExt2;

		out1.open(oufn1.str().c_str(), ofstream::out);

		if (!out1) {
			cerr << "\n***** Error opening ouput file " << oufn1.str().c_str() << " *****" << endl;
			exit(EXIT_FAILURE);
		}

		oufn2 << OFname << i << ".R2." << fileExt2;

		out2.open(oufn2.str().c_str(), ofstream::out);

		if (!out2) {
			cerr << "\n***** Error opening ouput file " << oufn2.str().c_str() << " *****" << endl;
			exit(EXIT_FAILURE);
		}

#if SIMILARITY
		oufn1.str("");

		oufn1 << OFname << "Sim.R1." << i << fileExt2;

		ofstream out1sim(oufn1.str().c_str(), ofstream::out);

		if (!out1sim) {
			cerr << "\n***** Error opening ouput file " << oufn1.str().c_str() << " *****" << endl;
			exit(EXIT_FAILURE);
		}

		oufn2.str("");

		oufn2 << OFname << "Sim.R2." << i << fileExt2;

		ofstream out2sim(oufn2.str().c_str(), ofstream::out);

		if (!out2sim) {
			cerr << "\n***** Error opening ouput file " << oufn2.str().c_str() << " *****" << endl;
			exit(EXIT_FAILURE);
		}
#endif

		int dec = 0, pdec = 0;

/*
 * XXX: assuming the two files have the same lengths.
 */
		while (!in1.eof()) {
			buffer1[0] = '\0';
			in1.getline(buffer1, MAXSIZE);
			num1 = in1.gcount();

			if (buffer1[num1 - 1] != '\0') {
				buffer1[num1] = '\0';
				num1++;
			}

			if ((buffer1[0] == '\0') || ((num1 > 0) && ((buffer1[0] == '@') || (buffer1[0] == '>')))) {
				strcpy(name1, buffer1);
			} else {
				cerr << "Error input format file " << infn1.str().c_str() << endl;
				exit(EXIT_FAILURE);
			}

			buffer1[0] = '\0';
			in1.getline(buffer1, MAXSIZE);
			num1 = in1.gcount();

			if (buffer1[num1 - 1] != '\0') {
				buffer1[num1] = '\0';
				num1++;
			}

			buffer2[0] = '\0';
			in2.getline(buffer2, MAXSIZE);
			num2 = in2.gcount();

			if (buffer2[num2 - 1] != '\0') {
				buffer2[num2] = '\0';
				num2++;
			}

			if ((buffer2[0] == '\0') || ((num2 > 0) && ((buffer2[0] == '@') || (buffer2[0] == '>')))) {
				strcpy(name2, buffer2);
			} else {
				cerr << "Error input format file " << infn2.str().c_str() << endl;
				exit(EXIT_FAILURE);
			}

			buffer2[0] = '\0';
			in2.getline(buffer2, MAXSIZE);
			num2 = in2.gcount();

			if (buffer2[num2 - 1] != '\0') {
				buffer2[num2] = '\0';
				num2++;
			}

			unsigned int countOk1 = 0, countOk2 = 0, countKo1 = 0, countKo2 = 0;

			if ((num1 >= window + 1) && ((buffer1[0] == 'A') || (buffer1[0] == 'G') || (buffer1[0] == 'C') || (buffer1[0] == 'T') || (buffer1[0] == 'N'))) {
#if !DEBUGOUTPUT
				searchAll(rd, buffer1, num1, out1, countOk1, countKo1);
#else
				searchEach(vetKmer1, buffer1, num1, out1, countOk1, countKo1);
#endif
			} else if (buffer1[0] != '\0') {
				if ((buffer1[0] == 'A') || (buffer1[0] == 'G') || (buffer1[0] == 'C') || (buffer1[0] == 'T') || (buffer1[0] == 'N')) {
					wrongreads++;
					in1.getline(buffer1, MAXSIZE);
					in1.getline(buffer1, MAXSIZE);
				} else {
					exit(EXIT_FAILURE);
				}
			}

			if ((countKo1 == 0) && (num2 >= window + 1) && ((buffer2[0] == 'A') || (buffer2[0] == 'G') || (buffer2[0] == 'C') || (buffer2[0] == 'T') || (buffer2[0] == 'N'))) {
#if !DEBUGOUTPUT
				searchAll(rd, buffer2, num2, out2, countOk2, countKo2);
#else
				searchEach(vetKmer2, buffer2, num2, out2, countOk2, countKo2);
#endif
			}

			if (num2 >= window + 1) {
				if ((countKo1 == 0) && ((buffer2[0] == 'A') || (buffer2[0] == 'G') || (buffer2[0] == 'C') || (buffer2[0] == 'T') || (buffer2[0] == 'N'))) {
#if !DEBUGOUTPUT
					searchAll(rd, buffer2, num2, out2, countOk2, countKo2);
#else
					searchEach(vetKmer2, buffer2, num2, out2, countOk2, countKo2);
#endif
				}
			} else if (buffer2[0] != '\0') {
				if ((buffer2[0] == 'A') || (buffer2[0] == 'G') || (buffer2[0] == 'C') || (buffer2[0] == 'T') || (buffer2[0] == 'N')) {
					wrongreads++;
					in2.getline(buffer2, MAXSIZE);
					in2.getline(buffer2, MAXSIZE);
				} else {
					exit(EXIT_FAILURE);
				}
			}

			if ((countKo1 != 0) || (countKo2 != 0)) {
				out1 << name1 << endl << buffer1 << endl;
				out2 << name2 << endl << buffer2 << endl;

				// +
				in1.getline(buffer1, MAXSIZE);
				in2.getline(buffer2, MAXSIZE);

#if !DEBUGOUTPUT
				out1 << buffer1 << endl;
				out2 << buffer2 << endl;
#endif

				// quality values
				in1.getline(buffer1, MAXSIZE);
				in2.getline(buffer2, MAXSIZE);

#if !DEBUGOUTPUT
				out1 << buffer1 << endl;
				out2 << buffer2 << endl;
#else
				out1 << "#K-mers: ";

				for (int i = 0; i <= num1 - window - 1; i++) {
					if (vetKmer1[i])
						out1 << "1";
					else
						out1 << "0";
				}

				out << "\n#(Wrong K-mers): " << countKo1 << endl;

				out2 << "#K-mers: ";

				for (int i = 0; i <= num2 - window - 1; i++) {
					if (vetKmer2[i])
						out2 << "1";
					else
						out2 << "0";
				}

				out2 << "\n#(Wrong K-mers): " << countKo2 << endl;
#endif

				pdec++;
			} else {
#if SIMILARITY
				out1sim << name1 << endl << buffer1 << endl;
				out2sim << name2 << endl << buffer2 << endl;
#endif

				// +
				in1.getline(buffer1, MAXSIZE);
				in2.getline(buffer2, MAXSIZE);

#if SIMILARITY
				out1sim << buffer1 << endl;
				out2sim << buffer2 << endl;
#endif

				// quality values
				in1.getline(buffer1, MAXSIZE);
				in2.getline(buffer2, MAXSIZE);

#if SIMILARITY
				out1sim << buffer1 << endl;
				out2sim << buffer2 << endl;
#endif

				dec++;
			}

			if (reads % 10000000 == 0) {
#pragma omp critical
			{
				cout << "Reads " << reads << " wrong reads " << wrongreads << endl;
			}
			}
		}

		in1.close();
		in2.close();
		out1.close();
		out2.close();

#if SIMILARITY
		out1sim.close();
		out2sim.close();
#endif

		endGlobal = clock();

#pragma omp critical
		{
			cout << "\n\n\tFound reads: " << dec << " not found reads: " << pdec << " wrong reads " << wrongreads << "\n";

			NotFReads[i]=pdec;

			cout << "\n\tTime to search overlappings in " << i <<": " << ((double) (endGlobal - startGlobal)) / CLOCKS_PER_SEC << "s." << endl << endl;
		}
	}
}
#endif


/**************************************************************/
/* NAME :  Class IOF search() */
/* DESCRIPTION : computing window signature. All the k-mer are stored in a hash table*/
/**************************************************************/
  void IOF::searchAll(class HASHCH<mybitset>& rd,const char buffer[],const int& num,ofstream& out, unsigned int& countok,unsigned int& countko){


   #pragma omp critical
        { 
	reads++;
	}
  class mybitset b,d;
  for(int k=0; k<=num-(window)-1; k++)
  {
    int z=0;
    int j=window*2-1;
    bool findN=false;
    for (int c = 0; c<window; c++) 
      {
	switch (buffer[c+k])
	{
	  case 'A':
	    b.set(z,0);
	    z++;
	    b.set(z,0);
	    z++;
	    d.set(j,1);
	    j--;
	    d.set(j,1);
	    j--;
	    break;
	  case 'C':
	    b.set(z,0);
	    z++;
	    b.set(z,1);
	    z++;
	    d.set(j,0);
	    j--;
	    d.set(j,1);
	    j--;
	    break;
	  case 'G':
	    b.set(z,1);
	    z++;
	    b.set(z,0);
	    z++;
	    d.set(j,1);
	    j--;
	    d.set(j,0);
	    j--;
	    break;
	  case 'T':
	    b.set(z,1);
	    z++;
	    b.set(z,1);
	    z++;
	    d.set(j,0);
	    j--;
	    d.set(j,0);
	    j--;
	    break;
	  default:
	    findN=true;
	    c=window;
	  
	}
	
      }

      if (!findN)
      {
	z=0;	
	if (d<b)
	{
	  rd.insert(d);
	  
	}
	else
	{
	  rd.insert(b);
	  
	}
	
      }
    
  }
  if (rd.size()!=0)
    rd.search(sh,countok,countko);
}



/**************************************************************/
/* NAME :  Class IOF search() */
/* DESCRIPTION : computing window signature.The k-mer are search step by step*/
/**************************************************************/
  void IOF::searchEach(bool vet[],const char buffer[],const int& num,ofstream& out, unsigned int& countok,unsigned int& countko){
  
  #pragma omp critical
  {//critical section for hash update
  reads++;
  }
  class mybitset b,d;
  for(int k=0; k<=num-(window)-1; k++)
  {
    int z=0;
    int j=window*2-1;
    bool findN=false;
    for (int c = 0; c<window; c++) 
      {
	switch (buffer[c+k])
	{
	  case 'A':
	    b.set(z,0);
	    z++;
	    b.set(z,0);
	    z++;
	    d.set(j,1);
	    j--;
	    d.set(j,1);
	    j--;
	    break;
	  case 'C':
	    b.set(z,0);
	    z++;
	    b.set(z,1);
	    z++;
	    d.set(j,0);
	    j--;
	    d.set(j,1);
	    j--;
	    break;
	  case 'G':
	    b.set(z,1);
	    z++;
	    b.set(z,0);
	    z++;
	    d.set(j,1);
	    j--;
	    d.set(j,0);
	    j--;
	    break;
	  case 'T':
	    b.set(z,1);
	    z++;
	    b.set(z,1);
	    z++;
	    d.set(j,0);
	    j--;
	    d.set(j,0);
	    j--;
	    break;
	  default:
	    findN=true;
	    c=window;
	  
	}
	
      }

      if (!findN)
      {
	z=0;	
	if (d<b)
	{
	  if (!(sh.search(d))) 
	  {
	    countko++;
	    vet[k]=false;
	  }
	  else
	  {
	    vet[k]=true;
	    countok++;
	    
	  }
	}
	else
	{
	  if (!(sh.search(b))) 
	  {
	    vet[k]=false;
	    countko++;
	    
	  }
	  else
	  {
	    vet[k]=true;
	    countok++;
	    
	  }
	  
	}
	
      }
    
  }

}




#if SIGNATURE


#if RAF
/**************************************************************/
/* NAME :  Class IOF search() */
/* DESCRIPTION : computing window signature x leukemia*/
/**************************************************************/
  void IOF::searchEachXSignature(const char buffer[],const int& num){

  #pragma omp critical
  {//critical section for hash update
  reads++;
  }
  class mybitset b,d;
  bool first=true;
  int size_tmp_string=0;
  char tmp_string[MAXSIZE];
  tmp_string[0]='\0';
  bool findN=false;
  for(int k=0; k<=num-(window)-1; k++)
  {
    int z=0;
    int j=window*2-1;
    for (int c = 0; c<window; c++) 
      {
	switch (buffer[c+k])
	{
	  case 'A':
	    b.set(z,0);
	    z++;
	    b.set(z,0);
	    z++;
	    d.set(j,1);
	    j--;
	    d.set(j,1);
	    j--;
	    break;
	  case 'C':
	    b.set(z,0);
	    z++;
	    b.set(z,1);
	    z++;
	    d.set(j,0);
	    j--;
	    d.set(j,1);
	    j--;
	    break;
	  case 'G':
	    b.set(z,1);
	    z++;
	    b.set(z,0);
	    z++;
	    d.set(j,1);
	    j--;
	    d.set(j,0);
	    j--;
	    break;
	  case 'T':
	    b.set(z,1);
	    z++;
	    b.set(z,1);
	    z++;
	    d.set(j,0);
	    j--;
	    d.set(j,0);
	    j--;
	    break;
	  default:
	    findN=true;
	    c=window;
	  
	}
	
      }

      if (!findN)
      {
	z=0;	
	if (d<b)
	{
	  if (!(sh.search(d))) 
	  {
	  first=true;
	  if (size_tmp_string>window+MINSIGNATURE){
	    insert_string(tmp_string,size_tmp_string,buffer);
	    tmp_string[0]='\0';
	    size_tmp_string=0;
	    }
	  }
	  else
	  {
	   if (first){
	     first=false;
	     strncpy(tmp_string,buffer+k,window);
	     tmp_string[window]='\0';
	     size_tmp_string=window;
	   }
	   else
	   {
	    tmp_string[size_tmp_string]= buffer[window+k-1];
	    ++size_tmp_string;
	    tmp_string[size_tmp_string]='\0'; 
	   }
	  }
	}
	else
	{
	  if (!(sh.search(b))) 
	  {
	  first=true;
	  if (size_tmp_string>window+MINSIGNATURE){
	    insert_string(tmp_string,size_tmp_string,buffer);
	    tmp_string[0]='\0';
	    size_tmp_string=0;
	    }
	  }
	  
	  else
	  {
	   if (first){
	     first=false;
	     strncpy(tmp_string,buffer+k,window);
	     tmp_string[window]='\0';
	     size_tmp_string=window;
	   }
	   else
	   {
	    tmp_string[size_tmp_string]= buffer[window+k-1];
	    ++size_tmp_string;
	    tmp_string[size_tmp_string]='\0';
	   }
	    
	  }
	  
	}
	
      }
  }
}
#else

/**************************************************************/
/* NAME :  Class IOF search() */
/* DESCRIPTION : computing window signature x leukemia*/
/**************************************************************/
  void IOF::searchEachXSignature(const char buffer[],const int& num){

  #pragma omp critical
  {//critical section for hash update
  reads++;
  }
  class mybitset b,d;
  int size=0;
  char tmp_string[MAXSIZE];
  memset(tmp_string,'\0',num);
  tmp_string[0]='\0';
  bool findN=false,fk=true;
  int last=0;
  //int last=0;
  //cout<<buffer<<endl;
  //cout<<"*************"<<endl;
  for(int k=0; k<=num-(window)-1; k++)
  {
    int z=0;
    int j=window*2-1;
    for (int c = 0; c<window; c++) 
      {
	switch (buffer[c+k])
	{
	  case 'A':
	    b.set(z,0);
	    z++;
	    b.set(z,0);
	    z++;
	    d.set(j,1);
	    j--;
	    d.set(j,1);
	    j--;
	    break;
	  case 'C':
	    b.set(z,0);
	    z++;
	    b.set(z,1);
	    z++;
	    d.set(j,0);
	    j--;
	    d.set(j,1);
	    j--;
	    break;
	  case 'G':
	    b.set(z,1);
	    z++;
	    b.set(z,0);
	    z++;
	    d.set(j,1);
	    j--;
	    d.set(j,0);
	    j--;
	    break;
	  case 'T':
	    b.set(z,1);
	    z++;
	    b.set(z,1);
	    z++;
	    d.set(j,0);
	    j--;
	    d.set(j,0);
	    j--;
	    break;
	  default:
	    findN=true;
	    c=window;
	  
	}
	
      }

      if (!findN)
      {
	z=0;	
	if (d<b){
	  if (!(sh.search(d))){
	  //tmp_string[k+window-1]='N';
	  fk=true;
	  }
	  else{
            int newk=k;
	    int newwindow=window;
	    if (fk){
	      if (k<last){
		newk=last;
	        newwindow=window+(k-last);		
		}
	      strncpy(tmp_string+size,buffer+newk,newwindow);
	      //cout<<tmp_string<<""<<endl;
	      size+=newwindow;
	      //cout<<newk<<" "<<newwindow<<" "<<size<<endl;
	      fk=false;
	    }
	    else{ 
	    tmp_string[size]=buffer[k+window-1];
	    //cout<<tmp_string<<" "<<buffer[k+window]<<" "<<size<<endl;
	    size++;
	    //cout<<k<<" 1 "<<size<<endl;
	    }
	   last=newk+newwindow; //qui
	  } 
	}
	else{
	  if (!(sh.search(b))){
	    //tmp_string[k+window-1]='N';
	    fk=true;
	  }
	  else{
	  int newk=k;
          int newwindow=window; 
	  if (fk){
	     if (k<last){
		newk=last;
		newwindow=window+(k-last);
	     }

	     strncpy(tmp_string+size,buffer+newk,newwindow);
	     //cout<<tmp_string<<endl;
	     size+=newwindow;
	     //cout<<newk<<" "<<newwindow<<" "<<size<<endl;
	     fk=false;
	   }
	   else{  
	    tmp_string[size]=buffer[k+window-1];
	     //cout<<tmp_string<<" "<<buffer[k+window]<<" "<<size<<endl;
	     size++;
	     //cout<<k<<" 1 "<<size<<endl;
	   }
	   last=newk+newwindow;
	  }
	}	
      }
    
  }
 //cout<<buffer<<endl;
 //cout<<tmp_string<<endl;
 //cout<<"*************"<<endl; 
 if (size>num){
   cout<<"Error "<<size<<endl;
   exit(1);
 }
 if (size>MINSIGNATURE){
  tmp_string[size]='\0'; 
  insert_string(tmp_string,size,buffer);
 } 
}



#endif




void IOF::insert_string(char tmp_string[],int size_tmp_string,const char buffer[]){
char r_tmp_string[MAXSIZE];
  if (size_tmp_string!=0){	    
    int i=0;
    for (int j=size_tmp_string-1;j>=0;--j){
      switch(tmp_string[j]){
	case 'A':
	 r_tmp_string[i]='T';
	break;
	case 'T':
	  r_tmp_string[i]='A';
	break;
	case 'C':
	  r_tmp_string[i]='G';
	break;     
	case 'G':
	  r_tmp_string[i]='C';
	break;
	case 'N':
	  r_tmp_string[i]='N';
	break;
      }
    ++i; 
    }
    r_tmp_string[size_tmp_string]='\0';

    for (i=0;(i<size_tmp_string)&&(tmp_string[i]-r_tmp_string[i]==0);++i);
#pragma omp critical
    {//critical section for hash update
    if ((i==size_tmp_string)||(tmp_string[i]-r_tmp_string[i]<0)){
    if (read_signature.find(tmp_string)!=read_signature.end())
    ++(read_signature[string(tmp_string)].frequency);
    else
    {
      Signature B(r_tmp_string,buffer,1);
      read_signature[string(tmp_string)]=B;
    }
    
    }  
    else
    {
    if (read_signature.find(r_tmp_string)!=read_signature.end())
       ++(read_signature[string(r_tmp_string)].frequency);//reverse compl
    else{
      Signature B(tmp_string,buffer,1);
      read_signature[string(r_tmp_string)]=B;;
      }
    }
    }//critical section for hash update
  }
}

#if RAF
/**************************************************************/
/* NAME :  Class IOF SignaturePrint() */
/* DESCRIPTION : printing the signature with their frequencies */
/**************************************************************/

void IOF::SignaturePrint(ostream& out){
 #pragma omp critical
  { 

//remove low frequency signature
  cout<<"Size before removing low frequency signatures: "<<read_signature.size()<<endl;
  auto it1=read_signature.begin();  
  while (it1!=read_signature.end()){
  if (it1->second.frequency<TFREQ){
  auto it2=it1;
  ++it1;
  read_signature.erase(it2);  
  }
  else
    it1++;
  }
  cout<<"Size after removing low frequency signatures: "<<read_signature.size()<<endl;
//remove low frequency signature    
    
//merge subsignature
  bool deleteit1=false;
  int i=0;
  it1=read_signature.begin();
  while (it1!=read_signature.end()){
    auto it2=it1;
    it2++;
    string r_T; 
    int size=it1->second.read.size();
    
    for(int i=size-1;i>=0;--i){
      switch(it1->second.read[i]){
	case 'A':
	 r_T.push_back('T');
	break;
	case 'T':
	  r_T.push_back('A');
	break;
	case 'C':
	  r_T.push_back('G');
	break;     
	case 'G':
	  r_T.push_back('C');
	break;
	case 'N':
	  r_T.push_back('N');
	break;
      }
    }
    while (it2!=read_signature.end()&&!(deleteit1)){ 
       if (it1->first.size()>it2->first.size()){
	if (((it1->first.find(it2->first)!=std::string::npos)||(it1->first.find(it2->second.r_signature)!=std::string::npos))){//&& (LevenshteinDistance(it1->second.read,it2->second.read)||(LevenshteinDistance(r_T,it2->second.read)))){
	  it1->second.frequency += it2->second.frequency;
	  auto it3=it2;
	  ++it2;
	  read_signature.erase(it3);
	}
	else{
	  ++it2;
	}
       }
      else{
	  if (((it2->first.find(it1->first)!=std::string::npos)||(it2->first.find(it1->second.r_signature)!=std::string::npos))){//&& (LevenshteinDistance(it1->second.read,it2->second.read)||LevenshteinDistance(r_T,it2->second.read))){
	  deleteit1=true;
	  }
	  else{
	++it2;
	}
      }
    }
    if (deleteit1){
      it2->second.frequency += it1->second.frequency;
      auto it3=it1;
      ++it1;
	read_signature.erase(it3);
	deleteit1=false;
    }
    else      
      ++it1;
    i++;
    /*
     if ((Edistance(it1->second.read,it2->second.read,r_T))){
	it1->second.frequency += it2->second.frequency;
	auto it3=it2;
	++it2;
	read_signature.erase(it3);
      }
      else
	++it2;
    }
    ++it1;
    i++;
    */
    if (i%1000==1)
      cout<<"Size:"<<read_signature.size()<<" Iter:"<<i<<endl;
	
  }
//merge subsignature 
  for (auto it=read_signature.begin(); it!=read_signature.end(); ++it){
    if (it->second.frequency>TFREQ)
      out<<it->first<<" "<<it->second.frequency<<" "<<it->second.read<<endl;
  }
  }
}

#else


void IOF::SignaturePrint(ostream& out){
 #pragma omp critical
  { 

//remove low frequency signature
  cout<<"Size before removing low frequency signatures: "<<read_signature.size()<<endl;
  auto it1=read_signature.begin();  
  while (it1!=read_signature.end()){
#if  Nnormalized 
//To cope with N    
  int Nnum=0; 
  int size=it1->second.read.size()-1;
  for(int i=size;i>=0;--i){
      switch(it1->second.read[i]){
	case 'N':
	 Nnum++;
	break;
      }
    }
  if (Nnum>1)  
    it1->second.frequencyNorm=it1->second.frequency/Nnum;
//To cope with N   
#endif  
  if (it1->second.frequency<TFREQ){
  auto it2=it1;
  ++it1;
  read_signature.erase(it2);  
  }
  else
    it1++;
  }
  cout<<"Size after removing low frequency signatures: "<<read_signature.size()<<endl;
//remove low frequency signature    
    
//merge subsignature
  bool deleteit1=false;
  int i=0;
  it1=read_signature.begin();
  while (it1!=read_signature.end()){
    auto it2=it1;
    it2++;
    string r_T; 
    int size=it1->second.read.size();
    for(int i=size-1;i>=0;--i){
      switch(it1->second.read[i]){
	case 'A':
	 r_T.push_back('T');
	break;
	case 'T':
	  r_T.push_back('A');
	break;
	case 'C':
	  r_T.push_back('G');
	break;     
	case 'G':
	  r_T.push_back('C');
	break;
	case 'N':
	  r_T.push_back('N');
	break;
      }
    }
    while (it2!=read_signature.end()&&!(deleteit1)){
      //if (/*(AminoDistance(it1->first,it2->first,it2->second.r_signature))&&*/((LevenshteinDistance(it1->first,it2->first))||(LevenshteinDistance(it1->first,it2->second.r_signature)))){//||(LevenshteinDistance(it1->second.read,it2->second.read)||(LevenshteinDistance(r_T,it2->second.read)))){
      if ((SmithWaterman(it1->first,it2->first)||SmithWaterman(it1->first,it2->second.r_signature))&&(SmithWaterman(it1->second.read,it2->second.read)||SmithWaterman(r_T,it2->second.read))){	
#if  Nnormalized  
//To cope with N   	
	if ( it1->second.frequencyNorm >= it2->second.frequencyNorm){	  
	  it1->second.frequencyNorm += it2->second.frequencyNorm;
//To cope with N   
#else
	  if ( it1->second.frequency >= it2->second.frequency){
#endif	  
	  it1->second.frequency += it2->second.frequency;
	  auto it3=it2;
	  ++it2;
	  read_signature.erase(it3);
	}
	else{
	  deleteit1=true;
	  
	}
	
      }
      else
	++it2;
    }
    if (deleteit1){
      it2->second.frequency += it1->second.frequency;

#if  Nnormalized        
//To cope with N         
      it2->second.frequencyNorm += it1->second.frequencyNorm;
//To cope with N 
#endif
      
      auto it3=it1;
      ++it1;
	read_signature.erase(it3);
	deleteit1=false;
    }
    else      
      ++it1;
    i++;
    if (i%10==1)
      cout<<"Size:"<<read_signature.size()<<" Iter:"<<i<<endl;
	
  }
//merge subsignature 
  for (auto it=read_signature.begin(); it!=read_signature.end(); ++it){
    if (it->second.frequency>TFREQ)
      out<<it->first<<" "<<it->second.frequency<<" "<<it->second.read<<endl;
  }
  }
}



#endif

bool IOF::AminoDistance(const string& T,const string& S,const string& r_T){
 
  int size=min3(S.size(),T.size(),r_T.size());
  int amini_freq[3][5] {0};
  for(int i=0;i<size;++i){
       switch(S[i]){
	case 'A':
	 ++amini_freq[0][0];
	break;
	case 'T':
	 ++amini_freq[0][1];
	break;
	case 'C':
	 ++amini_freq[0][2];
	break;     
	case 'G':
	 ++amini_freq[0][3];
	break;
	case 'N':
	 ++amini_freq[0][4];
	break;
      }
      switch(T[i]){
	case 'A':
	 ++amini_freq[1][0];
	break;
	case 'T':
	 ++amini_freq[1][1];
	break;
	case 'C':
	 ++amini_freq[1][2];
	break;     
	case 'G':
	 ++amini_freq[1][3];
	break;
	case 'N':
	 ++amini_freq[1][4];
	break;
      }
      switch(r_T[i]){
	case 'A':
	 ++amini_freq[2][0];
	break;
	case 'T':
	 ++amini_freq[2][1];
	break;
	case 'C':
	 ++amini_freq[2][2];
	break;     
	case 'G':
	 ++amini_freq[2][3];
	break;
	case 'N':
	 ++amini_freq[2][4];
	break;
      }
  }
  int max_1=0,max_2=0;
  for(int i=0;i<5;++i){
    if ( amini_freq[0][i]>amini_freq[1][i])
      max_1+=amini_freq[0][i]-amini_freq[1][i];
    if ( amini_freq[0][i]>amini_freq[2][i])
      max_2+=amini_freq[0][i]-amini_freq[2][i];
  }
  return (max_1<THEDIT||max_2<THEDIT)?true:false;
}    
   

bool IOF::SmithWaterman(const string& S,const string& T){

    // Declares a default Aligner
  StripedSmithWaterman::Aligner aligner;
  // Declares a default filter
  StripedSmithWaterman::Filter filter;
  // Declares an alignment that stores the result
  StripedSmithWaterman::Alignment alignment;
  // Aligns the query to the ref
  int size=0;
  if (S.size()>T.size()){
    size=T.size();
    if (size<S.size()*0.7) return false;
    aligner.Align(T.c_str(), S.c_str(), S.size(), filter, &alignment);
    
  }
  else{
     size=S.size();
     if (size<T.size()*0.7) return false;
     
     aligner.Align(S.c_str(), T.c_str(), T.size(), filter, &alignment);
     
  }
  
  return (alignment.sw_score>=(size*THEDIT)?true:false);
}




bool IOF::LevenshteinDistance(const string& s1,const string& s2){
  const size_t m(s1.size());
  const size_t n(s2.size());

  if (( m==0 ) ||( n==0 )){
      cerr << "\n*****Error size string for Edit Distance ("<<m<<","<<n<<")*****" << endl;
      cerr<<"*"<<s1<<"*"<<endl;
      cerr<<"*"<<s2<<"*"<<endl;
      exit(EXIT_FAILURE);
  }
 
  size_t costs[MAXSIZE];

 
  for( size_t k=0; k<=n; k++ ) costs[k] = k;
 
  size_t i = 0;
  for ( std::string::const_iterator it1 = s1.begin(); it1 != s1.end(); ++it1, ++i )
  {
    costs[0] = i+1;
    size_t corner = i;
 
    size_t j = 0;
    for ( std::string::const_iterator it2 = s2.begin(); it2 != s2.end(); ++it2, ++j )
    {
      size_t upper = costs[j+1];
      if( *it1 == *it2 )
      {
		  costs[j+1] = corner;
	  }
      else
	  {
		size_t t(upper<corner?upper:corner);
        costs[j+1] = (costs[j]<t?costs[j]:t)+1;
	  }
 
      corner = upper;
    }
  }
  int diff=abs((int)m-(int)n);
  if (diff > MINSIGNATURE)
    diff=0;
  return (costs[n]-diff)<THEDIT?true:false;
}
 


#endif


/**************************************************************/
/* NAME :  Class IOF FreqPrint() */
/* DESCRIPTION : printing the frequencies of the k-mers in the hash table*/
/**************************************************************/
void IOF::FreqPrint(){ 
  sh.FreqPrint();
  
}

/**************************************************************/
/* NAME :  Class IOF FreqKmerPrint() */
/* DESCRIPTION : printing the k-mer and their frequencies in an output file*/
/**************************************************************/
void IOF::FreqKmerPrint(){

    ofstream outfreq(string(OFname+".freq").c_str(),ofstream::out);
    if (!outfreq)
    {
      cerr << "\n*****Error opening ouput file "<< OFname.c_str()<<".freq *****" << endl;
      exit(EXIT_FAILURE);
    }
  sh.FreqKmerPrint(outfreq,window); 
}


/**************************************************************/
/* NAME :  Class IOF  ReadXKerm() */
/* DESCRIPTION :  it save in the output file the read kmer with their quality */
/**************************************************************/  
void IOF::ReadXKerm(int l, int u, PairedEnd p){
  
  using namespace general;
  clock_t startGlobal,endGlobal;
  char buffer[MAXSIZE],name[MAXSIZE],quality[MAXSIZE];
  
  ifstream in;
  ofstream out;
  reads=0;
  memset(count,0, files1*sizeof(unsigned long long));

  for (int i=l;i<=u;i++)
  {
    startGlobal=clock();
    ostringstream of;
    if (p == p1)
     of<<ReadFname1<<i<<"."<<fileExt1;
    else
     of<<ReadFname2<<i<<"."<<fileExt2;  
    in.open(of.str().c_str(),ifstream::in);
    if(!in) 
    {
      cerr << "\n*****Error opening input file "<<of.str().c_str() <<" *****" << endl;
      exit(EXIT_FAILURE);
      
    }
    of.str("");
    if (p == p1)
      of<<ReadFname1<<"Output"<<i<<"."<<fileExt1;
     else
     of<<ReadFname2<<"Output"<<i<<"."<<fileExt2;
     
    out.open(of.str().c_str(),ofstream::out);
  
    if(!out) 
    {
      cerr << "\n*****Error opening output file "<<of.str().c_str() <<" *****" << endl;
      exit(EXIT_FAILURE);
    }
    int hh=0;
    //For each read in the pool file
    while (!in.eof())
    {
      name[0]='\0';
      //Name
      in.getline(name,MAXSIZE);
      //Sequence
      buffer[0]='\0';
      in.getline(buffer,MAXSIZE);
      int num=in.gcount();
      if (buffer[num-1]!='\0')
      {
	buffer[num]='\0';
	num++;
      }
      //separator
      in.getline(quality,MAXSIZE);
      //Quality
      quality[0]='\0';
      in.getline(quality,MAXSIZE);
      if ((num>=window+1)&&((buffer[0]=='A')||(buffer[0]=='G')||(buffer[0]=='C')||(buffer[0]=='T')||(buffer[0]=='N')))
      {
	hh++;
	//insert window in the hash table
	writeKmer(name,buffer,quality,num,i,p,out);
	if(reads%10000000==0)
	{
	//#pragma omp critical
        //{
	  cout<<"Reads: "<<reads<<" windows: "<<count[i]<<" rejected  windows: "<<disc<<endl;
	//}  
	}
      }
      else
      {
	if ((buffer[0]!='\0'))
	{
        cerr<<"Error input format file "<<of.str().c_str()<<endl;
	exit(0);
	}
      }
      
    }
    in.close();
    out.close();
    endGlobal=clock();
    //#pragma omp critical
   //{
    cout<<"\n\n\tReads: "<<reads<<" windows: "<<count[i]<<"\n\tTime to read pool "<<i<<":"<<((double)(endGlobal-startGlobal))/CLOCKS_PER_SEC<<"s."<<endl;
    //}
  }

 
}





/**************************************************************/
/* NAME :  Class IOF write() */
/* DESCRIPTION : encoding and inserting windows */
/**************************************************************/				
void IOF::writeKmer(const char name[],const char buffer[],const char quality[],const int num,const int i,const PairedEnd p,ofstream &out){
  
// #pragma omp critical
//      { 
	reads++;
//	}
  count[i]+=(int)(num-1)/window;
  for(int k=0; k<=num-(window); k=k+window)
  {
    out<<name<<"SUB-"<<k<<"/";
    if (p == p1)
      out<<"1"<<endl;
    else
      out<<"2"<<endl;  
    out<<std::string(buffer+k,window)<<endl<<"+\n"<<std::string(quality+k,window)<<endl;
  }  
  
}

}

