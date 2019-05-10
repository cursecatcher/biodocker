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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef THR_H_
        #define THR_H_
        #include <pthread.h>
#endif

#ifndef CIO_H_
  #define CIO_H
  #include "classIO.h"
#endif

#ifndef TIM_H_
  #define TIM_H
  #include <sys/time.h>
#endif

#ifndef RES_H_
  #define RES_H
  #include <sys/resource.h>
#endif

#include <stdlib.h>

using namespace IOFn;

bool ContSeach=false;
unsigned int* NotFReads;
//int mybitset::ksize=0;

class IOF* p_ref; 

struct thread_data* td_array;

int main(int argc, char **argv) {
  clock_t startGlobal,endGlobal;
  time_t time1,time2;
  time1=time(NULL);
  startGlobal=clock();
  
  cout<<"\n\n =========================================================\n";
  cout<<"|                      HashChecker                       |\n";
  cout<<"|                     READs FILTER                       |\n";
  cout<<" =========================================================\n";
  cout<<"\n If you find any bug, send an email to beccuti@di.unito.it\n";

  if (argc!=13)
  {
    std::cerr<<"\n\nUSE: hashchecker <path/ReadFile1> <file_extension1> <input_file_number1> <path/ReadFile2>  <file_extension2>  <input_file_number2> <path/OutputFile> <k-mer_size>  <thread> <hash_size> <collision_list_size> <num_kmer>\n\n";
    exit(EXIT_FAILURE);
  }
  
 //Initialize input data 
 unsigned int qgram=atoi(argv[8]);
 //mybitset::setKmer(qgram*2);
 if (qgram>DIM)
 {
   std::cerr<<"\n\nThe windows size is greater than DIM. You have to update DIM in const.h\n\n";
   exit(EXIT_FAILURE);
 }
  

 unsigned int cores=atoi(argv[9]); 
 unsigned int files1=atoi(argv[3]);
 unsigned int files2=atoi(argv[6]);
 unsigned int hash=atoi(argv[10]);
 unsigned int list=atoi(argv[11]); 
 unsigned int numker=atoi(argv[12]); 
 //Initialize input data 
 
 //Initialize input/output files
 string ReadFname1(argv[1]);
 string ExtFile1(argv[2]);
  string ReadFname2(argv[4]);
 string ExtFile2(argv[5]);
 
 string IOread=std::string(argv[7]);
  //Initialize input/output files
 
 cout<<"\n\n___________________________________________________________\n\n\t\t\tInput parameters \n___________________________________________________________\n\n";
 
 cout<<"\tInput 1:";
  cout<<"\n\t\tfiles name: "<<ReadFname1;
  cout<<"\n\t\tInput files extension: "<<ExtFile1;
  cout<<"\n\t\t#Input files: "<<files1;
 cout<<"\n\n\tInput 2:";
  cout<<"\n\t\tfiles name: "<<ReadFname2;
  cout<<"\n\t\tfiles extension: "<<ExtFile2;
  cout<<"\n\t\t#Input files2: "<<files2;
 cout<<"\n\n\tOutput:";
  cout<<"\n\t\tfiles name: "<<IOread<<endl;
 cout<<"\n\tWindow's size: "<<qgram<<"\n";
 cout<<"\tThread numbers for mapping: "<<cores<<"\n";
 cout<<"\tHash table entries: "<<hash<<"\n";
 cout<<"\tList collision size: "<<list<<"\n";
 cout<<"\tNumber of found kmers to accept a read: "<<numker<<endl;
#if DEBUGOUTPUT
 cout<<"\tDEBUGOUTPUT is enable\n";
#endif
#if  DEBUGFREQ
 cout<<"\tDEBUGFREQ is enable\n";
#endif 
#if   SIMILARITY
 cout<<"\tSIMILARITY is enable\n";
#endif
#if  SIGNATURE
 cout<<"\tSIGNATURE is enable\n";
#endif
#if  PAIREDEND
 cout<<"\tPAIREDEND is enable\n";
#endif 
 cout<<"___________________________________________________________\n";

 cout<<"\n\nSTART FORMATTING...\n"<<endl;
 
 

 td_array = (thread_data*) malloc (cores*sizeof(thread_data));
 unsigned int last=0,dim=files1/cores;
 unsigned int tmpcore=cores;
 //!Check if the number of files is lower than the assigned cores.
 if (cores>files1)
  {
  cores=files1;
  dim=1;
  cout<<"Reset the number of cores used to "<<cores<<endl;
  }

   
  
 
 
 
 //!create variable for deriving the memory utilization  
 int who = RUSAGE_SELF;  
 struct rusage usage;
 //!create variable for deriving the memory utilization  
 
 class IOF ref(ReadFname1,ExtFile1,ReadFname2,ExtFile2,IOread,qgram,files1,files2,numker);
 
 //Init threads
 p_ref=&ref;
 
 startGlobal=clock();
 //!Initialize  hash tables 
 bool storedhash=ref.ReadB();
 if (!storedhash)
 {
  ref.AllocHash(hash,list);
    //!Assign the pools to each available core
  ref.ResetReads();  
  ref.Read(0,files1-1);
 }
 
 
 
 //compute memory utilization
 getrusage(who,&usage); 
 //compute memory utilization
 
 cout<<"\n\n\tTotal memory used: "<<usage.ru_maxrss<<"KB"<<endl;
 cout<<"\n\nEND FORMATTING\n"<<endl;

 
 endGlobal=clock();
 time2=time(NULL);
 cout<<"\n=========================== TIME ===========================\n\n";
  cout<<"\tReal Time to  create HASH TABLE:: "<<time2-time1<<"s."<<endl;
 cout<<"\tTime to create HASH TABLE: "<<((double)(endGlobal-startGlobal))/CLOCKS_PER_SEC<<"s."<<endl;
 cout<<"\n============================================================\n\n";
 
 
//!Save hash tables
 if (!storedhash)
 {
   cout<<"\n---------------------- SAVE HASH TABLE ---------------------\n\n";
   ref.WriteB();
   cout<<"\n------------------------------------------------------------\n\n";
  }

  
#if DEBUGFREQ == 1
  {  
  ref.FreqPrint();
  exit(0);
  }
#elif DEBUGFREQ == 2  
  {  
  ref.FreqKmerPrint();
  exit(0);
  }
#endif 
  

 time1=time(NULL);
 startGlobal=clock();
 cout<<"\n\nSTART MAPPING...\n"<<endl;

 cores=tmpcore;
 dim=files2/cores;
 last=0;
 //!Check if the number of files is lower than the assigned cores.
 if (cores>files2)
  {
   cores=files2;
   dim=1;
   cout<<"Reset the number of cores used to "<<files2<<endl;
  } 
 omp_set_num_threads(cores);

 //!Assign the pools to each available core
 for (unsigned int i=0;i<cores;i++)
 {
   td_array[i].l=last;
   td_array[i].u=last+dim-1;
   if ((td_array[i].u>files2-1)||(i==cores-1))
     td_array[i].u=files2-1;
   else
     last=last+dim;
   cout<<"\t Thread "<<i<<" works on the pools from "<<td_array[i].l<<" to "<<td_array[i].u<<endl;
 }
 cout<<endl;

 //!Start checking in parallel
  ref.ResetReads(); 
  
  NotFReads=(unsigned int*)malloc(files2*sizeof(unsigned int));
  
 #pragma omp parallel 
 {
   int id = omp_get_thread_num();
   p_ref->Search(td_array[id].l,td_array[id].u);
 }
 
 
 time2=time(NULL);
 //End threads
 cout<<"\n\nEND MAPPING\n"<<endl;
 endGlobal=clock();
 
 cout<<"\n=========================== TIME ===========================\n\n";
 cout<<"\tReal Time to find overlappings: "<<time2-time1<<"s."<<endl;
 cout<<"\tClock Time to find overlappings: "<<((double)(endGlobal-startGlobal))/CLOCKS_PER_SEC<<"s."<<endl;
 cout<<"\n============================================================\n\n";
 
#if !SIGNATURE
 cout<<"\n====================== FINAL STATISTIC =====================\n\n";
 for (unsigned int i=0;i<files2;i++)
  cout<<"File "<<i<<"\n\tTot. not found reads: "<<NotFReads[i]<<endl;
 cout<<"\n============================================================\n\n";
 //free memory
 free(NotFReads);
#endif  
 exit(EXIT_SUCCESS);
}
