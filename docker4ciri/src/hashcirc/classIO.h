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

#ifndef __FSTREAM__
	#define __FSTREAM__
	#include <fstream>
#endif

#ifndef __CON_H__
	#define __CON_H__
	#include "const.h"
#endif



#include <omp.h>	// OMP multithreading


#ifndef __SSTREAM__
	#define __SSTREAM__
	#include <sstream>
#endif

#ifndef __IOS_H__
	#define __IOS_H__
	#include <iostream>
#endif 


#ifndef __STDL__
	#define __STDL__
	#include <stdlib.h>
#endif

#ifndef __STR_H__
	#define __STR_H__
	#include <string.h>

#endif


 
#ifndef __TIM_H__
	#define __TIM_H__
	#include <time.h>

#endif


#ifndef __CMB_H__
  #define __CMB_H__
  #include "classMybitset.h"
#endif

#ifndef HASH_H__
  #define HASH_H__
  #include "classHashChaining.h"
#endif

#ifndef __GNR_H__
	#define __GNR_H__
	#include "general.h"
#endif

#include "conf.h"

#if SIGNATURE
  #ifndef __MAP_H__
	#define __MAP_H__
	#include <map>
  #endif
  #ifndef __SSW_H__
      #define __SSW_H__
      #include "ssw_cpp.h"
      #include "ssw.h"
      #include <vector>
  #endif
#endif


#define DEBUGBAC 0
#define DEBUG 0


namespace IOFn
{
using namespace std;
using namespace MYBIT;
using namespace Cl_HASH;
using namespace general;

enum PairedEnd {p1,p2};

#if SIGNATURE
class Signature 
{
  public:
  //It encodes the reverse complement of the signature   
  string r_signature;
  //It stores the representative sample read
  string read;
  //It encodes the signature frequency 
  unsigned int frequency;
#if Nnormalized  
   //It encodes the signature frequency normalized by the number of N bases
  unsigned int frequencyNorm {0};
#endif  
  Signature():r_signature{""},read{""},frequency{0}{};
  Signature(string r_sig, string r, unsigned int freq):r_signature{r_sig},read{r},frequency{freq}{};
   Signature(char r_sig[], char r[], unsigned int freq):r_signature{string(r_sig)},read{string(r)},frequency{freq}{};
  Signature(const class Signature& B):r_signature{B.r_signature},read{B.read},frequency{B.frequency}{};
};
#endif

class IOF
{
  //!Input file1, this name is concatenated with an integer, so we assume the input file name [A-Ba-b0-9]+[0-1]+
  string ReadFname1;
  //!Input file1, this name is concatenated with an integer, so we assume the input file name [A-Ba-b0-9]+[0-1]+
  string ReadFname2;
  //!Input file1 extension.
  string  fileExt1;
  //!Input file1 extension.
  string  fileExt2;
  //!Output file
  string OFname;
  //!Window size
  int window;
  //!Read number
  unsigned int reads;
  //!Read counting
  unsigned long long *count;
  //!Reject windows number
  int disc;
  //!number of input files1
  int files1;
  //!number of input files1
  int files2;
  //! HASH TABLE for windows
  HASHCH<mybitset> sh;
  //! minimum found kmer numbers used to select the read
  unsigned int numker;
  
#if SIGNATURE
  //! SET containing the reads' signatures for leukemia
  map<string, class Signature> read_signature;
  //! It is used for edit distance.
  double  d[2][MAXSIZE];
#endif
  
public:
      //!Construtor. takes in input: 1)read file name; 2)read file name extention;  2)window size;  3)the number of input files.
      IOF(std::string ReadFname1,std::string fileExt1,std::string ReadFname2,std::string fileExt2,const std::string& OFname,const int& window,const int& file1,const int& file2,const unsigned int& numker);
      
      //!It takes in input the hash table size (number of buckets and max  number of elements) and allocates it.
      void AllocHash(const unsigned int& hash,const unsigned int&list);
      
      //!It resets the number of reads
      void ResetReads();
      
      //!It reads the input files containing reads and stores the their data into the Hash table 
      void Read(int l, int u);
      
      //!It reads the input files containing reads and generates for each read its k-mers saving in a output file.
      void ReadXKerm(int l, int u, PairedEnd p);
      
      //!It reads the input files (.count, .hash and .collist),  storing in binary format the hash tables 
      bool ReadB();
      
      //!It writes the hash tables  in three output files (.count, .hash and .collist).
      void WriteB();
	
      //!It searches in the hash table the windows of all reads, and derives  for each  read its signature. Its inputs are the interval  of pools on which is executed
      void Search(int l, int u);
	
      //!It prints the pool hash table
      void  print(){ sh.print();};
      
      //!It prints the frequencies of the k-mers in the hash table
      void FreqPrint();
      
      //It prints in the Output file hte k-mers and their frequencies in a file (namely OFname.freq)
      void FreqKmerPrint(void);

#if SIGNATURE     
      //!It prints the signatures with their frequencies
      void SignaturePrint(ostream& out);
#endif      
      //!Descrutor
      ~IOF(){
	free(count);
      }
      
private:
      //!It encodes a buffer[num] and the pool number  on two bitsets, and  it stores them into the hash table.
      inline void  insert(const char buffer[],const int& num,const int& pool);
      //!It encodes a buffer[num] in bitsets, and it searches  into the hash table. If the window is presented in the hash table then  It returns its signature
      void searchAll(class HASHCH<mybitset>& rd,const char buffer[],const int& num,ofstream& out,unsigned int& countok,unsigned int& countko);
       //!It encodes a buffer[num] in bitsets, and it searches  into the hash table. If the window is presented in the hash table then  It returns its signature
      void searchEach(bool vetKerm[],const char buffer[],const int& num,ofstream& out,unsigned int& countok,unsigned int& countko);
      //!It writes in the output file the k-mers of a read.
      void writeKmer(const char name[],const char buffer[],const char quality[],const int num,const int i,const PairedEnd p,ofstream &out);
#if SIGNATURE
      //!It computes and stores the read signature.
      //!Read signatures is the greatest sequence of continuous k-mers present in the Hash Table for the read.
      void searchEachXSignature(const char buffer[],const int& num);
      //!It inserts the read signature in the map.
      void insert_string(char tmp_string[],int size_tmp_string,const char buffer[]);
      //!It computes the  Levenshtein Distance and returns true if its values is lower than a threshold. 
      bool LevenshteinDistance(const string& S,const string& T);
       //!It computes the  mismatch Distance and returns true if its values is lower than a threshold. 
      bool AminoDistance(const string& S, const string& T, const string& r_T);
      //!It implements Smith-Waterman algorithm
      bool SmithWaterman(const string& S,const string& T);
#endif      
};

struct thread_data{
   unsigned int l;
   unsigned int u;
};

}
