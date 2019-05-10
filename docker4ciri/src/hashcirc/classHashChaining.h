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


#ifndef __CON_H__
	#define __CON_H__
	#include "const.h"
#endif



#ifndef __FSTREAM__
	#define __FSTREAM__
	#include <fstream>
#endif

#ifndef __IOS_H__
	#define __IOS_H__
	#include <iostream>
#endif 


#ifndef __STDL__
	#define __STDL__
	#include <stdlib.h>
#endif

#ifndef __VCT_H__
	#define __VCT_H__
	#include <vector>
#endif


#ifndef __LMT_H__
	#define __LMT_H__
	#include <limits.h>
#endif




namespace Cl_HASH
{
using namespace  MYBIT;
using namespace std;


template <class A_Type> 

class HASHCH{

private:
//!Hash table
unsigned int * HashVec;
//!Collisions' list
A_Type * ListCol;
//!Hash table size 
unsigned int SizeHash;
//!Whole Collisions' list size.
//unsigned int mod;
unsigned int SizeList;
//!Number of elements in the hash table
unsigned int NumEl;

public:
  	//!Empty Constructor.
	HASHCH();
	
	//!Constructor. It takes in input the Hash table size in terms of number of buckets and max number of elements.
	HASHCH(const unsigned int sizeHash, const unsigned int sizeList);
	
	//!It updates the hash size
	void UpdateHash(const unsigned int& sizeHash,const unsigned int& sizeList);
	 
	//!it doubles  the hash table size  
	void IncreaseHash();
	
	//!Hash function. It takes in input a integer value.
	inline unsigned int hash(const unsigned long long& value);

 	//!It inserts the input value in the hash table, using chaining policy.  
	inline void insert(A_Type& b);
	
	//!It resets the hash table from ``index'' position and sets NumEl to 0
	inline void resetHT(unsigned int index);
	
	//!It searches the input elements into the hash table. It returns read signature. 
	void search(HASHCH<mybitset>& sh,unsigned int& countOk,unsigned int& countKo);
	
  	//!It searches  a value in the hash table and if it is found then it returns true,  false otherwise.
	inline bool search(A_Type& b);
	
	//!It searches  a value in the hash table and  if it is found then it updates its frequency
	inline bool searchIncFreq(A_Type& b);
	
	//!It prints the collision list of HashTable[row] into a file
	void printList(int row,ofstream& out);

	//!It prints the linked list of HashTable[row]
	void printList(int row);

	//!It prints the size of the linked list of  HashTable[row]
	void printSize(int row);

	//!It stores the hash table  into  a file 
	void print(ofstream& out);

	//!It prints the hash table
	void print();
	
	//!It returns the hash table size
	unsigned int size();
	
	//!It prints some statistical information on the hash table
	void info();
	
	//!It takes in input the output file name and saves on it (in binary format) the hash table.
	void Write(string output);

	//!It takes in input the intput file name and creates the hash table
	bool Read(string output);
	
	//!It prints the frequencies of the k-mers in the hash table
	void FreqPrint(void);
	
	//!It prints in a file the k-mers with their frequencies
	void FreqKmerPrint(ofstream& fout,const int window);
	
	//!It removes the allocated memory
	~HASHCH(){ free(HashVec);  free(ListCol); };
};
  
}