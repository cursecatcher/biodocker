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

#ifndef __BIT_H__
        #define __BIT_H__
        #include <bitset>

#endif
#ifndef __MTH_H__
        #define __MTH_H__
        #include <math.h>
#endif

#ifndef __SET_H__
        #define __SET_H__
        #include <list>
#endif

#ifndef __CON_H__
        #define __CON_H__
        #include "const.h"
#endif

#ifndef __STDL__
        #define __STDL__
        #include <stdlib.h>
#endif

#ifndef __STR_H__
        #define __STR_H__
        #include <string.h>

#endif

#ifndef __FSTREAM__
        #define __FSTREAM__
        #include <fstream>
#endif

#ifndef __IOS_H__
        #define __IOS_H__
        #include <iostream>
#endif 

 
namespace MYBIT{

using namespace std;


class mybitset{

  
public:
        //!Sequence  in binary format.
        bitset <DIM*2> cc;

        //!K-mer frequency
        unsigned int freqK;
	
    
        //!Pointer to the next element in the collision list
        unsigned int Inext; 

        //!Empty Constructor.
        mybitset();

	//!Copy operator
	mybitset& operator=(const mybitset& b) {cc=b.cc; freqK=b.freqK; Inext=b.Inext; return *this;};

        //!It sets to val  the ith position  of the sequence's bit vector. 
        inline void set(const int& i,int val){ cc[i]=val; };


        //!It returns the value of ith position of the sequence's bit vector.
        inline int get(const int& i){ return cc[i]; };

        //!It increments the kmer frequency  
        inline void incFreq(){ freqK++;};

        //!It return  the kmer frequency
        inline int getFreq(){ return freqK; };

        //!It inverts a k-mer encoded in the sequence's bit-vector.
        inline void invert();

        //!It compares  two k-mers  p1 and p2. Returns: a)0 when p1=p2; b)1 when p1>p2; c)2 when p1<p2. 
        inline int compare(const mybitset& p2);


        //!It returns the k-mer encoded in the bit-vector  (buffer). 
        inline void bit2buffer(char buffer[]);

        //!It compacts an unsigned integer into a char
        inline void store_compact(unsigned int nval, char& c);

        //!It un-compacts a char into an unsigned integer
        void load_compact(char& c,unsigned int& nval)const;

        //!It encodes a bitset in an unsigned long int
        unsigned long long int trans();

};

}
