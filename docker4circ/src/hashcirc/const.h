/***************************************************************************
 *   Copyright (C) 2013 by Marco Beccuti   *
 *   beccuti@di.unito.it   *
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



//!k-mer size.
#define DIM  32

//!Threshold frequencies of k-mers 
#define THFREQ 0


//!Hash function alpha coefficient.
#define A 48271
//!Hash function prime number
#define L 2147483647

//!Max read size
// #define READSIZE 2048
#define READSIZE 262144


//!Default value used for  null pointer
#define DEFAULTP  4294967295U

//!max buffer's size.
// #define MAXSIZE 2048
#define MAXSIZE 262144

//!Insertion cost for edit distance
#define INS 1
//!Deletion cost for edit distance
#define DEL 1
//!Substitution cost for edit distance
#define SUB 0.5

//!Threshold value for edit distance x reads (LEUKEMIA)
#define THEDIT 1.47
//!Minimum number of k-mers for signature (must be >)
#define MINSIGNATURE 50
//!Threshold value for frequencies (LEUKEMIA)
#define TFREQ 25

//!min function for edit distance
#define min3(a,b,c) (a < b ? (a < c ? a : c) : (b < c ? b : c))

//!minimum size to perform alignment using SW in IGH
#define MINSIZEALIGN 5