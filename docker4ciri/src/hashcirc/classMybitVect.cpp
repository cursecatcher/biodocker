/***************************************************************************
 *   Copyright (C) 2013 by Marco Beccuti  				   *
 *   beccuti@di.unito.it   						   *
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



using namespace MYBIT;

/**************************************************************/
/* NAME :  Class mybitset*/
/* DESCRIPTION : Empty Constructor*/
/**************************************************************/  
mybitset::mybitset(){
  freqK=1;
  Inext = DEFAULTP;
}

/**************************************************************/
/* NAME :  Class mybitset invert()*/
/* DESCRIPTION : invert operator */
/**************************************************************/
inline void  mybitset::invert(){
  
  bool inv[DIM];
  int i=0;
  for (int jj=ksize-1;jj>=0;jj=jj-2)
  {
    if (cc[jj]==cc[jj-1])
      if(cc[jj]==1)
      {
	inv[i]=inv[i+1]=0;
      }
      else
      {
      	inv[i]=inv[i+1]=1;
      }
      else
      {
	inv[i]=cc[jj];
        inv[i+1]=cc[jj-1];
      }
    i=i+2;
  }
  for (int i=0;i<ksize; i++)
    cc[i]=inv[i];
  
}

/**************************************************************/
/* NAME :  Class mybitset compare()*/
/* DESCRIPTION : comparing two k-mers */
/**************************************************************/	
inline int  mybitset::compare(const mybitset& p2){
  
	int i=0;
 	while (((this->cc[i]==p2.cc[i]))&&(i<=ksize))
 			{
 			i++;		
 			}
 	if (i>ksize)
 			{
 			return 0;
 			}
 	if (this->cc[i]==1)
 			return 1;
 		else		
 			return 2;
	}

	
/**************************************************************/
/* NAME :  Class mybitset bit2buffer()*/
/* DESCRIPTION : returning k-mers */
/**************************************************************/	
	
inline void mybitset:: bit2buffer(char buffer[]){
  
	int zz=0;
	zz=0;
	for (int ii=0;ii<ksize;ii++)
		{
		if ((cc[ii]==0)&& (cc[ii+1]==0))
				buffer[zz]='A';
		else
			if ((cc[ii]==0)&&(cc[ii+1]==1))
				buffer[zz]='C';
			else
				if ((cc[ii]==1)&&(cc[ii+1]==0))
					buffer[zz]='G';
				else
					buffer[zz]='T';
			ii++;
			zz++;
		}
	buffer[zz]='\0';
	};

/**************************************************************/
/* NAME :  Class mybitset store_compact()*/
/* DESCRIPTION : encoding an integer into a char */
/**************************************************************/	
inline void mybitset::store_compact(unsigned int nval, char& c){
	register unsigned char cc;
	cc = (unsigned char)(0xFF & nval); 
	c=cc;
	};	
	
	
/**************************************************************/
/* NAME :  Class mybitset load_compact()*/
/* DESCRIPTION : un-encoding an integer  */
/**************************************************************/
void mybitset::load_compact(char& c,unsigned int& nval)const{
	register unsigned char cc0 = c;
	register unsigned int uu = (unsigned int)(cc0 & 0xFF);
	nval=uu;
	};

/**************************************************************/
/* NAME :  Class mybitset trans()*/
/* DESCRIPTION : encoding a bit-vector into a long   */
/**************************************************************/
unsigned long long int mybitset::trans(){
  
	register unsigned long long i=0;
	
	if (cc[(ksize)-1]==1)
		i++;
	for (int ii=(ksize)-2;ii>=0;ii--)
		{
		i=i<<1;
		if (cc[ii]==1)
			i++;
		}
	return i;
	};	
	
	
/**************************************************************/
/* NAME :  function operator<<*/
/* DESCRIPTION : operator<< */
/**************************************************************/
ostream& operator<<(ostream& out, const mybitset& p){
  
	int zz=0;
	string buffer="";
	zz=0;
	int ksize=mybitset::getKmer();
	for (int ii=0;ii<ksize;ii++)
		{
		if ((p.cc[ii]==0)&& (p.cc[ii+1]==0))
				buffer+='A';
		else
			if ((p.cc[ii]==0)&&(p.cc[ii+1]==1))
				buffer+='C';
			else
				if ((p.cc[ii]==1)&&(p.cc[ii+1]==0))
					buffer+='G';
				else
					buffer+='T';
			ii++;
			zz++;
		}
	
	out<<"\t"<<buffer<<"%";
	out<<p.freqK<<endl;	
#if DEBUG 
	if (p.Inext==DEFAULTP)
		out<<"\nNext: NULL"<<endl;
	else
		out<<"\nNext: "<<p.Inext<<endl;
	
#endif
	return out;
	};

	
/**************************************************************/
/* NAME :  function operator==*/
/* DESCRIPTION : operator==  */
/**************************************************************/
bool operator==(mybitset& p1,mybitset& p2){
        int i=0;
	int ksize=mybitset::getKmer();
  	while ((i<ksize)&&((p1.cc[i]==p2.cc[i])))
	{
	  i++;
	}
	if (i==ksize)
	  return true;
	else
	  return false;
	};

/**************************************************************/
/* NAME :  function operator<*/
/* DESCRIPTION : operator<  */
/**************************************************************/	

bool operator<(const mybitset& p1, const mybitset& p2){
  
  	int i=0;
	int ksize=mybitset::getKmer();
  	while (((p1.cc[i]==p2.cc[i]))&&(i<=ksize))
		{
		i++;
		}
	if (i>ksize)
		{
		return false;
		}
	if (p1.cc[i]==1)
		return false;
	else		
		return true;
	};


		