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

using namespace  Cl_HASH;

 /**************************************************************/
/* NAME :  Class HASHCH*/
/* DESCRIPTION : Empty Constructor*/
/**************************************************************/ 
template <class A_Type>
inline HASHCH<A_Type>::HASHCH(){
  
  HashVec=NULL;
  SizeHash=0;
  SizeList=0;
  NumEl=0;
}

 /**************************************************************/
/* NAME :  Class HASHCH*/
/* DESCRIPTION : Constructor*/
/**************************************************************/ 
template <class A_Type>
inline HASHCH<A_Type>::HASHCH(const unsigned int sizeHash, const unsigned int sizeList){

  this->SizeHash=sizeHash;
  this->SizeList=sizeList;
  NumEl=0;
  HashVec=(unsigned int*)malloc(SizeHash*sizeof(unsigned int));
  for (unsigned int i=0;i<SizeHash;i++)
  {
    HashVec[i]=DEFAULTP; 
  }
  ListCol=(A_Type*)malloc(SizeList*sizeof(A_Type));
}


/**************************************************************/
/* NAME :  Class HASHCH UpdateHash() */
/* DESCRIPTION :  re-size hash table*/
/**************************************************************/ 
template <class A_Type>
void  HASHCH<A_Type>::UpdateHash(const unsigned int& sizeHash,const unsigned int& sizeList){
        this->SizeHash=sizeHash;
        this->SizeList=sizeList;
        HashVec=(unsigned int*)malloc(SizeHash*sizeof(unsigned int));
        for (unsigned int i=0;i<SizeHash;i++)
                {
                HashVec[i]=DEFAULTP;
                }
        //memset(HashVec,DEFAULTP, SizeHash*sizeof(unsigned long) );
        ListCol=(A_Type*)malloc(SizeList*sizeof(A_Type));
        };

/**************************************************************/
/* NAME :  Class HASHCH IncreaseHash() */
/* DESCRIPTION :  it doubles the hash table size */
/**************************************************************/ 	
template <class A_Type>
void  HASHCH<A_Type>::IncreaseHash(){
      SizeList*=2;
      A_Type * tmpListCol=(A_Type*)realloc(ListCol,(SizeList*sizeof(A_Type)));
      if (!tmpListCol){
	cerr << "\n*****Error the collision list cannot be increased *****" << endl;
	exit(EXIT_FAILURE);
      }
      else
	ListCol=tmpListCol;
      };

 /**************************************************************/
/* NAME :  Class HASHCH hash()*/
/* DESCRIPTION :  hash function implementation*/
/**************************************************************/ 
template <class A_Type>
inline unsigned int  HASHCH<A_Type>::hash(const unsigned long long& value){
        //return (H(value)) % SizeHash;
        register unsigned  long long u=value;
        register unsigned  long long v=value>>32;
        //return (((u*2+v)*A)%L)%SizeHash;
        return ((u*2+v)*A)%SizeHash;
        //return (value%L) % SizeHash;
        };


/**************************************************************/
/* NAME :  Class HASHCH insert()*/
/* DESCRIPTION :  insertion into the hash table*/
/**************************************************************/
template <class A_Type>
inline void  HASHCH<A_Type>::insert(A_Type& b){
        unsigned int probe = hash(b.trans());
//      #pragma omp critical
//      {
        unsigned int  tmpInext = HashVec[probe];
        HashVec[probe] = NumEl;
        ListCol[NumEl]=b;
        ListCol[NumEl].Inext = tmpInext;
        NumEl++;
//      }
        if (NumEl+1>SizeList)
                {
                cerr << "\n*****Error the number of elements ("<<NumEl <<") is bigger then list entries *****" << endl;
                cerr << "\n***** The size of the collision list will be increased *****" << endl;
		this->IncreaseHash();
                }
        };


/**************************************************************/
/* NAME :  Class HASHCH resetHT() */
/* DESCRIPTION :  resetting  hash table  */
/**************************************************************/ 
template <class A_Type>
inline void  HASHCH<A_Type>::resetHT(unsigned int index){
        NumEl=0;
        for (unsigned i=index;i<SizeHash;i++)
                HashVec[i]=DEFAULTP;
        }


/**************************************************************/
/* NAME :  Class HASHCH search() */
/* DESCRIPTION :  searching an element in the hash table  */
/**************************************************************/ 
template <class A_Type>
inline bool HASHCH<A_Type>::search(A_Type& b){
        unsigned int tmpInode = HashVec[hash(b.trans())];
        while(tmpInode !=DEFAULTP)
                {
                if(ListCol[tmpInode]==b) 
                        {
                        return ((ListCol[tmpInode].getFreq()> THFREQ)? true : false);
                        }
                tmpInode= ListCol[tmpInode].Inext;
                }
        return false;
        };

/**************************************************************/
/* NAME :  Class HASHCH searchSetPool() */
/* DESCRIPTION :  searching an element in the hash table  */
/**************************************************************/ 
template <class A_Type>
inline bool HASHCH<A_Type>::searchIncFreq(A_Type& b){

        unsigned int tmpInode = HashVec[hash(b.trans())];
        while(tmpInode !=DEFAULTP)
                {
                if(ListCol[tmpInode]==b) 
                        {
                        //#pragma omp critical
                        //{
                        ListCol[tmpInode].incFreq();
                        //}
                        return true;
                        }
                tmpInode = ListCol[tmpInode].Inext;
                }
        return false;
        };


/**************************************************************/
/* NAME :  Class HASHCH printList() */
/* DESCRIPTION : storing a collision list into a file */
/**************************************************************/ 
template <class A_Type>
void HASHCH<A_Type>::printList(int row,ofstream& out){
        unsigned int tmpInode = HashVec[row];
        while(tmpInode!=DEFAULTP)
                {
                out<<ListCol[tmpInode]<<endl;
                tmpInode = ListCol[tmpInode].Inext;
                }
        };


/**************************************************************/
/* NAME :  Class HASHCH printList() */
/* DESCRIPTION : printing a collision list  */
/**************************************************************/ 
template <class A_Type>
void HASHCH<A_Type>::printList(int row){
        unsigned int  tmpInode = HashVec[row];
        if (tmpInode!=DEFAULTP)
                cout<<"\nHASH["<<row<<"]:\n"; 
        while(tmpInode!=DEFAULTP)
                {
                cout<<ListCol[tmpInode]<<" ";
                tmpInode = ListCol[tmpInode].Inext;
                }
        };


/**************************************************************/
/* NAME :  Class HASHCH printSize() */
/* DESCRIPTION : printing  collision list size  */
/**************************************************************/ 
template <class A_Type>
void  HASHCH<A_Type>::printSize(int row){
        unsigned int tmpInode = HashVec[row];
        int i=0;
        if (tmpInode!=DEFAULTP)
                cout<<"\nHASH["<<row<<"]:\t"; 
        while(tmpInode!=DEFAULTP)
                {
                i++;
                tmpInode = ListCol[tmpInode].Inext;
                }
        if (i>0)
                cout<<i<<endl;
        };


/**************************************************************/
/* NAME :  Class HASHCH print() */
/* DESCRIPTION : storing the hash table into a file */
/**************************************************************/ 
template <class A_Type>
void HASHCH<A_Type>::print(ofstream& out){
        for(unsigned int i=0;i<SizeHash;i++)
                {
                printList(i,out);
                }
        };


/**************************************************************/
/* NAME :  Class HASHCH print() */
/* DESCRIPTION : printing the hash table */
/**************************************************************/ 
template <class A_Type>
void HASHCH<A_Type>::print(){
        for(unsigned int i=0;i<SizeHash;i++)
                {
//#if DEBUG1
                printList(i);
//#endif
//              printSize(i);
                }
        cout<<endl;
        };



/**************************************************************/
/* NAME :  Class HASHCH size() */
/* DESCRIPTION : hash table size */
/**************************************************************/ 
template <class A_Type>
unsigned int HASHCH<A_Type>::size(){
        return NumEl;
        };

/**************************************************************/
/* NAME :  Class HASHCH info() */
/* DESCRIPTION : hash table information */
/**************************************************************/ 
template <class A_Type>
void HASHCH<A_Type>::info(){
  
        unsigned int hashentry=0;
        for (unsigned int i=0;i<SizeHash;i++)
                {
                if (HashVec[i]!=DEFAULTP)
                        hashentry++;
                }
        cout<<"Hash table buckets: "<<hashentry<<endl;
        cout<<"Avergage size of collisions list: "<<(double)NumEl/hashentry<<endl;
        }

/**************************************************************/
/* NAME :  Class HASHCH Write() */
/* DESCRIPTION : hash table writer*/
/**************************************************************/ 
template <class A_Type>
void HASHCH<A_Type>::Write(string output){
  
  
        FILE * pFile;
        string tmp=output+".hash";
        pFile = fopen ( tmp.c_str(), "w" );
        if (pFile==NULL)
                {
                cerr<<"\n*****Error opening output file .hash *****" << endl;
                exit(EXIT_FAILURE);
                }
        //write hash size
        fwrite (&SizeHash, sizeof(SizeHash),1 , pFile );
        //write int vector
        fwrite (HashVec, sizeof(unsigned int), SizeHash, pFile );
        fclose (pFile);
        cout<<"Hash table saved in:\n\t"<<tmp<<endl;
        tmp=output+".collist";
        pFile = fopen ( tmp.c_str(), "w" );
        //write collisions list size
        fwrite (&SizeList, sizeof(SizeList), 1, pFile );
        //write number of elements
        fwrite (&NumEl, sizeof(NumEl), 1, pFile );
        //write myclassbit vector
        fwrite (ListCol, sizeof(A_Type), SizeList, pFile );
	/*
	for (unsigned int i=0;i<NumEl;i++)
	{
	 fwrite(ListCol[i].cc,sizeof(bool),(ListCol[0].ksize),pFile);
	} */ 
        fclose (pFile);
        cout<<"Collisions' list saved in:\n\t"<<tmp<<endl;
        };

/**************************************************************/
/* NAME :  Class HASHCH Read() */
/* DESCRIPTION : hash table reader*/
/**************************************************************/ 
template <class A_Type>
bool HASHCH<A_Type>::Read(string output){
  
        FILE * pFile;
        string tmp=output+".hash";
        pFile = fopen ( tmp.c_str(), "r" );
        if (pFile==NULL)
                {
                cerr<<"\n*****Error opening input file .hash *****" << endl;
                cout<<"\nHash Table binary encoding not found!"<<endl;
                return false;
                }
        //write hash size
        if (fread (&SizeHash, sizeof(SizeHash),1 , pFile )!=1)
                {
                cerr<<"Reading error"<<endl; 
                exit (EXIT_FAILURE);
                }
        HashVec=(unsigned int*)malloc(SizeHash*sizeof(unsigned int));
        //write int vector
        if (fread (HashVec, sizeof(unsigned int), SizeHash, pFile )!=SizeHash)
                {
                cerr<<"Reading error"<<endl; 
                exit (EXIT_FAILURE);
                }
        fclose (pFile);
        tmp=output+".collist";
        pFile = fopen ( tmp.c_str(), "r" );
        if (pFile==NULL)
                {
                cerr<<"\n*****Error opening input file .hash *****" << endl;
                cout<<"\nHash Table binary encoding not found!"<<endl;
                return false;
                }
        //write collisions list size
        if (fread (&SizeList, sizeof(SizeList), 1, pFile )!=1)
                {
                cerr<<"Reading error"<<endl; 
                exit (EXIT_FAILURE);
                }
        //write number of elements
        if (fread (&NumEl, sizeof(NumEl), 1, pFile )!=1)
                {
                cerr<<"Reading error"<<endl; 
                exit (EXIT_FAILURE);
                }
        //write myclassbit vector
        ListCol=(A_Type*)malloc(SizeList*sizeof(A_Type));
        if (fread (ListCol, sizeof(A_Type), SizeList, pFile )!=SizeList)
                {
                cerr<<"Reading error"<<endl; 
                exit (EXIT_FAILURE);
                }          
       /* for (unsigned int i=0;i<NumEl;i++)
	{
	 ListCol[i].cc=(bool*)malloc(sizeof(bool)* mybitset::ksize); 
	 if (fread(ListCol[i].cc,sizeof(bool),(ListCol[0].ksize),pFile)!=(unsigned int)(ListCol[0].ksize))
                {
                cerr<<"Reading error"<<endl; 
                exit (EXIT_FAILURE);
                }
	}*/          
//      for (int j=0;j<SizeList;j++)
//              cout<<"["<<j<<"]"<<ListCol[j]<<endl;
        fclose (pFile);
        cout<<"\n\tHash Table encoding found";
        return true;

}

/**************************************************************/
/* NAME :  Class HASHCH search() */
/* DESCRIPTION :  searching in the hash table*/
/**************************************************************/ 
template <class A_Type> void HASHCH<A_Type>::search(HASHCH<mybitset>& sh,unsigned int& countOk,unsigned int& countKo){

        countKo=countOk=0;
        for (unsigned int i=0;i<SizeHash;i++)
                {
                unsigned int tmpInext = HashVec[i];
                while(tmpInext !=DEFAULTP)
                        {
                         if (!(sh.search(ListCol[tmpInext]))) 
                         {
                           countKo++;
                         }
                         else
                         {
                           countOk++;
                         }
                        tmpInext= ListCol[tmpInext].Inext;
                        }
                HashVec[i]=DEFAULTP;
                }
        NumEl=0;
        };

/**************************************************************/
/* NAME :  Class HASHCH FreqPrint() */
/* DESCRIPTION :  printing the frequencies of the k-mers*/
/**************************************************************/ 

template <class A_Type> void HASHCH<A_Type>::FreqPrint(void){

unsigned int freq[1000];
memset(freq,0,1000*sizeof(unsigned int));

for (unsigned int i=0;i<NumEl;i++)
{
  if (ListCol[i].getFreq()<1000)
    freq[ListCol[i].getFreq()]++;
  else
    freq[999]++;
}
cout<<"=================================\n";
cout<<"Freq.\tNum.\t(Freq.)*(Num)\n";
for (int i=0;i<1000;i++)
{
  if (freq[i]!=0)
  cout<<i<<"\t"<< freq[i]<<"\t"<<freq[i]*i<<endl;
}
cout<<"=================================\n\n";
}


/**************************************************************/
/* NAME :  Class HASHCH FreqPrint() */
/* DESCRIPTION :  printing the frequencies of the k-mers*/
/**************************************************************/ 

template <class A_Type> void HASHCH<A_Type>::FreqKmerPrint(ofstream& fout,const int window){

char buffer[MAXSIZE];  
for (unsigned int i=0;i<NumEl;i++)
{
  ListCol[i].bit2buffer(buffer);
  if (DIM>window)
    buffer[window]='\0';
  //if (ListCol[i].getFreq()> TFREQ)
    fout<<buffer<<"\t"<<ListCol[i].getFreq()<<endl;
}

}
