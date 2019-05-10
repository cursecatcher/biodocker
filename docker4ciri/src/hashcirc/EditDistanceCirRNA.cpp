#include <string>
#include <iostream>
#include <vector>
#include <map>    
#include <fstream>

#include "ssw_cpp.h"
#include "ssw.h"

#ifndef __CON_H__
	#define __CON_H__
	#include "const.h"
#endif

#ifndef __GEN_H__
	#define __GEN_H__
	#include "general.h"
#endif


using namespace std;
using std::string;
using std::cout;
using std::endl;
 
const std::string* idRef {nullptr};
const std::string* idRead  {nullptr};
 
double  d[2][MAXSIZE];
int input_matches {0};

#define  BLAST_OUTPUT 0
#define  BEST_MATCH 1

#if BLAST_OUTPUT
static void ssw_write (StripedSmithWaterman::Alignment& a,
			const string& ref_seq,
			const string& read_seq) {
	static int count=0;
	cout<<"\n\n\n________________________________________________________________________________\n";
	cout<<"===================================MAPPING "<<++count<<"====================================\n"<<endl;
	
	cout<<"Ref:\n";
	cout<<*idRef;
	for (unsigned int i=0;i<ref_seq.size();++i){
	  if (i%80==0)
	    cout<<"\n";
	  cout<<ref_seq[i];
	}
	cout<<"\n\nRead:\n";
	cout<<*idRead;
	for (unsigned int i=0;i<read_seq.size();++i){
	  if (i%80==0)
	    cout<<"\n";
	  cout<<read_seq[i];
	}
	       
	cout<<"\n\n================================Info alignment==================================\n"<<endl;
	fprintf(stdout, "\toptimal_alignment_score: %d\tsub-optimal_alignment_score: %d\t\n", a.sw_score, a.sw_score_next_best);
	if (a.ref_begin + 1) fprintf(stdout, "\ttarget_begin: %d\t", a.ref_begin + 1);
	fprintf(stdout, "target_end: %d\t\n", a.ref_end + 1);
	if (a.query_begin + 1) fprintf(stdout, "\tquery_begin: %d\t", a.query_begin + 1);
	fprintf(stdout, "\tquery_end: %d\n\n", a.query_end + 1);
	cout<<"\n===================================Alignment====================================\n"<<endl;
	if (a.cigar.size()>0) {
		int32_t c = 0, left = 0, e = 0, qb = a.ref_begin, pb = a.query_begin;
		uint32_t i;
		int cigarLen=a.cigar.size();
		while (e < cigarLen || left > 0) {
			int32_t count = 0;
			int32_t q = qb;
			int32_t p = pb;
			fprintf(stdout, "Target: %8d    ", q + 1);
			for (c = e; c < cigarLen; ++c) {
				char letter = cigar_int_to_op(a.cigar[c]);
				uint32_t length = cigar_int_to_len(a.cigar[c]);
				uint32_t l = (count == 0 && left > 0) ? left: length;
				for (i = 0; i < l; ++i) {
					if ((letter == 'I')||(letter == 'S')) fprintf(stdout, "-");
					else {
						fprintf(stdout, "%c", ref_seq[q]);
						++ q;
					}
					++ count;
					if (count == 60) goto step2;
				}
			}
step2:
			fprintf(stdout, "    %d\n                    ", q);
			q = qb;
			count = 0;
			for (c = e; c < cigarLen; ++c) {
				char letter = cigar_int_to_op(a.cigar[c]);
				uint32_t length = cigar_int_to_len(a.cigar[c]);
				uint32_t l = (count == 0 && left > 0) ? left: length;
				for (i = 0; i < l; ++i){
					if (letter == '=') {
						fprintf(stdout, "|");
						++q;
						++p;
					} else {
						
						fprintf(stdout, "*");
						if ((letter == 'I')||(letter == 'S')) ++p;
						else ++q;
					}
					++ count;
					if (count == 60) {
						qb = q;
						goto step3;
					}
				}
			}
step3:
			p = pb;
			fprintf(stdout, "\nQuery:  %8d    ", p + 1);
			count = 0;
			for (c = e; c < cigarLen; ++c) {
				char letter = cigar_int_to_op(a.cigar[c]);
				uint32_t length = cigar_int_to_len(a.cigar[c]);
				uint32_t l = (count == 0 && left > 0) ? left: length;
				for (i = 0; i < l; ++i) {
					if ((letter == 'D')||(letter == 'S')) fprintf(stdout, "-");
					else {
						fprintf(stdout, "%c", read_seq[p]);
						++p;
					}
					++ count;
					if (count == 60) {
						pb = p;
						left = l - i - 1;
						e = (left == 0) ? (c + 1) : c;
						goto end;
					}
				}
			}
			e = c;
			left = 0;
end:
			fprintf(stdout, "    %d\n\n", p);
		}
	}
        cout<<"________________________________________________________________________________\n";
	cout<<"================================================================================"<<endl;
}

#endif 


    // Declares a default Aligner
  StripedSmithWaterman::Aligner aligner(2,4,100,1);
  // Declares a default filter
  StripedSmithWaterman::Filter filter;
  // Declares an alignment that stores the result
  StripedSmithWaterman::Alignment alignment;
  // Aligns the query to the ref

double SmithWaterman(const string& S,const string& T,double& matches){


  
  aligner.Align(T.c_str(), S.c_str(), S.size(), filter, &alignment);

  //compute matches
  int csize= alignment.cigar.size();
  matches=0;
  double tot=0;
  for (int i=0;i<csize;++i){
    if(cigar_int_to_op(alignment.cigar[i])=='='){
      matches+=(double)cigar_int_to_len(alignment.cigar[i]);
      tot+=(double)cigar_int_to_len(alignment.cigar[i]);
      }
   // tot+=(double)cigar_int_to_len(alignment.cigar[i]);
    else
      if(cigar_int_to_op(alignment.cigar[i])=='X')
	tot+=(double)cigar_int_to_len(alignment.cigar[i]);
	
  }

  
#if BLAST_OUTPUT
  if (matches>input_matches)  
    ssw_write(alignment,S,T);
#endif 
  
  return alignment.sw_score;
}


 inline void reverse(string& S,string& r_S){
 int size=S.size();   
 for(int i=size-1;i>=0;--i){
      switch(S[i]){
	case 'A':
	case 'a':  
	 r_S.push_back('T');
	break;
	case 'T':
	case 't':  
	  r_S.push_back('A');
	break;
	case 'C':
	case 'c': 
	  r_S.push_back('G');
	break;     
	case 'G':
	case 'g':  
	  r_S.push_back('C');
	break;
	case 'N':
	case 'n':  
	  r_S.push_back('N');
	break;
	default:
	  cerr<<"READ \n"<<S<<"\n contains an character not allowed\n";
	  exit(EXIT_FAILURE);
      }
    }

}
    
    
    int main(int argc, char **argv)
    {
        cout<<"\n\n ================================================================================\n";
        cout<<"|	      	         ALIGNMENT TOOL  For CirRNA          		         |\n";
        cout<<"|	      	      		                                                 |\n";
        cout<<" ================================================================================\n";
        cout<<"\n If you find any bug, send an email to beccuti@di.unito.it\n";
        
        
        if (argc<5)
        {
            std::cerr<<"\n\nUSE: AlignmentBlast <junction_file> <read_file_1> <matches> <output_file>\n\n"<<endl;
            std::cerr<<"\tNB: Input files must be in FastQ format\n\n";
            exit(EXIT_FAILURE);
        }
        
        
        
    
        clock_t startGlobal,endGlobal;
        startGlobal=clock();
        
        map < string, pair<string,int> > ref;
    
        ifstream f_ref(argv[1],ifstream::in); 
         if(!f_ref) {
                cerr << "\n*****Error opening input file "<< argv[1] <<" *****" << endl;
                exit(EXIT_FAILURE);
            }
        std::string id("");

        while (!f_ref.eof())
        { 
                std::string buffer("");
                getline(f_ref,buffer);
               // cout<<buffer<<endl;
                if ((buffer[0]=='>')||(buffer[0]=='@')||(buffer!="")){
                    id=buffer; 
                    getline(f_ref,buffer);
                 //   cout<<buffer<<endl;
                    if (buffer[0]=='A'||buffer[0]=='a' ||buffer[0]=='C'||buffer[0]=='c'||buffer[0]=='G'||buffer[0]=='g'||buffer[0]=='T'||buffer[0]=='t')
                    {
                        ref[id]=(make_pair(buffer,0));
                        getline(f_ref,buffer); //+
                        getline(f_ref,buffer); //quality
                    }
                    else
                    {
                        
                        cerr << "\n*****Error in the format of reference file "<< argv[1] <<" *****" << endl;
                        exit(EXIT_FAILURE);
                    }
                }
        }
        f_ref.close();
        cout<<"\nTot. junction reads: "<<ref.size()<<endl;
   
        cout<<"\nReading  read file: "<< argv[2]<<endl; 
        ifstream f_clone(argv[2],ifstream::in);
        if(!f_clone) {
            cerr << "\n*****Error opening input file "<< argv[2] <<" *****" << endl;
            exit(EXIT_FAILURE);
        } 
        double max_matches=atoi(argv[3]);
        cout<<"\nStart Mapping: "<<endl; 
                int proc_reads=0;
        while (!f_clone.eof()){//for all reads
            std::string id("");
            getline(f_clone,id);
            std::string clone("");
            getline(f_clone,clone);
            if ((id!="")&&(clone!="")&&  ((id[0]=='>')||(id[0]=='@'))&& (clone[0]=='A'||clone[0]=='a' ||clone[0]=='C'||clone[0]=='c'||clone[0]=='G'||clone[0]=='g'||clone[0]=='T'||clone[0]=='t')){
                std::string cloneR(""); 
                reverse(clone,cloneR); 
                //scores
               // double Maxscore =0, 
                //matches
                
                input_matches=max_matches;
                //positions
                map< string, pair<string,int> >::iterator max_position = ref.end();
               // set< pair<string,string> >::iterator Maxit[num_file];
                //set< pair<string,string> >::iterator MaxitR[num_file];
                
#if BEST_MATCH                
                for (auto it=ref.begin();(it!=ref.end());++it  ){
#else
                bool exitFor=false;
                for (auto it=ref.begin();((it!=ref.end())&&(!exitFor));++it  ){   
#endif                    
                        double matches=0,matchesR=0;
                        double score=0, scoreR=0;
#if BLAST_OUTPUT                        
                        idRef=&(id);
                        idRead=&(it->first);
#endif                        
                        score=SmithWaterman(clone,it->second.first,matches);
                        scoreR=SmithWaterman(cloneR,it->second.first,matchesR);
                      
                        if (matchesR > matches){
                        //    score=scoreR;
                            matches=matchesR;
                        }
                     //    cout<<"Matches "<<matches<<endl;    
                        if (matches>input_matches)
                        {
                            //Maxscore=score;//*100/(2*(min(clone.size(),it->second.size())));
                            max_position=it;
                            input_matches=matches;
#if !BEST_MATCH                           
                            exitFor=true;
#endif                            
                        } 
                    }                    
                if (max_position !=ref.end()){
                    ++(max_position->second.second);
                    //cout<<id<<" "<<max_position->first<<endl;
                }
                if ((++proc_reads%1000==1)){
                    cout<<"\tProcessed reads: " <<proc_reads<<endl;
                    cout<<"\tTime for read: "<<((double)(clock()-startGlobal))/(CLOCKS_PER_SEC*proc_reads)<<endl<<endl;
#if DEBUG
                    int sum=0;          
                    for (auto it=ref.begin();it!=ref.end();++it  ){
                        if (it->second.second>0)
                        {
                            cout<<it->first<<endl<<"\t"<<it->second.second<<endl;
                            sum+=it->second.second;
                        }
                       
                    }
                 cout<<"\tSummary:"<<sum<<endl; 
#endif                 
                }
            }
        }
        cout<<"\nTot.  reads: "<<proc_reads<<endl;
        f_clone.close();
        cout<<"\n================================================================================\n\n";
        
         ofstream out(argv[4],ofstream::out); 
         if(!out) {
                cerr << "\n*****Error opening output file "<< argv[4] <<" *****" << endl;
                exit(EXIT_FAILURE);
            }
          cout<<"\n\n ===============================================================================";  
                cout<<"\n|		      	    Reference Counting	                                |\n"; 
                cout<<" ===============================================================================\n\n";  
                for (auto it=ref.begin();it!=ref.end();++it  ){
                     cout<<it->first<<endl<<"\t"<<it->second.second<<endl;
                     out<<it->first<<"\t"<<it->second.second<<endl;
                 }
        out.close();         
        endGlobal=clock();
        
        cout<<"\n=====================================TIME=======================================\n\n";
        cout<<"\tTime to identify alignment: "<<((double)(endGlobal-startGlobal))/CLOCKS_PER_SEC<<"s."<<endl;
        cout<<"\n================================================================================\n\n";
        
    }
