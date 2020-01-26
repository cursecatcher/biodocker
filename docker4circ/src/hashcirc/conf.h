
//!If 1 then it  sets the output in DEBUG format
//#define DEBUGOUTPUT 0

//!If 1 then it prints frequencies of the k-mer in the hash table. If 2 then  it prints k-mers and their frequencies (X leukemia)
//#define DEBUGFREQ 0


//!if 1 then it prints the similar reads
//#define SIMILARITY 0


//! if 1 the it returns the signatures shared by a set of reads (i.e. the greatest sequence of continuous k-mers present in the Hash Table).
//#define SIGNATURE 0

//#if SIGNATURE==0
//! if 1 then it checks the presence of the two mates in Hash Table. It returns the paired-end reads with at least one mate not present in the Hash Table
  //#define PAIREDEND 0
//#endif
//! if 1 then output od CompCheckerRead is generated for the GUI
#define GUI 1

//! if we want to normalized the frequencies according to the number of N bases present in the read. It is used only in the selection of candidate for each clone.
//! Then the choice is done on the normalized frequencies.
#define  Nnormalized 1  