CIRI-AS is a de novo detection tool for circRNA internal components and alternative splicing events based on spliced junction signatures, especially from back-spliced junction read pairs. Original codes and detailed manuals are included in the downloading URL: https://sourceforge.net/projects/ciri/files/CIRI-AS/. 
If you use CIRI-AS, please cite the following paper: 
Yuan Gao†, Jinfeng Wang†, Yi Zheng†, Jinyang Zhang, Shuai Chen and Fangqing Zhao*. Comprehensive identification of internal structure and alternative splicing events in circular RNAs. Nature Communications, 2016, DOI:10.1038/ncomms12060.

Q&A: 
Q: An error "Please make sure the SAM is the same one used for circRNA detection by CIRI without modification!" when running CIRI-AS after CIRI2?
A: Because CIRI-AS was developed before CIRI2, it cannot process sequencing reads with different lengths whereas CIRI2 can. It seems that your sequencing reads were trimmed to different lengths before mapping, so the output of CIRI2 based on these reads cannot be recognized by CIRI-AS. My suggestion is to use raw reads directly, or trim all reads to the same length before BWA, CIRI2 and CIRI-AS.