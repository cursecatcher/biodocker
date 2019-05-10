OPTION =   -O3 -Wall  -fmessage-length=0 -fopenmp  -std=c++11
CC = gcc
C++ = g++


all:	 AlignmentCirRNA HashCheckerFilter

clean:
	rm   AlignmentCirRNA HashCheckerFilter

install:
	cp AlignmentCirRNA HashCheckerFilter /usr/bin


AlignmentCirRNA:	  EditDistanceCirRNA.cpp ssw_cpp.cpp ssw.c  ssw_cpp.h ssw.h const.h
	$(C++)  -o AlignmentCirRNA EditDistanceCirRNA.cpp ssw_cpp.cpp ssw.c  $(OPTION)

HashCheckerFilter:	main.cpp classIO.cpp general.cpp classIO.h const.h conf.h classHashChaining.h classHashChaining.cpp classMybitset.h classMybitset.cpp
	$(C++) -D SIMILARITY=1 -o  HashCheckerFilter classIO.cpp  main.cpp  general.cpp     $(OPTION)


doc:
	doxygen  DoxyFile
