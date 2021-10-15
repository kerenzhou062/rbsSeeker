CXXC=g++
LIBS=-lm -lz
CFLAGS = -O3 -g
HG_DEFS = -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_GNU_SOURCE
HG_WARN=-Wformat -Wreturn-type
UTILITIES_DIR = ../../thirdUtils
BIN_DIR = ../../bin
BIO_DIR = ../../bioUtils
INCLUDES = -I$(UTILITIES_DIR)/BamTools/include \
           -I$(UTILITIES_DIR)/BamTools/include/api \
           -I$(UTILITIES_DIR)/cdflib \
           -I$(UTILITIES_DIR)/alglibsrc \
           -I$(UTILITIES_DIR)/RNAfoldLib \
           -I$(UTILITIES_DIR)/RNAshapesLib \
           -I$(BIO_DIR)
BIO_LIBS   = -L$(UTILITIES_DIR)/BamTools/lib/ -lbamtools \
             -L$(UTILITIES_DIR)/RNAfoldLib/ -lRNAfold \
             -L$(BIO_DIR)/ -lbiotools \
             -L$(UTILITIES_DIR)/cdflib/ -lcdf \

all:
	make ../../thirdUtils/RNAfoldLib
	make ../../thirdUtils/BamTools
	make ../../thirdUtils/cdflib
	make rbsSeeker

rbsSeeker: rbsSeeker.o rbsSeekerMain.o
	$(CXXC) $(CFLAGS) ${HG_DEFS} ${HG_WARN} $(INCLUDES) -o ${BIN_DIR}/rbsSeeker rbsSeekerMain.o rbsSeeker.o \
	$(BIO_LIBS) $(LIBS) 

rbsSeeker.o: rbsSeeker.cpp rbsSeeker.h
	$(CXXC) $(CFLAGS) ${HG_DEFS} ${HG_WARN} $(INCLUDES) -c rbsSeeker.cpp
	
rbsSeekerMain.o: rbsSeekerMain.cpp rbsSeeker.h
	$(CXXC) $(CFLAGS) ${HG_DEFS} ${HG_WARN} $(INCLUDES) -c rbsSeekerMain.cpp
	
clean:
	rm -f *.o
