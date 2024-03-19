CXXC=g++
LIBS=-lm -lz
CFLAGS = -O3 -g
HG_DEFS = -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_GNU_SOURCE
HG_WARN=-Wformat -Wreturn-type
UTILITIES_DIR = ./thirdUtils
BIN_DIR = ./bin
BIO_DIR = ./bioUtils
BAM_DIR = ./thirdUtils/BamTools
RNAFOLD_DIR = ./thirdUtils/RNAfoldLib
INCLUDES = -I$(UTILITIES_DIR)/BamTools/include \
           -I$(UTILITIES_DIR)/BamTools/include/api \
           -I$(UTILITIES_DIR)/RNAfoldLib \
           -I$(BIO_DIR)
BIO_LIBS   = -L$(UTILITIES_DIR)/BamTools/lib/ -lbamtools \
             -L$(UTILITIES_DIR)/RNAfoldLib/ -lRNAfold \
             -L$(BIO_DIR)/ -lbiotools

all:
	cd $(BAM_DIR); make api; make
	cd $(BIO_DIR); make
	cd $(RNAFOLD_DIR); make
	make rbsSeeker

clean:
	cd $(BAM_DIR); make clean_api
	cd $(BIO_DIR); make clean
	cd $(RNAFOLD_DIR); make clean
	rm -f *.o

rbsSeeker: rbsSeeker.o rbsSeekerMain.o
	$(CXXC) $(CFLAGS) ${HG_DEFS} ${HG_WARN} $(INCLUDES) -o ${BIN_DIR}/rbsSeeker rbsSeekerMain.o rbsSeeker.o \
	$(BIO_LIBS) $(LIBS) 

rbsSeeker.o: rbsSeeker.cpp rbsSeeker.h
	$(CXXC) $(CFLAGS) ${HG_DEFS} ${HG_WARN} $(INCLUDES) -c rbsSeeker.cpp
	
rbsSeekerMain.o: rbsSeekerMain.cpp rbsSeeker.h
	$(CXXC) $(CFLAGS) ${HG_DEFS} ${HG_WARN} $(INCLUDES) -c rbsSeekerMain.cpp
	
