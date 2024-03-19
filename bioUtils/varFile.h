#ifndef VARIATION_HEAD_H
#define VARIATION_HEAD_H

struct variationInfo
{
	struct variationInfo *next; // for list
	char *chrom;	// chromosome name
	uint32_t chromStart; // chromosome start
	uint32_t chromEnd; // chromosome end
	double readNum; // normalized read number
	char strand; // strand
	uint32_t mutNum; // mutation number
	uint32_t *position; // muations positions
	uint8_t *mutation; // muations
};

typedef std::pair<uint32_t, uint8_t> varPair;

typedef struct variationInfo Variation;

typedef vector<varPair> variationVector;

typedef vector<Variation *> varVector;

typedef map<string, varVector> chromVarMap;

char revComBase(char ch, char strand);

void getOneVariation(char *chrom, char *chromSeq, char *seq, char *cigar, Variation *var);

uint8_t encodeVariation(char a, char b);

void freeChromVarMap(chromVarMap &varHash) /*free var map */;

void freeVarVector(varVector &varList);

void freeVariationItem(Variation *var);

bool compareVar(const Variation *oVar, const Variation *tVar);

#endif /* End VARIATION_HEAD_H */
