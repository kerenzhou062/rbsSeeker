/********************************************************************
 * rbsSeeker: identify RBP binding sites at single-base resolution
 * Author: jianhua yang
 * Email: yangjh7@mail.sysu.edu.cn
 * Copyright: School of Life Sciences, Sun Yat-sen University
 * $ 2024/12/09
 ********************************************************************/
#ifndef rbsSeeker_HEAD_H
#define rbsSeeker_HEAD_H

#ifndef ENCODENUM
#define ENCODENUM 9
#endif

#define MINVAL 1E-300
#define MIN_RATIO 0.00001
#define EXTEND_LEN 10
#define SEEKER_TYPE 7
#define CLIP_TYPE 6
#define rpmFactor 1000000
#define mergePeakDist -50
#define mergeVarDist  10
#define pseudoCount 1
#define baseNum 5
#define nuclNum 15
#define peakExtLen 20
#define mutDist 3
#define peakExtWindowLen 500
#define minCtrlNum 10

typedef std::pair<uint32_t, uint32_t> exonPair;
typedef vector<exonPair> exonVector;

typedef vector<int> nuclVector;

typedef map<string, double*> varProfileMap;

struct rbsSumTypeInfo {
  double totalReadNum;
  double totalHeight;
  double varRatio;
  long int totalLen;
};

typedef struct rbsSumTypeInfo rbsSumType;

typedef map<string, rbsSumType*> rbsSumTypeMap;

struct rbsSiteInfo {
  //char *chrom;
  uint16_t chromCode;
  uint32_t start;
  uint32_t end;
  char     strand;
  double   tVarNum;
  double   tHeight;
  double   cVarNum;
  double   cHeight;
  double   mfold;
  double   ufold;
  double   dfold;
  double   pval;
  double   qval;
};

typedef struct rbsSiteInfo RBSite;

typedef vector<RBSite *> siteVector;

typedef map<string, siteVector> siteTypeMap;

typedef vector<double> windowVector;

typedef struct parameterInfo
{
  double minReadNum;
  double minHeight;
  double minRpm;
  double minT2cNum;
  double minMutNum;
  double pval;
  double qval;
  double mfold;
  double pfold;
  double minRatio;
  double maxRatio;
  double totalTreatNum;
  double totalCtrlNum;
  double treatVsCtrlRatio;
  double dropRatio;
  int    clipType;
  int    maxLocusNum;
  int    maxReadDist;
  int    minClusterLen;
  int    maxClusterLen;
  int    verbose;
  int    norm;
  int    PCR;
  int    rmSeMutation;
  int    peakModel;
  int    minReadLen;
  int    maxReadLen;
  int    bam;
  int    fullLength;
  int    barcodeLen;
  int    skipSplice;
  int    windowLen;
  int    rnafold;
  uint8_t cvsCode;
  long   genomeSize;
  long   transcriptomeSize;
  char  *cvs;
  char  *primerSeq;
  char  *motifSeq;
} parameterInfo;

void seekRbsSites(FILE *gfp, FILE *faifp, char *outdir, struct parameterInfo *paraInfo, char *treatFile, char *controlFile);

void callPeakAndMutation(struct parameterInfo *paraInfo, chromVarMap &treatVarHash,
                         chromVarMap &controlVarHash, FILE *gfp, faidxMap &faiHash, char *outdir);

void initiateRbsSumTypeMap(rbsSumTypeMap &rbsSumHash, int clipType);

void freeRbsSumTypeMap(rbsSumTypeMap &rbsSumHash);

void allocVarProfileMap(varProfileMap &varProHash, int span, int clipType);

void freeVarProfileMap(varProfileMap &varProHash);

double varToProfile(struct parameterInfo *paraInfo, varVector &varList, varProfileMap &varProHash,
                    char *chromSeq, int chromLen, char strand);

int getVariation(struct parameterInfo *paraInfo, varProfileMap &varProHash, Variation *var, int clipType, char strand);

void identifyRbsSites(struct parameterInfo *paraInfo, varProfileMap &treatVarProHash, varProfileMap &ctrlVarProHash,
                      siteTypeMap &allSiteTypeHash, char *chrom, int chromLen, char strand,
                      rbsSumTypeMap &treatRbsSumHash, rbsSumTypeMap &ctrlRbsSumHash);

int identifyPeaks(struct parameterInfo *paraInfo, double *treatHgt, double *ctrlHgt,
                  uint16_t chromCode, int chromLen, int peakStart, int peakEnd, char strand, 
                  double maxCov, int maxPos, double ctrlRatio, siteVector &siteList);

void calculatePval(struct parameterInfo *paraInfo, siteTypeMap &siteTypeHash,
                   rbsSumTypeMap &treatRbsSumHash, rbsSumTypeMap &ctrlRbsSumHash);

void calculateGenomeVariation(struct parameterInfo *paraInfo, rbsSumTypeMap &treatRbsSumHash);

void outputAllTypeRbsSite(struct parameterInfo *paraInfo, siteTypeMap &treatSiteTypeHash,
                          FILE *gfp, faidxMap &faiHash, char *outdir);

int outputPeakTypeSites(struct parameterInfo * paraInfo, FILE * gfp, faidxMap & faiHash,
                        siteVector & siteArray, char *outdir, char *siteType);

void mergeAllPeakSite(struct parameterInfo * paraInfo, siteVector & siteArray, siteVector & mergeSiteArray, char strand);

int outputOnePeakSite(struct parameterInfo * paraInfo, FILE * outfp, FILE * gfp, faidxMap & faiHash,
                      RBSite * rbs, int siteIdx, char *siteType);

void mergeAllVariationSite(struct parameterInfo * paraInfo, siteVector & siteArray, siteVector & mergeSiteArray, char strand);

int outputOneVarSite(struct parameterInfo * paraInfo, FILE * outfp, FILE * gfp, faidxMap & faiHash,
                     RBSite *rbs, int siteIdx, char *siteType);

void copyRBSiteNoChrom(RBSite * tRbs, RBSite * oRbs);

void copyRBSite(RBSite * tRbs, RBSite * oRbs);

int outputVariationTypeSites(struct parameterInfo * paraInfo, FILE * gfp, faidxMap & faiHash,
                             siteVector & siteArray, char *outdir, char *siteType);

void freeSiteTypeMap(siteTypeMap &siteTypeHash) /*free site map */;

RBSite *allocRBSite(uint16_t chromCode, int start, int end, char strand);

char revComChar(char ch, char strand);

int searchMotif(char *motif, char *seq);

void freeSite(RBSite *rbs);

void freeSiteVector(siteVector &sList);

void freeParameters(parameterInfo *paraInfo);

void getSiteBhFdr(siteVector &BHhash);

int cmpChromLocus(const RBSite *x, const RBSite *y);

int encodeNucleotide(char base);

void decodeNucleotide(char base, nuclVector &nuclList);

#endif /* End rbsSeeker_HEAD_H */

