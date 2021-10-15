/* bed head file */

#ifndef rbsSeeker_HEAD_H
#define rbsSeeker_HEAD_H

#ifndef ENCODENUM
#define ENCODENUM 9
#endif

#define MINVAL 1E-300

#define EXTEND_LEN 10

extern double transMatrix[ENCODENUM][ENCODENUM];
extern double transMatrixNum[ENCODENUM][ENCODENUM];

struct peakInfo {
  char   *chrom;
  int    start;
  int    end;
  char   strand;
  double readNum;
  double totalHeight;
  double totalConversionNum;
  double totalConversionHeight;
  double totalMutationNum;
  double totalMutationHeight;
  double totalDeletionNum;
  double totalDeletionHeight;
  double totalTruncationNum;
  double totalTruncationHeight;
  double totalInsertionNum;
  double totalInsertionHeight;
  double totalEndNum;
  double totalEndHeight;
  double t2cRatio;
  double mutRatio;
  double truRatio;
  double delRatio;
  double insRatio;
  double endRatio;
  double pval;
  double qval;
  double mfold;
  double maxHeight;
  double heightRpm;
  double ratio;
  int readLen; // sum of all read length
  int maxHeightPos;
  char *seq;
  double **profile;
};

typedef struct peakInfo Cluster;

typedef vector<Cluster *> clusterVector;

struct genomePeakInfo {
  double totalReadNum;
  double totalHeight;
  double totalConversionNum;
  double totalConversionHeight;
  double totalMutationNum;
  double totalMutationHeight;
  double totalDeletionNum;
  double totalDeletionHeight;
  double totalTruncationNum;
  double totalTruncationHeight;
  double totalInsertionNum;
  double totalInsertionHeight;
  double totalEndNum;
  double totalEndHeight;
  double t2cRatio;
  double mutRatio;
  double truRatio;
  double delRatio;
  double insRatio;
  double endRatio;
  long int totalLen;
};

typedef struct genomePeakInfo genomeCluster;

struct clusterSiteInfo {
  char   *chrom;
  int    start;
  int    end;
  char   strand;
  double heightRpm;
  double readNum;
  double height;
  double ratio;
  double mfold;
  double pval;
  double qval;
  char *type;
};

typedef struct clusterSiteInfo RBSite;

typedef vector<RBSite *> siteVector;

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
  double minRatio;
  double maxRatio;
  int    maxLocusNum;
  int    maxReadDist;
  int    minClusterLen;
  int    verbose;
  int    PCR;
  int    rmSeMutation;
  int    peakModel;
  int    minReadLen;
  int    maxReadLen;
  int    bam;
  int    fullLength;
  int    barcodeLen;
  long   genomeSize;
  long   transcriptomeSize;
  char  *cvs;
  char  *primerSeq;
} parameterInfo;

struct bedName
{
  struct bedName *next;
  char *readName;
};

void seekRbsSites(FILE *gfp, FILE *faifp, char *outdir, struct parameterInfo *paraInfo, char *bamFile);

double normalizedReads(struct parameterInfo *paraInfo, chromSamMap &samHash);

void callPeaks(struct parameterInfo *paraInfo, chromSamMap &bedHash,
               clusterVector &contigVector, FILE *gfp, faidxMap &faiHash);

int outputTypeSites(struct parameterInfo *paraInfo, char *outdir, FILE *gfp, faidxMap &faiHash, siteVector &siteArray);

int outputPeakSites(struct parameterInfo *paraInfo, char *outdir, FILE *gfp, faidxMap &faiHash, clusterVector &contigVector);

void outputSiteInfo(FILE *gfp, faidxMap &faiHash, char *outdir, struct parameterInfo *paraInfo, double totalReadNum, clusterVector &contigVector);

double getGenomeClusterInfo(struct parameterInfo *paraInfo, clusterVector &contigVector, genomeCluster *gCluster);

double **getProfile(struct parameterInfo *paraInfo, samVector &samList,
                    int start, int end, char strand, Cluster *pCluster);

void getAllMutations(struct parameterInfo *paraInfo, char strand, int start, int end,
                     samVector &samList, clusterVector &contigVector, FILE *gfp, faidx *fai);

int getPeakInfo(struct parameterInfo *paraInfo, char strand, samVector &samList,
                clusterVector &contigVector, FILE *gfp, faidx *fai);

double **alloc2array(int nrows, int ncolumns);

void free2array(double **array, int nrows);

int encodeAlignment(char ch, char strand);

void getCIMS(double **readProfile, char *seq, char *cigar, double readNum, char strand);

char decodeAlignment(int i);

void identifyClusterInfo(char *clusterSeq, char strand, int start, int end, double **profile, Cluster *pCluster, struct parameterInfo *paraInfo);

int identifyClusterSites(struct parameterInfo * paraInfo, FILE * gfp,
                         faidxMap & faiHash, double coverageRead,
                         Cluster * pCluster, genomeCluster * gCluster,
                         siteVector & conversionArray, siteVector & mutationArray,
                         siteVector & truncationArray, siteVector & peakArray,
                         siteVector & endArray, siteVector & delArray, siteVector & insArray);

int outputPeakHeightTypeSites(struct parameterInfo *paraInfo, char *outdir, FILE *gfp, faidxMap &faiHash, siteVector &siteArray);

int outputStrandRbsSite(struct parameterInfo *paraInfo, FILE *outfp, 
  FILE *gfp, faidxMap &faiHash, siteVector &siteArray, int siteIdx, char strand);

int outputRbsSite(struct parameterInfo *paraInfo, FILE *outfp, FILE *gfp, faidxMap &faiHash, RBSite *rbs, int siteIdx);

int searchMotif(char *seq);

void freeCluster(Cluster *clust);

void freeClusterVector(clusterVector &cList);

void copyRBSite(RBSite *tRbs, RBSite *oRbs);

void freeSite(RBSite *rbs);

void freeSiteVector(siteVector &sList);

void freeParameters(parameterInfo *paraInfo);

void getSiteBhFdr(siteVector &BHhash);

void getClusterBhFdr(clusterVector &BHhash);

int cmpChromLocus(const RBSite *x, const RBSite *y);

double removeMisprimingReads(struct parameterInfo* paraInfo,
                             FILE* genomefp,
                             faidxMap &faiHash,
                             chromSamMap &samHash,
                             chromSamMap &newSamHash);

void initTransMatrix(void);

void outputTransMatrix(struct parameterInfo *paraInfo, char *outdir);

void removeSeMutations(double **profile, double **readProfile, CSam *sam, int start, char strand);

#endif /* End rbsSeeker_HEAD_H */
