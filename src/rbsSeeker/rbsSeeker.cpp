/* API for sam format */
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include<assert.h>
#include<math.h>
#include<time.h>
#include<limits.h>
extern "C" {
#include "fold.h"
#include "fold_vars.h"
#include "utils.h"
#include "pair_mat.h"
}
#include "BamReader.h"
#include "BamAux.h"
using namespace BamTools;
#include <map>
#include <algorithm>
#include <ios>
#include <iostream>
#include <string>
#include <ostream>
#include <fstream>
#include <iomanip>
#include <locale>
#include <sstream>
#include <vector>

using namespace std;

#include "bioUtils.h"
#include "faiFile.h"
#include "bedFile.h"
#include "samFile.h"
#include "homer_statistics.h"
#include "statistic.h"
#include "rbsSeeker.h"

int cmpSite(const RBSite *x, const RBSite *y)
{
  return x->pval < y->pval;
}

int cmpChromLocus(const RBSite *x, const RBSite *y)
{
  if (strcmp(x->chrom, y->chrom) != 0)
    return strcmp(x->chrom, y->chrom) < 0;

  return x->start < y->start;
}

int cmpCluster(const Cluster *x, const Cluster *y)
{
  return x->pval < y->pval;
}

void getClusterBhFdr(clusterVector &BHhash) {
  sort(BHhash.begin(), BHhash.end(), cmpCluster);
  int rank = 1;
  int increaseRank = 1;
  double lastVal = 1000;
  int arrayNum = BHhash.size();
  for (clusterVector::iterator curr = BHhash.begin(); curr != BHhash.end(); curr++)
  {
    Cluster *clust = *curr;
    double pval = clust->pval;
    if (lastVal != pval)
    {
      rank = increaseRank;
    }
    increaseRank++;
    lastVal = pval;
    double FDR = pval + log10Val(arrayNum) - log10Val(rank);
    if (FDR > 0) {
      FDR = 0;
    }
    clust->qval = FDR;
  }
}

void getSiteBhFdr(siteVector &BHhash) {
  sort(BHhash.begin(), BHhash.end(), cmpSite);
  int rank = 1;
  int increaseRank = 1;
  double lastVal = 1000;
  int arrayNum = BHhash.size();
  for (siteVector::iterator curr = BHhash.begin(); curr != BHhash.end(); curr++) {
    RBSite *rbs = *curr;
    double pval = rbs->pval;
    if (lastVal != pval)
    {
      rank = increaseRank;
    }
    increaseRank++;
    lastVal = pval;
    double FDR = pval + log10Val(arrayNum) - log10Val(rank);
    if (FDR > 0) {
      FDR = 0;
    }
    rbs->qval = FDR;
  }
}


template<typename T>
string NumberToString(T Number)
{
  ostringstream ss;
  ss << Number;
  return ss.str();
}

void seekRbsSites(FILE *gfp, FILE *faifp, char *outdir, struct parameterInfo *paraInfo, char *bamFile)
{
  int mapLineNum = 0;
  int skipSplice = 1;
  double totalNum = 0;
  int keepDup = 1;
  if (paraInfo->PCR) keepDup = 0;

  faidxMap faiHash;
  chromSamMap samHash;
  chromSamMap newSamHash;
  clusterVector contigVector;

  BamReader reader;
  openBamFile(bamFile, reader);

  time_t start = 0, end = 0, begin = 0;
  time(&start);
  time(&begin);

  fprintf(stderr, "initialized the transmition matrix\n");
  initTransMatrix();

  fprintf(stderr, "read genome fai file\n");
  paraInfo->genomeSize = readFai(faifp, faiHash);
  if (paraInfo->transcriptomeSize > paraInfo->genomeSize) paraInfo->transcriptomeSize = paraInfo->genomeSize;
  time(&end);
  fprintf(stderr, "elasped %.2f seconds for reading genome fai file\n", (double)(end - start));

  time(&start);
  fprintf(stderr, "read bam to sam list\n");
  totalNum = readBamToSamMapNew(reader, samHash, keepDup, paraInfo->maxLocusNum, skipSplice);
  time(&end);
  fprintf(stderr, "elasped %.2f seconds for reading %.0f reads to bed list\n", (double)(end - start), totalNum);

  if (paraInfo->primerSeq != NULL)
  {
    time(&start);
    fprintf(stderr, "remove mispriming reads from sequencing data\n");
    totalNum = removeMisprimingReads(paraInfo, gfp, faiHash, samHash, newSamHash);
    fprintf(stderr, "elasped %.2f seconds for removing the mispriming reads, and remain reads=%.f \n", (double)(end - start), totalNum);
    time(&end);

    time(&start);
    fprintf(stderr, "normalized the read list\n");
    totalNum = normalizedReads(paraInfo, newSamHash);
    time(&end);
    fprintf(stderr, "elasped %.2f seconds for normalizing %.0f reads\n", (double)(end - start), totalNum);

    time(&start);
    fprintf(stderr, "identify contigs\n");
    callPeaks(paraInfo, newSamHash, contigVector, gfp, faiHash);
    time(&end);
    fprintf(stderr, "elasped %.2f seconds for identifying contigs\n", (double)(end - start));
  }
  else
  {
    time(&start);
    fprintf(stderr, "normalized the read list\n");
    totalNum = normalizedReads(paraInfo, samHash);
    time(&end);
    fprintf(stderr, "elasped %.2f seconds for normalizing %.0f reads\n", (double)(end - start), totalNum);

    time(&start);
    fprintf(stderr, "identify contigs\n");
    callPeaks(paraInfo, samHash, contigVector, gfp, faiHash);
    time(&end);
    fprintf(stderr, "elasped %.2f seconds for identifying contigs\n", (double)(end - start));
  }

  time(&start);
  fprintf(stderr, "output binding sites\n");
  outputSiteInfo(gfp, faiHash, outdir, paraInfo, totalNum, contigVector);
  time(&end);
  fprintf(stderr, "elasped %.2f seconds for outputing binding sites\n", (double)(end - start));

  reader.Close();

  fprintf(stderr, "free all memories\n");
  freeClusterVector(contigVector);
  freeChromSamMap(samHash);
  time(&end);
  fprintf(stderr, "total elasped %.2f seconds for whole workflow\n", (double)(end - begin));
}

double removeMisprimingReads(struct parameterInfo* paraInfo,
                             FILE* genomefp,
                             faidxMap &faiHash,
                             chromSamMap &samHash,
                             chromSamMap &newSamHash)
{
  int i = 0;
  int j = 0;
  int primerNum = 0;
  double totalPrimerNum = 0;
  double totalNum = 0;
  int primerSeqLen = strlen(paraInfo->primerSeq);
  chromSamMap::iterator it;

  for (it = samHash.begin(); it != samHash.end(); ++it)
  {
    samVector samList = it->second;
    char* chrom = (char*)it->first.c_str();
    string chromStr(chrom);
    if (faiHash.find(chromStr) == faiHash.end())
    {
      fprintf(stderr, "can't not find the chromosome %s, skip it.\n", chrom);
      continue;
    }
    faidx* fai = faiHash[chromStr];
    if (samList.size() >= 1)
    {
      for (samVector::iterator vecItr = samList.begin(); vecItr != samList.end(); vecItr++)
      {
        CSam *sam = *vecItr;
        int start = sam->chromStart;
        int end = sam->chromEnd;
        int primerStart = end + paraInfo->barcodeLen;
        int primerEnd = primerStart + primerSeqLen;
        if (sam->strand == '-')
        {
          primerEnd = start - paraInfo->barcodeLen;
          primerStart = primerEnd - primerSeqLen;
        }
        if (primerStart < 0) primerStart = 0;
        if (primerEnd < 0) primerEnd = 0;
        int matchNum = 0;
        char* primerSeq = faidxFetchSeq(genomefp, fai, primerStart, primerEnd, sam->strand);
        for (i = 0; i < primerSeqLen; i++)
        {
          if (primerSeq[i] == paraInfo->primerSeq[i])
          {
            matchNum++;
          }
        }
        if (matchNum >= primerSeqLen - 1)
        {
          primerNum++;
          totalPrimerNum += sam->readNum;
          if (paraInfo->verbose) fprintf(stderr, "mispriming read%d: %s %c\n", primerNum, sam->readName, sam->strand);
        }
        else {
          newSamHash[chromStr].push_back(sam);
          totalNum += sam->readNum;
        }
        safeFree(primerSeq);
      } // for end
    }
  } // for sam hash
  fprintf(stderr, "remove the mispriming reads: uniqueRead=%d and allRead=%.0f\n", primerNum, totalPrimerNum);
  return totalNum;
}

double normalizedReads(struct parameterInfo *paraInfo, chromSamMap &samHash)
{
  double totalNum = 0;
  int maxReadNum = 0;
  chromSamMap::iterator it;
  map<string, int> mapReads;

  for (it = samHash.begin(); it != samHash.end(); ++it)
  {
    samVector samList = it->second;
    if (samList.size() >= 1)
    {
      for (samVector::iterator vecItr = samList.begin(); vecItr != samList.end(); vecItr++)
      {
        CSam *sam = *vecItr;
        mapReads[sam->readName] = 0;
      } // for end
    }
  } // for bed hash

  for (it = samHash.begin(); it != samHash.end(); ++it)
  {
    samVector samList = it->second;
    if (samList.size() >= 1)
    {
      for (samVector::iterator vecItr = samList.begin(); vecItr != samList.end(); vecItr++)
      {
        CSam *sam = *vecItr;
        mapReads[sam->readName]++;
      } // for end
    }
  } // for bed hash

  for (it = samHash.begin(); it != samHash.end(); ++it)
  {
    samVector samList = it->second;
    if (samList.size() >= 1)
    {
      for (samVector::iterator vecItr = samList.begin(); vecItr != samList.end(); vecItr++)
      {
        CSam *sam = *vecItr;
        int locusNum = mapReads[sam->readName];
        if (locusNum > paraInfo->maxLocusNum)
        {
          sam->readNum = 0;
          sam->strand = '#';
          maxReadNum++;
        }
        else {
          sam->readNum = sam->readNum / (double)locusNum;
          totalNum += sam->readNum;
        }
      } // for end
    }
  } // for bed hash

  fprintf(stderr, "remove %d reads with locus large than %d\n", maxReadNum, paraInfo->maxLocusNum);
  return totalNum;
}

void callPeaks(struct parameterInfo *paraInfo, chromSamMap &samHash,
               clusterVector &contigVector, FILE *gfp, faidxMap &faiHash)
{
  chromSamMap::iterator it;
  faidx *fai = NULL;
  for (it = samHash.begin(); it != samHash.end(); ++it)
  {
    char *chromName =  (char *)it->first.c_str();
    if (faiHash.find(it->first) == faiHash.end())
    {
      fprintf(stderr, "can't find the chromosome %s, skip it.\n", chromName);
      continue;
    }
    fai = faiHash[it->first];
    if (skipChrom(chromName))
    {
      fprintf(stderr, "skip the unusal chromosome %s.\n", chromName);
      continue;
    }
    samVector samList = it->second;
    if (samList.size() >= 1)
    {
      if (paraInfo->verbose) fprintf(stderr, "sort the chromosome: %s with %d reads\n", chromName, (int)samList.size());
      sort(samList.begin(), samList.end(), compareSam);
      if (paraInfo->verbose) fprintf(stderr, "get peak + strand\n");
      getPeakInfo(paraInfo, '+', samList, contigVector, gfp, fai); // plus strand
      if (paraInfo->verbose) fprintf(stderr, "get peak - strand\n");
      getPeakInfo(paraInfo, '-', samList, contigVector, gfp, fai); // minus strand
    }
  } // for bed hash
}

int getPeakInfo(struct parameterInfo *paraInfo, char strand, samVector &samList,
                clusterVector &contigVector, FILE *gfp, faidx *fai)
{
  CSam *subSam   = NULL;
  CSam *sam      = NULL;
  CSam *nextSam  = NULL;
  peakInfo *peak = NULL;
  int preStart   = 0;
  int preEnd     = 0;
  int tagNum     = 0;
  int peakNum    = 0;
  samVector subSamList;
  int peakStart = samList[0]->chromStart;
  int peakEnd   = samList[0]->chromEnd;
  if (paraInfo->verbose) fprintf(stderr, "get peaks from %d reads\n", (int)samList.size());
  for (samVector::iterator vecItr = samList.begin(); vecItr != samList.end(); vecItr++)
  {
    CSam *sam = *vecItr;
    if (sam->strand != strand) // skip other strand
    {
      continue;
    }
    if (overlapLength(peakStart, peakEnd, sam->chromStart, sam->chromEnd) > paraInfo->maxReadDist)
    {
      subSam = (CSam *)safeMalloc(sizeof(CSam));
      copySam(subSam, sam);
      subSamList.push_back(subSam);
      peakStart = MIN(peakStart, sam->chromStart);
      peakEnd = MAX(peakEnd, sam->chromEnd);
      tagNum++;
    } // if overlap
    else
    {
      if (tagNum >= 1)
      {
        getAllMutations(paraInfo, strand, peakStart, peakEnd, subSamList, contigVector, gfp, fai);
        if (paraInfo->verbose) fprintf(stderr, "%s %d %d %c\n", samList[0]->chrom, peakStart, peakEnd, strand);
        peakNum++;
      }
      freeSamVector(subSamList);
      subSam = (CSam *)safeMalloc(sizeof(CSam));
      copySam(subSam, sam);
      subSamList.push_back(subSam);
      tagNum    = 1;
      peakStart = sam->chromStart;
      peakEnd   = sam->chromEnd;
    } // else end for overlap
  } // for end
  // last time
  if (subSamList.size() > 0)
  {
    if (tagNum >= 1)
    {
      getAllMutations(paraInfo, strand, peakStart, peakEnd, subSamList, contigVector, gfp, fai);
      if (paraInfo->verbose) fprintf(stderr, "%s %d %d %c\n", samList[0]->chrom, peakStart, peakEnd, strand);
      peakNum++;
    }
    freeSamVector(subSamList);
  }
  return peakNum;
}


double **getProfile(struct parameterInfo *paraInfo, samVector &samList,
                    int start, int end, char strand, Cluster *pCluster)
{
  int i = 0;
  int j = 0;
  int index = 0;
  double **profile = NULL;
  double readNum = 0;
  int  readLen = 0;
  int span = end - start;
  profile  = alloc2array(span, ENCODENUM);
  for (samVector::iterator vecItr = samList.begin(); vecItr != samList.end(); vecItr++)
  {
    CSam *sam = *vecItr;
    readNum += sam->readNum;
    double **readProfile = NULL;

    readProfile = alloc2array(sam->chromEnd - sam->chromStart, ENCODENUM);
    getCIMS(readProfile, sam->readSeq, sam->cigar, sam->readNum, strand);
    //fprintf(stderr, "mdz:%s\t%s\t%c\n", sam->mdz, sam->readSeq, strand);
    for (i = sam->chromStart; i < sam->chromEnd; i++)
    {
      if (i >= start && i < end)
      {
        index = i - start;
        for (j = 0; j < ENCODENUM; j++)
          profile[index][j] = profile[index][j] + readProfile[i - sam->chromStart][j];
      }
    }

    if (strand == '+')
    {
      int startIdx = sam->chromStart - start;
      int endIdx   = sam->chromEnd - start - 1;
      if (startIdx >= 0 && startIdx < span)
        profile[startIdx][encodeAlignment('S', strand)] = profile[startIdx][encodeAlignment('S', strand)] + sam->readNum;
      if (endIdx >= 0 && endIdx < span)
        profile[endIdx][encodeAlignment('E', strand)] = profile[endIdx][encodeAlignment('E', strand)] + sam->readNum;
    }
    else
    {
      int startIdx = sam->chromEnd - start - 1;
      int endIdx   = sam->chromStart - start;
      if (startIdx >= 0 && startIdx < span)
        profile[startIdx][encodeAlignment('S', strand)] = profile[startIdx][encodeAlignment('S', strand)] + sam->readNum;
      if (endIdx >= 0 && endIdx < span)
        profile[endIdx][encodeAlignment('E', strand)] = profile[endIdx][encodeAlignment('E', strand)] + sam->readNum;
    }

    if (paraInfo->rmSeMutation)
    {
      removeSeMutations(profile, readProfile, sam, start, strand);
    }

    // free memory
    free2array(readProfile, sam->chromEnd - sam->chromStart);

  }
  pCluster->readNum = readNum;
  return profile;
}

void removeSeMutations(double **profile, double **readProfile, CSam *sam, int start, char strand)
{
  int i, j;
  int mdzLen = strlen(sam->mdz);
  int strMut = 0;
  int endMut = 0;
  int startIdx = 0;
  int endIdx = 0;
  int samStart = 0;
  int samEnd = 0;

  if (strand == '+')
  {
    if (sam->mdz[0] == '0' && isalpha(sam->mdz[1])) strMut = 1;
    if (sam->mdz[mdzLen - 1] == '0' && isalpha(sam->mdz[mdzLen - 2])) endMut = 1;
    startIdx = sam->chromStart - start;
    endIdx   = sam->chromEnd - start - 1;
    samStart = sam->chromStart - sam->chromStart;
    samEnd = sam->chromEnd - sam->chromStart - 1;
  }
  else
  {
    if (sam->mdz[0] == '0' && isalpha(sam->mdz[1])) endMut = 1;
    if (sam->mdz[mdzLen - 1] == '0' && isalpha(sam->mdz[mdzLen - 2])) strMut = 1;
    startIdx = sam->chromEnd - start - 1;
    endIdx = sam->chromStart - start;
    samStart = sam->chromEnd - sam->chromStart - 1;
    samEnd   = sam->chromStart - sam->chromStart;
  }
  if (strMut)
  {
    //fprintf(stderr, "startMutation: start=%d end=%d", startIdx, endIdx);
    profile[startIdx][encodeAlignment('S', strand)] -= sam->readNum;
    for (j = 0; j < ENCODENUM; j++)
      profile[startIdx][j] -= readProfile[samStart][j];
  }
  if (endMut)
  {
    //fprintf(stderr, "endMutation: start=%d end=%d", startIdx, endIdx);
    profile[endIdx][encodeAlignment('E', strand)] -= sam->readNum;
    for (j = 0; j < ENCODENUM; j++)
      profile[endIdx][j] -= readProfile[samEnd][j];
  }

}

void getAllMutations(struct parameterInfo *paraInfo, char strand, int start, int end, samVector &samList, clusterVector &contigVector, FILE *gfp, faidx *fai)
{
  int i = 0;
  int j = 0;
  int span = end - start;
  Cluster *pCluster = NULL;
  if (paraInfo->verbose) fprintf(stderr, "get all mutations: chrom:%s start:%d end:%d strand:%c\n", samList[0]->chrom, start, end, strand);
  char *chrom   = samList[0]->chrom;
  pCluster      = (Cluster *)safeMalloc(sizeof(Cluster));
  double **profile = getProfile(paraInfo, samList, start, end, strand, pCluster);
  char *clusterSeq   = faidxFetchSeq(gfp, fai, start, end, strand);
  identifyClusterInfo(clusterSeq, strand, 0, span, profile, pCluster, paraInfo);
  pCluster->profile  = profile;
  pCluster->chrom    = strClone(chrom);
  pCluster->start    = start;
  pCluster->end      = end;
  pCluster->strand   = strand;
  pCluster->seq      = strClone(clusterSeq);
  contigVector.push_back(pCluster);
  fflush(stderr);
  //safeFree(clusterSeq);
  //free2array(profile, end - start);
}

void identifyClusterInfo(char *clusterSeq, char strand, int start, int end, double **profile, Cluster *pCluster, struct parameterInfo *paraInfo)
{
  int i, j;
  double totalConversionNum         = 0;
  double totalConversionHeight      = 0;
  double totalMutationNum           = 0;
  double totalMutationHeight        = 0;
  double totalTruncationNum         = 0;
  double totalTruncationHeight      = 0;
  double totalDeletionNum           = 0;
  double totalDeletionHeight        = 0;
  double totalInsertionNum          = 0;
  double totalInsertionHeight       = 0;
  double totalEndNum                = 0;
  double totalEndHeight             = 0;
  double totalHeight                = 0;

  int span = strlen(clusterSeq);
  for (i = start; i < end; i++)
  {
    double height       = 0;
    int truncatedIdx    = -1;
    int mutationIdx     = -1;
    int deletionIdx     = -1;
    int conversionIdx   = -1;
    int insertionIdx    = -1;
    int endIdx          = -1;
    for (j = 0; j < ENCODENUM; j++)
    {
      if (strand == '-')
      {
        int m = span - i - 1;
        if (j != encodeAlignment('I', strand) && j != encodeAlignment('S', strand) && j != encodeAlignment('E', strand)) // discard insert reads
        {
          height += profile[m][j];
        }
        if (paraInfo->cvs != NULL)
        {
          if (clusterSeq[i] == paraInfo->cvs[0] && decodeAlignment(j) == paraInfo->cvs[1])
          {
            totalConversionNum += profile[m][j];
            if (profile[m][j] > 0) conversionIdx = i;
          }
        }
        // truncations
        if (decodeAlignment(j) == 'S')
        {
          totalTruncationNum += profile[m][j];
          if (profile[m][j] > 0) truncatedIdx = i;
        }
        // ends
        if (decodeAlignment(j) == 'E')
        {
          totalEndNum += profile[m][j];
          if (profile[m][j] > 0) endIdx = i;
        }
        // substitutions
        if (clusterSeq[i] != decodeAlignment(j) && j != encodeAlignment('D', strand) && j != encodeAlignment('I', strand) && j != encodeAlignment('S', strand) && j != encodeAlignment('E', strand))
        {
          totalMutationNum  += profile[m][j];
          if (profile[m][j] > 0)
          {
            mutationIdx = i;
          }
        }
        // deletions
        if (decodeAlignment(j) == 'D')
        {
          totalDeletionNum += profile[m][j];
          //totalMutationNum += profile[m][j];
          if (profile[m][j] > 0)
          {
            deletionIdx = i;
            //mutationIdx = i;
          }
        }
        // insertions
        if (decodeAlignment(j) == 'I')
        {
          totalInsertionNum += profile[m][j];
          //totalMutationNum += profile[m][j];
          if (profile[m][j] > 0)
          {
            insertionIdx = i;
            //mutationIdx = i;
          }
        }
      } // minus strand end
      else // plus strand start
      {
        if (j != encodeAlignment('I', strand) && j != encodeAlignment('S', strand) && j != encodeAlignment('E', strand)) // discard insert reads
        {
          height += profile[i][j];
        }
        if (paraInfo->cvs != NULL)
        {
          if (clusterSeq[i] == paraInfo->cvs[0] && decodeAlignment(j) == paraInfo->cvs[1])
          {
            totalConversionNum += profile[i][j];
            if (profile[i][j] > 0) conversionIdx = i;
          }
        }
        // truncations
        if (decodeAlignment(j) == 'S')
        {
          totalTruncationNum += profile[i][j];
          if (profile[i][j] > 0) truncatedIdx = i;
        }
        // ends
        if (decodeAlignment(j) == 'E')
        {
          totalEndNum += profile[i][j];
          if (profile[i][j] > 0) endIdx = i;
        }
        // substitutions
        if (clusterSeq[i] != decodeAlignment(j) && j != encodeAlignment('D', strand) && j != encodeAlignment('I', strand) && j != encodeAlignment('S', strand) && j != encodeAlignment('E', strand))
        {
          totalMutationNum  += profile[i][j];
          if (profile[i][j] > 0)
          {
            mutationIdx = i;
          }
        }
        // deletions
        if (decodeAlignment(j) == 'D')
        {
          totalDeletionNum += profile[i][j];
          //totalMutationNum += profile[i][j];
          if (profile[i][j] > 0)
          {
            //mutationIdx = i;
            deletionIdx = i;
          }
        }
        // insertions
        if (decodeAlignment(j) == 'I')
        {
          totalInsertionNum += profile[i][j];
          //totalMutationNum += profile[m][j];
          if (profile[i][j] > 0)
          {
            insertionIdx = i;
            //mutationIdx = i;
          }
        }
      } // else end
    } // ENCODENUM
    totalHeight += height;
    if (conversionIdx != -1)     totalConversionHeight   += height;
    if (truncatedIdx != -1)      totalTruncationHeight   += height;
    if (mutationIdx != -1)       totalMutationHeight     += height;
    if (deletionIdx != -1)       totalDeletionHeight     += height;
    if (insertionIdx != -1)      totalInsertionHeight    += height;
    if (endIdx != -1)            totalEndHeight          += height;
  } //span

  if (paraInfo->cvs != NULL ) pCluster->totalConversionNum    = totalConversionNum;
  if (paraInfo->cvs != NULL ) pCluster->totalConversionHeight = totalConversionHeight;
  pCluster->totalHeight               = totalHeight;
  pCluster->totalMutationNum          = totalMutationNum;
  pCluster->totalMutationHeight       = totalMutationHeight;
  pCluster->totalTruncationNum        = totalTruncationNum;
  pCluster->totalTruncationHeight     = totalTruncationHeight;
  pCluster->totalDeletionNum          = totalDeletionNum;
  pCluster->totalDeletionHeight       = totalDeletionHeight;
  pCluster->totalInsertionNum         = totalInsertionNum;
  pCluster->totalInsertionHeight      = totalInsertionHeight;
  pCluster->totalEndNum               = totalEndNum;
  pCluster->totalEndHeight            = totalEndHeight;

  if (paraInfo->cvs != NULL ) pCluster->t2cRatio = pCluster->totalConversionNum / pCluster->totalConversionHeight;
  pCluster->mutRatio = pCluster->totalMutationNum / pCluster->totalMutationHeight;
  pCluster->truRatio = pCluster->totalTruncationNum / pCluster->totalTruncationHeight;
  pCluster->delRatio = pCluster->totalDeletionNum / pCluster->totalDeletionHeight;
  pCluster->insRatio = pCluster->totalInsertionNum / pCluster->totalInsertionHeight;
  pCluster->endRatio = pCluster->totalEndNum / pCluster->totalEndHeight;
}

double getGenomeClusterInfo(struct parameterInfo *paraInfo, clusterVector &contigVector, genomeCluster *gCluster)
{
  double initVal        = 0.000001;
  gCluster->totalReadNum             = initVal;
  gCluster->totalConversionNum       = initVal;
  gCluster->totalConversionHeight    = initVal;
  gCluster->totalMutationNum         = initVal;
  gCluster->totalMutationHeight      = initVal;
  gCluster->totalTruncationNum       = initVal;
  gCluster->totalTruncationHeight    = initVal;
  gCluster->totalDeletionNum         = initVal;
  gCluster->totalDeletionHeight      = initVal;
  gCluster->totalInsertionNum        = initVal;
  gCluster->totalInsertionHeight     = initVal;
  gCluster->totalEndNum              = initVal;
  gCluster->totalEndHeight           = initVal;

  gCluster->t2cRatio = 0;
  gCluster->mutRatio = 0;
  gCluster->truRatio = 0;
  gCluster->delRatio = 0;
  gCluster->insRatio = 0;
  gCluster->endRatio = 0;
  gCluster->totalLen = 0;

  for (clusterVector::iterator vecItr = contigVector.begin(); vecItr != contigVector.end(); vecItr++)
  {
    Cluster *pCluster  = *vecItr;
    gCluster->totalReadNum            += pCluster->readNum;
    if (paraInfo->cvs != NULL ) gCluster->totalConversionNum     += pCluster->totalConversionNum;
    if (paraInfo->cvs != NULL ) gCluster->totalConversionHeight  += pCluster->totalConversionHeight;
    gCluster->totalHeight             += pCluster->totalHeight;
    gCluster->totalMutationNum        += pCluster->totalMutationNum;
    gCluster->totalMutationHeight     += pCluster->totalMutationHeight;
    gCluster->totalTruncationNum      += pCluster->totalTruncationNum;
    gCluster->totalTruncationHeight   += pCluster->totalTruncationHeight;
    gCluster->totalDeletionNum        += pCluster->totalDeletionNum;
    gCluster->totalDeletionHeight     += pCluster->totalDeletionHeight;
    gCluster->totalInsertionNum       += pCluster->totalInsertionNum;
    gCluster->totalInsertionHeight    += pCluster->totalInsertionHeight;
    gCluster->totalEndNum             += pCluster->totalEndNum;
    gCluster->totalEndHeight          += pCluster->totalEndHeight;
    gCluster->totalLen                += abs(pCluster->end - pCluster->start);
  }
  if (paraInfo->cvs != NULL ) gCluster->t2cRatio = gCluster->totalConversionNum / gCluster->totalConversionHeight;
  gCluster->mutRatio = gCluster->totalMutationNum / gCluster->totalMutationHeight;
  gCluster->truRatio = gCluster->totalTruncationNum / gCluster->totalTruncationHeight;
  gCluster->delRatio = gCluster->totalDeletionNum / gCluster->totalDeletionHeight;
  gCluster->insRatio = gCluster->totalInsertionNum / gCluster->totalInsertionHeight;
  gCluster->endRatio = gCluster->totalEndNum / gCluster->totalEndHeight;
  if (paraInfo->cvs != NULL ) fprintf(stderr, "the number of %c->%c mutations is %.0f, total number is %.0f, ratio is %.5f\n", paraInfo->cvs[0], paraInfo->cvs[1], gCluster->totalConversionNum, gCluster->totalConversionHeight, gCluster->t2cRatio);
  fprintf(stderr, "the number of mutations is %.0f, total number is %.0f, ratio is %.5f\n", gCluster->totalMutationNum, gCluster->totalMutationHeight, gCluster->mutRatio);
  fprintf(stderr, "the number of insertions is %.0f, total number is %.0f, ratio is %.5f\n", gCluster->totalInsertionNum, gCluster->totalInsertionHeight, gCluster->insRatio);
  fprintf(stderr, "the number of truncations is %.0f, total number is %.0f, ratio is %.5f\n", gCluster->totalTruncationNum, gCluster->totalTruncationHeight, gCluster->truRatio);
  fprintf(stderr, "the number of ends is %.0f, total number is %.0f, ratio is %.5f\n", gCluster->totalEndNum, gCluster->totalEndHeight, gCluster->endRatio);
  fprintf(stderr, "the number of deletions is %.0f, total number is %.0f, ratio is %.5f\n", gCluster->totalDeletionNum, gCluster->totalDeletionHeight, gCluster->delRatio);
  fprintf(stderr, "the total reads is %.5f\n", gCluster->totalReadNum);
  return gCluster->totalReadNum;
}

void outputSiteInfo(FILE *gfp, faidxMap &faiHash, char *outdir, struct parameterInfo *paraInfo, double totalReadNum, clusterVector &contigVector)
{
  siteVector conversionArray;
  siteVector mutationArray;
  siteVector truncationArray;
  siteVector peakArray;
  siteVector endArray;
  siteVector delArray;
  siteVector insArray;

  genomeCluster *gCluster = (genomeCluster *)safeMalloc(sizeof(genomeCluster));
  getGenomeClusterInfo(paraInfo, contigVector, gCluster);
  //double coverageRead = (gCluster->totalReadLen/gCluster->totalReadNum)*gCluster->totalReadNum / (double)(gCluster->totalLen);
  if (paraInfo->transcriptomeSize > gCluster->totalLen) gCluster->totalLen = paraInfo->transcriptomeSize;
  double coverageRead = gCluster->totalHeight / (double)gCluster->totalLen;
  fprintf(stderr, "the total height is %.5f, total length is %ld, mean coverage is %.5f\n", gCluster->totalHeight, gCluster->totalLen, coverageRead);
  for (clusterVector::iterator vecItr = contigVector.begin(); vecItr != contigVector.end(); vecItr++)
  {
    Cluster *pCluster  = *vecItr;
    identifyClusterSites(paraInfo, gfp, faiHash, coverageRead, pCluster, gCluster,
                         conversionArray, mutationArray, truncationArray,
                         peakArray, endArray, delArray, insArray);

    double s = pCluster->maxHeight;
    double lambda = coverageRead;
    //double localLambda = pCluster->totalHeight / (double)(pCluster->end - pCluster->start);
    //if (localLambda < lambda) lambda = localLambda;
    double mfold = s / lambda;
    double logPoissonP = ilogCumulativePoisson(s, lambda);
    pCluster->pval = log10Lnp(logPoissonP);
    pCluster->mfold = mfold;
    pCluster->heightRpm = pCluster->maxHeight / gCluster->totalReadNum * 1000000;
  }
  if (conversionArray.size() > 0)
    outputTypeSites(paraInfo, outdir, gfp, faiHash, conversionArray);
  if (mutationArray.size() > 0)
    outputTypeSites(paraInfo, outdir, gfp, faiHash, mutationArray);
  if (truncationArray.size() > 0)
    outputTypeSites(paraInfo, outdir, gfp, faiHash, truncationArray);
  if (endArray.size() > 0)
    outputTypeSites(paraInfo, outdir, gfp, faiHash, endArray);
  if (delArray.size() > 0)
    outputTypeSites(paraInfo, outdir, gfp, faiHash, delArray);
  if (insArray.size() > 0)
    outputTypeSites(paraInfo, outdir, gfp, faiHash, insArray);
  if (peakArray.size() > 0)
    outputPeakHeightTypeSites(paraInfo, outdir, gfp, faiHash, peakArray);
  if (contigVector.size() > 0)
    outputPeakSites(paraInfo, outdir, gfp, faiHash, contigVector);

  outputTransMatrix(paraInfo, outdir);

  safeFree(gCluster);
  freeSiteVector(conversionArray);
  freeSiteVector(mutationArray);
  freeSiteVector(truncationArray);
  freeSiteVector(peakArray);
  freeSiteVector(endArray);
  freeSiteVector(delArray);
  freeSiteVector(insArray);
}

void outputTransMatrix(struct parameterInfo *paraInfo, char *outdir)
{

  int i, j;
  double pseudoCount = 0.0005;
  double mutReadNum = pseudoCount;
  double mutSiteNum = pseudoCount;
  double delReadNum = pseudoCount;
  double delSiteNum = pseudoCount;
  double insReadNum = pseudoCount;
  double insSiteNum = pseudoCount;
  double strReadNum = pseudoCount;
  double strSiteNum = pseudoCount;
  double endReadNum = pseudoCount;
  double endSiteNum = pseudoCount;
  const char *mutationType[] = {"Deletion", "Insertion", "Truncation", "End", "Mutation", "Unknown"};
  const char *varType = mutationType[5];
  char outfile[512];
  strcpy(outfile, outdir);
  strcat(outfile, "_");
  strcat(outfile, "transMatrix");
  strcat(outfile, ".txt");
  FILE *outfp = NULL;
  outfp = (FILE *) fopen(outfile, "w");
  if (outfp == NULL)
  {
    fprintf(stderr, "ERROR: Can't open %s\n", outfile);
    exit(1);
  }

  for (i = 0; i < ENCODENUM; i++)
  {
    for (j = 0; j < ENCODENUM; j++)
    {
      char c = decodeAlignment(j);
      double siteNum = transMatrixNum[i][j];
      double readNum = transMatrix[i][j];
      switch (c)
      {
      case 'D':
        delSiteNum += siteNum;
        delReadNum += readNum;
        break;
      case 'I':
        insSiteNum += siteNum;
        insReadNum += readNum;
        break;
      case 'S':
        strSiteNum += siteNum;
        strReadNum += readNum;
        break;
      case 'E':
        endSiteNum += siteNum;
        endReadNum += readNum;
        break;
      case 'A':
      case 'C':
      case 'G':
      case 'T':
      case 'N':
        mutSiteNum += siteNum;
        mutReadNum += readNum;
        break;
      default:
        fprintf(stderr, "unknow base %c\n", c);
        break;
      }
    }
  }

  fprintf(outfp, "refBase\tvarBase\tvarType\tvarSiteNum\tvarReadNum\tvarSitePercent\tvarReadPercent\n");
  for (i = 0; i < ENCODENUM; i++)
  {
    for (j = 0; j < ENCODENUM; j++)
    {
      char c = decodeAlignment(j);
      double siteNum = transMatrixNum[i][j];
      double readNum = transMatrix[i][j];
      double siteRatio = 0;
      double readRatio = 0;
      switch (c)
      {
      case 'D':
        siteRatio = siteNum / delSiteNum;
        readRatio = readNum / delReadNum;
        varType = mutationType[0];
        break;
      case 'I':
        siteRatio = siteNum / insSiteNum;
        readRatio = readNum / insReadNum;
        varType = mutationType[1];
        break;
      case 'S':
        siteRatio = siteNum / strSiteNum;
        readRatio = readNum / strReadNum;
        varType = mutationType[2];
        break;
      case 'E':
        siteRatio = siteNum / endSiteNum;
        readRatio = readNum / endReadNum;
        varType = mutationType[3];
        break;
      case 'A':
      case 'C':
      case 'G':
      case 'T':
      case 'N':
        siteRatio = siteNum / mutSiteNum;
        readRatio = readNum / mutReadNum;
        varType = mutationType[4];
        break;
      default:
        fprintf(stderr, "unknow base %c\n", c);
        break;
      }
      if (readNum > 0)
        fprintf(outfp, "%c\t%c\t%s\t%.0f\t%.0f\t%.5f\t%.5f\n",
                decodeAlignment(i), decodeAlignment(j),  varType, siteNum, readNum, siteRatio * 100, readRatio * 100);
    }
  }
}

int outputPeakSites(struct parameterInfo * paraInfo, char *outdir, FILE * gfp, faidxMap & faiHash, clusterVector & contigVector)
{
  int extendLen = EXTEND_LEN;
  int siteIdx = 1;
  char outfile[512];
  strcpy(outfile, outdir);
  strcat(outfile, "_");
  strcat(outfile, "Peak");
  strcat(outfile, ".bed");
  FILE *outfp = NULL;
  int i = 0;
  outfp = (FILE *) fopen(outfile, "w");
  if (outfp == NULL)
  {
    fprintf(stderr, "ERROR: Can't open %s\n", outfile);
    exit(1);
  }
  getClusterBhFdr(contigVector);
  fprintf(outfp, "#chrom\tchromStart\tchromEnd\tname\tscore\tstrand\t");
  fprintf(outfp, "extendSeq\tmotifPos\ttype\tlog10(p-value)\tlog10(q-value)\treadNum\theight\theightRpm\tmfold\tratio\n");
  for (clusterVector::iterator vecItr = contigVector.begin(); vecItr != contigVector.end(); vecItr++)
  {
    Cluster *pCluster = *vecItr;
    if (pCluster->pval < paraInfo->pval
        && pCluster->qval < paraInfo->qval
        && pCluster->mfold >= paraInfo->mfold
        && pCluster->maxHeight >= paraInfo->minHeight
        && pCluster->heightRpm >= paraInfo->minRpm)
    {
      char *chrom       = pCluster->chrom;
      int chromStart    = pCluster->start + pCluster->maxHeightPos;
      int chromEnd      = chromStart + 1;
      if (paraInfo->fullLength)
      {
        chromStart = pCluster->start;
        chromEnd   = pCluster->end;
        extendLen  = 0;
      }
      char strand = pCluster->strand;
      faidx *fai  = NULL;
      string chromStr(chrom);
      fai = faiHash[chromStr];
      int extendStart = chromStart - extendLen;
      int extendEnd = chromEnd + extendLen;
      if (extendStart < 0) extendStart = 0;
      if (extendEnd > fai->len) extendEnd = fai->len;
      char *extendSeq   = (char *)faidxFetchSeq(gfp, fai, extendStart, extendEnd, strand);
      int motifPos = searchMotif(extendSeq);
      pCluster->ratio = pCluster->maxHeight / pCluster->readNum;
      fprintf(outfp, "%s\t%d\t%d\trbsSeeker_peak_%d\t%.0f\t%c\t%s\t%d\t",
              chrom, chromStart, chromEnd, siteIdx, pCluster->maxHeight, strand, extendSeq, motifPos);
      fprintf(outfp, "peak\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\n",
              pCluster->pval, pCluster->qval, pCluster->readNum, pCluster->maxHeight, pCluster->maxHeight, pCluster->heightRpm, pCluster->ratio);
      siteIdx++;
      safeFree(extendSeq);
    }
  }
  fclose(outfp);
  fprintf(stderr, "get %d\tpeaks\tfrom %lu sites\n", siteIdx, contigVector.size());
  return siteIdx;
}

int outputPeakHeightTypeSites(struct parameterInfo * paraInfo, char *outdir, FILE * gfp, faidxMap & faiHash, siteVector & siteArray)
{
  int extendLen = EXTEND_LEN;
  int siteIdx = 1;
  char outfile[512];
  strcpy(outfile, outdir);
  strcat(outfile, "_");
  strcat(outfile, siteArray[0]->type);
  strcat(outfile, ".bed");
  FILE *outfp = NULL;
  int i = 0;
  outfp = (FILE *) fopen(outfile, "w");
  if (outfp == NULL)
  {
    fprintf(stderr, "ERROR: Can't open %s\n", outfile);
    exit(1);
  }
  getSiteBhFdr(siteArray);
  sort(siteArray.begin(), siteArray.end(), cmpChromLocus);
  fprintf(outfp, "#chrom\tchromStart\tchromEnd\tname\tscore\tstrand\t");
  fprintf(outfp, "extendSeq\tmotifPos\ttype\tlog10(p-value)\tlog10(q-value)\treadNum\theight\theightRpm\tmfold\tratio\n");
  fflush(outfp);
  siteIdx = outputStrandRbsSite(paraInfo, outfp, gfp, faiHash, siteArray, siteIdx, '+');
  siteIdx = outputStrandRbsSite(paraInfo, outfp, gfp, faiHash, siteArray, siteIdx, '-');
  fclose(outfp);
  fprintf(stderr, "get %d\t%ss\tfrom %lu sites\n", siteIdx, siteArray[0]->type, siteArray.size());
  return siteIdx;
}

int outputStrandRbsSite(struct parameterInfo * paraInfo, FILE * outfp, FILE * gfp, faidxMap & faiHash, siteVector & siteArray, int siteIdx, char strand)
{
  // output + strand
  RBSite *rbsPeak = NULL;
  int firstTag = 0;
  int maxDist = -50; // the maximum distance between two sites
  for (siteVector::iterator vecItr = siteArray.begin(); vecItr != siteArray.end(); vecItr++)
  {
    RBSite *rbs = *vecItr;
    if (rbs->strand != strand) continue;
    if (rbs->pval < paraInfo->pval && rbs->qval < paraInfo->qval && rbs->mfold >= paraInfo->mfold)
    {
      if (firstTag == 0)
      {
        firstTag = 1;
        rbsPeak = (RBSite *)safeMalloc(sizeof(RBSite));
        copyRBSite(rbsPeak, rbs);
        continue;
      }
      if (strcmp(rbsPeak->chrom, rbs->chrom) == 0 && overlapLength(rbsPeak->start, rbsPeak->end, rbs->start, rbs->end) >= maxDist)
      {
        rbsPeak->start = MIN(rbsPeak->start, rbs->start);
        rbsPeak->end   = MAX(rbsPeak->end, rbs->end);
        rbsPeak->readNum  = MAX(rbsPeak->readNum, rbs->readNum);
        rbsPeak->height   = MAX(rbsPeak->height, rbs->height);
        rbsPeak->ratio    = MAX(rbsPeak->ratio, rbs->ratio);
        rbsPeak->mfold    = MAX(rbsPeak->mfold, rbs->mfold);
        rbsPeak->pval     = MIN(rbsPeak->pval, rbs->pval);
        rbsPeak->qval     = MIN(rbsPeak->qval, rbs->qval);
        rbsPeak->heightRpm = MAX(rbsPeak->heightRpm, rbs->heightRpm);
      }
      else {
        siteIdx = outputRbsSite(paraInfo, outfp, gfp, faiHash, rbsPeak, siteIdx);
        freeSite(rbsPeak);
        rbsPeak = NULL;
        rbsPeak = (RBSite *)safeMalloc(sizeof(RBSite));
        copyRBSite(rbsPeak, rbs);
      }
    }
  }// for site vectors
  if (rbsPeak != NULL)
  {
    siteIdx = outputRbsSite(paraInfo, outfp, gfp, faiHash, rbsPeak, siteIdx);
    freeSite(rbsPeak);
  }
  return siteIdx;
}

int outputRbsSite(struct parameterInfo * paraInfo, FILE * outfp, FILE * gfp, faidxMap & faiHash, RBSite * rbs, int siteIdx)
{
  char *chrom     = rbs->chrom;
  int chromStart  = rbs->start;
  int chromEnd    = rbs->end;
  char strand     = rbs->strand;
  faidx *fai = NULL;
  string chromStr(chrom);
  fai = faiHash[chromStr];
  int extendStart = chromStart;
  int extendEnd = chromEnd;
  if (extendStart < 0) extendStart = 0;
  if (extendEnd > fai->len) extendEnd = fai->len;
  char *extendSeq   = (char *)faidxFetchSeq(gfp, fai, extendStart, extendEnd, strand);
  int motifPos = searchMotif(extendSeq);
  if (rbs->pval < paraInfo->pval
      && rbs->qval < paraInfo->qval
      && rbs->mfold >= paraInfo->mfold
      && rbs->heightRpm >= paraInfo->minRpm)
  {
    fprintf(outfp, "%s\t%d\t%d\trbsSeeker_%s_%d\t%.0f\t%c\t%s\t%d\t", chrom, chromStart, chromEnd, rbs->type, siteIdx, rbs->height, strand, extendSeq, motifPos);
    fprintf(outfp, "%s\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\n", rbs->type, rbs->pval, rbs->qval, rbs->readNum, rbs->height, rbs->heightRpm, rbs->mfold, rbs->ratio);
  }
  safeFree(extendSeq);
  siteIdx++;
  fflush(outfp);
  return siteIdx;
}

void copyRBSite(RBSite * tRbs, RBSite * oRbs)
{
  tRbs->chrom = strClone(oRbs->chrom);
  tRbs->start = oRbs->start;
  tRbs->end   = oRbs->end;
  tRbs->strand = oRbs->strand;
  tRbs->readNum = oRbs->readNum;
  tRbs->height = oRbs->height;
  tRbs->heightRpm = oRbs->heightRpm;
  tRbs->ratio = oRbs->ratio;
  tRbs->mfold = oRbs->mfold;
  tRbs->pval = oRbs->pval;
  tRbs->qval = oRbs->qval;
  tRbs->type = strClone(oRbs->type);
}

int outputTypeSites(struct parameterInfo * paraInfo, char *outdir, FILE * gfp, faidxMap & faiHash, siteVector & siteArray)
{
  int extendLen = EXTEND_LEN;
  int siteIdx = 0;
  char outfile[512];
  strcpy(outfile, outdir);
  strcat(outfile, "_");
  strcat(outfile, siteArray[0]->type);
  strcat(outfile, ".bed");
  FILE *outfp = NULL;
  int i = 0;
  outfp = (FILE *) fopen(outfile, "w");
  if (outfp == NULL)
  {
    fprintf(stderr, "ERROR: Can't open %s\n", outfile);
    exit(1);
  }
  getSiteBhFdr(siteArray);
  fprintf(outfp, "#chrom\tchromStart\tchromEnd\tname\tscore\tstrand\t");
  fprintf(outfp, "extendSeq\tmotifPos\ttype\tlog10(p-value)\tlog10(q-value)\treadNum\theight\theightRpm\tmfold\tratio\n");
  for (siteVector::iterator vecItr = siteArray.begin(); vecItr != siteArray.end(); vecItr++)
  {
    RBSite *rbs = *vecItr;
    if (rbs->mfold >= paraInfo->mfold
        && rbs->heightRpm >= paraInfo->minRpm
        && rbs->readNum >= paraInfo->minMutNum
        && rbs->pval < paraInfo->pval
        && rbs->qval < paraInfo->qval
        && rbs->ratio >= paraInfo->minRatio
        && rbs->ratio <= paraInfo->maxRatio)
    {
      char *chrom       = rbs->chrom;
      int chromStart    = rbs->start;
      int chromEnd      = rbs->end;
      char strand       = rbs->strand;
      faidx *fai = NULL;
      string chromStr(chrom);
      fai = faiHash[chromStr];
      int extendStart = chromStart - extendLen;
      int extendEnd = chromEnd + extendLen;
      if (extendStart < 0) extendStart = 0;
      if (extendEnd > fai->len) extendEnd = fai->len;
      char *extendSeq   = (char *)faidxFetchSeq(gfp, fai, extendStart, extendEnd, strand);
      int motifPos = searchMotif(extendSeq);
      fprintf(outfp, "%s\t%d\t%d\trbsSeeker_%s_%d\t%.0f\t%c\t%s\t%d\t", chrom, chromStart, chromEnd, rbs->type, siteIdx, rbs->height, strand, extendSeq, motifPos);
      fprintf(outfp, "%s\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\n", rbs->type, rbs->pval, rbs->qval, rbs->readNum, rbs->height, rbs->heightRpm, rbs->mfold, rbs->ratio);
      siteIdx++;
      safeFree(extendSeq);
    }
  }
  fclose(outfp);
  fprintf(stderr, "get significant %d\t%ss\tfrom %lu sites\n", siteIdx, siteArray[0]->type, siteArray.size());
  return siteIdx;
}

/*statistic transmatrix */
double transMatrix[ENCODENUM][ENCODENUM];
double transMatrixNum[ENCODENUM][ENCODENUM];

void initTransMatrix(void)
{
  int i, j;
  for (i = 0; i < ENCODENUM; i++)
  {
    for (j = 0; j < ENCODENUM; j++)
    {
      transMatrix[i][j] = 0;
      transMatrixNum[i][j] = 0;
    }
  }
}

int identifyClusterSites(struct parameterInfo * paraInfo, FILE * gfp,
                         faidxMap & faiHash, double coverageRead,
                         Cluster * pCluster, genomeCluster * gCluster,
                         siteVector & conversionArray, siteVector & mutationArray,
                         siteVector & truncationArray, siteVector & peakArray,
                         siteVector & endArray, siteVector & delArray, siteVector & insArray)
{
  int i, j;
  char *clusterSeq  = pCluster->seq;
  char *chrom       = pCluster->chrom;
  int chromStart    = pCluster->start;
  int chromEnd      = pCluster->end;
  char strand       = pCluster->strand;
  double **profile  = pCluster->profile;
  int span          = strlen(clusterSeq);
  const char *mutationType[] = {"TC", "Mutation", "Truncation", "Deletion", "PeakHeight", "End", "Insertion"};
  int start = 0;
  int end   = span;
  int pos   = start;
  double maxHeight = 0;
  double maxHeightPos = 0;
  double minMutNum = 1.0;
  for (i = start; i < end; i++)
  {
    double height               = 0;
    double totalConversionNum   = 0;
    double totalTruncationNum   = 0;
    double totalMutationNum     = 0;
    double totalDeletionNum     = 0;
    double totalInsertionNum    = 0;
    double totalEndNum          = 0;
    int idx = encodeAlignment(clusterSeq[i], '+');

    int findSiteTag = 0;
    for (j = 0; j < ENCODENUM; j++)
    {
      if (strand == '-')
      {
        int m = span - i - 1;
        pos = m;
        if (j != encodeAlignment('I', strand) && j != encodeAlignment('S', strand) && j != encodeAlignment('E', strand)) // discard insert reads
        {
          height += profile[m][j];
        }
        if (paraInfo->cvs != NULL)
        {
          if (clusterSeq[i] == paraInfo->cvs[0] && decodeAlignment(j) == paraInfo->cvs[1])
          {
            totalConversionNum += profile[m][j];
          }
        }
        // truncations
        if (decodeAlignment(j) == 'S')
        {
          totalTruncationNum += profile[m][j];
        }
        // ends
        if (decodeAlignment(j) == 'E')
        {
          totalEndNum += profile[m][j];
        }
        // substitutions
        if (clusterSeq[i] != decodeAlignment(j) && j != encodeAlignment('D', strand) && j != encodeAlignment('I', strand) && j != encodeAlignment('S', strand) && j != encodeAlignment('E', strand))
        {
          totalMutationNum  += profile[m][j];
        }
        // deletions
        if (decodeAlignment(j) == 'D')
        {
          totalDeletionNum += profile[m][j];
          //totalMutationNum += profile[m][j];
        }
        // insertions
        if (decodeAlignment(j) == 'I')
        {
          totalInsertionNum += profile[m][j];
        }
        // transmatrix
        if (clusterSeq[i] != decodeAlignment(j))
        {
          transMatrix[idx][j]  += profile[m][j];
          if (profile[m][j] > 0) transMatrixNum[idx][j] += 1;
        }
      } // minus strand end
      else // plus strand
      {
        pos = i;
        if (j != encodeAlignment('I', strand) && j != encodeAlignment('S', strand) && j != encodeAlignment('E', strand)) // discard insert reads
        {
          height += profile[i][j];
        }
        if (paraInfo->cvs != NULL)
        {
          if (clusterSeq[i] == paraInfo->cvs[0] && decodeAlignment(j) == paraInfo->cvs[1])
          {
            totalConversionNum += profile[i][j];
          }
        }
        // truncations
        if (decodeAlignment(j) == 'S')
        {
          totalTruncationNum += profile[i][j];
        }
        // ends
        if (decodeAlignment(j) == 'E')
        {
          totalEndNum += profile[i][j];
        }
        // substitutions
        if (clusterSeq[i] != decodeAlignment(j) && j != encodeAlignment('D', strand) && j != encodeAlignment('I', strand) && j != encodeAlignment('S', strand) && j != encodeAlignment('E', strand))
        {
          totalMutationNum += profile[i][j];
        }
        // deletions
        if (decodeAlignment(j) == 'D')
        {
          totalDeletionNum += profile[i][j];
          //totalMutationNum += profile[i][j];
        }
        // insertions
        if (decodeAlignment(j) == 'I')
        {
          totalInsertionNum += profile[i][j];
        }
        //
        if (clusterSeq[i] != decodeAlignment(j))
        {
          transMatrix[idx][j]  += profile[i][j];
          if (profile[i][j] > 0) transMatrixNum[idx][j] += 1;
        }
      } // else end
    } // ENCODENUM
    if (height > maxHeight)
    {
      maxHeight = height;
      maxHeightPos = pos;
    }
    double heightRpm = height / gCluster->totalReadNum * 1000000;
    if (height >= paraInfo->minHeight)
    {
      double lambda = coverageRead;
      double mfold = height / lambda;
      double ratio = height / pCluster->readNum;
      double logPoissonP = ilogCumulativePoisson(height, lambda);
      double pval = log10Lnp(logPoissonP);
      if (mfold > 0)
      {
        RBSite *rbs   = (RBSite *)safeMalloc(sizeof(RBSite));
        rbs->readNum  = pCluster->readNum;
        rbs->height   = height;
        rbs->mfold    = mfold;
        rbs->ratio    = ratio;
        rbs->pval     = pval;
        rbs->chrom    = strClone(chrom);
        rbs->start    = chromStart + pos;
        rbs->end      = rbs->start + 1;
        rbs->strand   = strand;
        rbs->heightRpm = heightRpm;
        rbs->type     = strClone((char *)mutationType[4]);
        peakArray.push_back(rbs);
        findSiteTag   = 1;
      }
    } // peak hight
    if (paraInfo->cvs != NULL && totalConversionNum >= minMutNum && height >= paraInfo->minHeight)
    {
      double pval = loghypergeoD(gCluster->totalConversionHeight, gCluster->totalConversionNum, height, totalConversionNum);
      double ratio = totalConversionNum / height;
      double mfold = ratio / gCluster->t2cRatio;
      pval = log10Lnp(pval);
      if (mfold > 0)
      {
        RBSite *rbs   = (RBSite *)safeMalloc(sizeof(RBSite));
        rbs->readNum  = totalConversionNum;
        rbs->height   = height;
        rbs->mfold    = mfold;
        rbs->ratio    = ratio;
        rbs->pval     = pval;
        rbs->chrom    = strClone(chrom);
        rbs->start    = chromStart + pos;
        rbs->end      = rbs->start + 1;
        rbs->strand   = strand;
        rbs->heightRpm = heightRpm;
        rbs->type     = strClone((char *)paraInfo->cvs);
        conversionArray.push_back(rbs);
        findSiteTag   = 1;
      }
    } // conversion
    // truncation
    if (totalTruncationNum >= minMutNum && height >= paraInfo->minHeight)
    {
      double pval = loghypergeoD(gCluster->totalTruncationHeight, gCluster->totalTruncationNum, height, totalTruncationNum);
      double ratio = totalTruncationNum / height;
      double mfold = ratio / gCluster->truRatio;
      pval = log10Lnp(pval);
      if (mfold > 0)
      {
        RBSite *rbs   = (RBSite *)safeMalloc(sizeof(RBSite));
        rbs->readNum  = totalTruncationNum;
        rbs->height   = height;
        rbs->mfold    = mfold;
        rbs->ratio    = ratio;
        rbs->pval     = pval;
        rbs->chrom    = strClone(chrom);
        rbs->start    = chromStart + pos;
        rbs->end      = rbs->start + 1;
        rbs->strand   = strand;
        rbs->heightRpm = heightRpm;
        rbs->type     = strClone((char *)mutationType[2]);
        truncationArray.push_back(rbs);
        findSiteTag   = 1;
      }
    }// truncation
    //end
    if (totalEndNum >= minMutNum && height >= paraInfo->minHeight)
    {
      double pval = loghypergeoD(gCluster->totalEndHeight, gCluster->totalEndNum, height, totalEndNum);
      double ratio = totalEndNum / height;
      double mfold = ratio / gCluster->endRatio;
      pval = log10Lnp(pval);
      if (mfold > 0)
      {
        RBSite *rbs   = (RBSite *)safeMalloc(sizeof(RBSite));
        rbs->readNum  = totalEndNum;
        rbs->height   = height;
        rbs->mfold    = mfold;
        rbs->ratio    = ratio;
        rbs->pval     = pval;
        rbs->chrom    = strClone(chrom);
        rbs->start    = chromStart + pos;
        rbs->end      = rbs->start + 1;
        rbs->strand   = strand;
        rbs->heightRpm = heightRpm;
        rbs->type     = strClone((char *)mutationType[5]);
        endArray.push_back(rbs);
        findSiteTag   = 1;
      }
    }// ends
    // mutations
    if (totalMutationNum >= minMutNum && height >= paraInfo->minHeight)
    {
      double pval = loghypergeoD(gCluster->totalMutationHeight, gCluster->totalMutationNum, height, totalMutationNum);
      double ratio = totalMutationNum / height;
      double mfold = ratio / gCluster->mutRatio;
      pval = log10Lnp(pval);
      if (mfold > 0)
      {
        RBSite *rbs   = (RBSite *)safeMalloc(sizeof(RBSite));
        rbs->readNum  = totalMutationNum;
        rbs->height   = height;
        rbs->mfold    = mfold;
        rbs->ratio    = ratio;
        rbs->pval     = pval;
        rbs->chrom    = strClone(chrom);
        rbs->start    = chromStart + pos;
        rbs->end      = rbs->start + 1;
        rbs->strand   = strand;
        rbs->heightRpm = heightRpm;
        rbs->type     = strClone((char *)mutationType[1]);
        mutationArray.push_back(rbs);
        findSiteTag   = 1;
      }
    }// mutation
    // totalDeletionNum
    if (totalDeletionNum >= minMutNum && height >= paraInfo->minHeight)
    {
      double pval = loghypergeoD(gCluster->totalDeletionHeight, gCluster->totalDeletionNum, height, totalDeletionNum);
      double ratio = totalDeletionNum / height;
      double mfold = ratio / gCluster->delRatio;
      pval = log10Lnp(pval);
      if (mfold > 0)
      {
        RBSite *rbs   = (RBSite *)safeMalloc(sizeof(RBSite));
        rbs->readNum  = totalDeletionNum;
        rbs->height   = height;
        rbs->mfold    = mfold;
        rbs->ratio    = ratio;
        rbs->pval     = pval;
        rbs->chrom    = strClone(chrom);
        rbs->start    = chromStart + pos;
        rbs->end      = rbs->start + 1;
        rbs->strand   = strand;
        rbs->heightRpm = heightRpm;
        rbs->type     = strClone((char *)mutationType[3]);
        delArray.push_back(rbs);
        findSiteTag   = 1;
      }
    }// deletion
    // totalInsertionNum
    if (totalInsertionNum >= minMutNum && height >= paraInfo->minHeight)
    {
      double pval = loghypergeoD(gCluster->totalInsertionHeight, gCluster->totalInsertionNum, height, totalInsertionNum);
      double ratio = totalInsertionNum / height;
      double mfold = ratio / gCluster->insRatio;
      pval = log10Lnp(pval);
      if (mfold > 0)
      {
        RBSite *rbs    = (RBSite *)safeMalloc(sizeof(RBSite));
        rbs->readNum   = totalInsertionNum;
        rbs->height    = height;
        rbs->mfold     = mfold;
        rbs->ratio     = ratio;
        rbs->pval      = pval;
        rbs->chrom     = strClone(chrom);
        rbs->start     = chromStart + pos;
        rbs->end       = rbs->start + 1;
        rbs->strand    = strand;
        rbs->heightRpm = heightRpm;
        rbs->type      = strClone((char *)mutationType[6]);
        insArray.push_back(rbs);
        findSiteTag   = 1;
      }
    }// deletion
  } //span
  pCluster->maxHeight = maxHeight;
  pCluster->maxHeightPos = maxHeightPos;
  return 1;
}

int searchMotif(char *seq)
{
  int pos = -1;
  int seqLen = strlen(seq);
  int i = 0;
  int motifLen = 4;
  int midLen = EXTEND_LEN;
  for (i = 0; i < seqLen - motifLen; i++)
  {
    if ((seq[i] == 'G' || seq[i] == 'A' || seq[i] == 'T')
        && (seq[i + 1] == 'G' || seq[i + 1] == 'A')
        && (seq[i + 2] == 'A')
        && (seq[i + 3] == 'C')
       )
    {
      pos = i;
      if (abs(pos - midLen) < 5) // next to variation site
        break;
    }
  }
  return pos;
}

void freeCluster(Cluster * clust)
{
  safeFree(clust->chrom);
  safeFree(clust->seq);
  int profileLen = clust->end - clust->start;
  free2array(clust->profile, profileLen);
  safeFree(clust);
}

void freeClusterVector(clusterVector & cList)
{
  for (clusterVector::iterator vecItr = cList.begin(); vecItr != cList.end(); vecItr++)
  {
    Cluster *clust = *vecItr;
    freeCluster(clust);
  }
  cList.clear();
}

void freeSite(RBSite * rbs)
{
  safeFree(rbs->chrom);
  safeFree(rbs->type);
  safeFree(rbs);
}

void freeSiteVector(siteVector & sList)
{
  for (siteVector::iterator vecItr = sList.begin(); vecItr != sList.end(); vecItr++)
  {
    RBSite *rbs = *vecItr;
    freeSite(rbs);
  }
  sList.clear();
}

void freeParameters(parameterInfo * paraInfo)
{
  if (paraInfo->cvs != NULL)
    safeFree(paraInfo->cvs);
  safeFree(paraInfo);
}

void getCIMS(double **readProfile, char *seq, char *cigar, double readNum, char strand)
{
  const char *samcigar = "MIDNSHP=X";
  int start = 0;
  int end   = 0;
  int i = 0;
  int j = 0;
  int readLen  = 0;
  int alignLen = 0;
  int profileLen = 0;
  for (i = 0; i < strlen(cigar); i++)
  {
    char c = cigar[i];
    if (!isdigit(c))
    {
      cigar[i] = '\0';
      int len = atoi(cigar + start);
      switch ( c )
      {
      //increase end position on CIGAR chars [DMXN=]
      case 'D' :
        for (j = 0; j < len; j++)
        {
          readProfile[profileLen + j][encodeAlignment('D', strand)] = readNum;
        }
        alignLen += len;
        profileLen += len;
        break;
      case 'I' :
        readProfile[profileLen][encodeAlignment('I', strand)] = readNum;
        alignLen += len;
        readLen  += len;
        break;
      case 'M' :
      case '=' :
        for (j = 0; j < len; j++)
        {
          readProfile[profileLen + j][encodeAlignment(seq[readLen + j], strand)] = readNum;
        }
        alignLen += len;
        readLen  += len;
        profileLen += len;
        break;
      }
      cigar[i] = c;
      start = i + 1;
    }
  }
}

char decodeAlignment(int ch)
{
  if (ch == 0)
    return 'A';
  else if (ch == 1)
    return 'C';
  else if (ch == 2)
    return 'G';
  else if (ch == 3)
    return 'T';
  else if (ch == 4)
    return 'D';
  else if (ch == 5)
    return 'I';
  else if (ch == 6)
    return 'S';
  else if (ch == 7)
    return 'E';
  else
    return 'N';
}

int encodeAlignment(char ch, char strand)
{
  ch = toupper(ch);
  if (strand == '-')
  {
    if (ch == 'A')
      ch = 'T';
    else if (ch == 'C')
      ch = 'G';
    else if (ch == 'G')
      ch = 'C';
    else if (ch == 'T' || ch == 'U')
      ch = 'A';
  }
  if (ch == 'A')
    return 0;
  else if (ch == 'C')
    return 1;
  else if (ch == 'G')
    return 2;
  else if (ch == 'T' || ch == 'U')
    return 3;
  else if (ch == 'D')
    return 4;
  else if (ch == 'I')
    return 5;
  else if (ch == 'S')
    return 6;
  else if (ch == 'E')
    return 7;
  else
    return 8;
}

double **alloc2array(int nrows, int ncolumns)
{
  int i = 0;
  int j = 0;
  double **array = NULL;
  array = (double **)safeMalloc(nrows * sizeof(double *));
  if (array == NULL)
  {
    fprintf(stderr, "out of memory\n");
    exit(1);
  }
  for (i = 0; i < nrows; i++)
  {
    array[i] = (double *)safeMalloc(ncolumns * sizeof(double));
    if (array[i] == NULL)
    {
      fprintf(stderr, "out of memory\n");
      exit(1);
    }
  }

  for (i = 0; i < nrows; i++)
  {
    for (j = 0; j < ncolumns; j++)
      array[i][j] = 0;
  }
  return array;
}

void free2array(double **array, int nrows)
{
  int i = 0;
  for (i = 0; i < nrows; i++)
    safeFree(array[i]);
  safeFree(array);
}
