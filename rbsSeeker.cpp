/********************************************************************
 * rbsSeeker: identify RBP binding sites at single-base resolution
 * Author: jianhua yang
 * Email: yangjh7@mail.sysu.edu.cn
 * Copyright: School of Life Sciences, Sun Yat-sen University
 * $ 2024/12/09
 ********************************************************************/
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

#include "BamReader.h"
#include "BamAux.h"

using namespace std;
using namespace BamTools;

#include "bioUtils.h"
#include "faiFile.h"
#include "bedFile.h"
#include "varFile.h"
#include "samFile.h"
#include "homer_statistics.h"
#include "statistic.h"
#include "rbsSeeker.h"

const string allSeekerType[SEEKER_TYPE] = {"Peak", "Conversion", "Mutation", "Deletion", "Insertion", "Truncation", "End"};
const string allClipType[CLIP_TYPE] = {"HITS", "PAR", "eCLIP", "iCLIP", "irCLIP", "miCLIP"};
chromCodeMap chromCodeHash;
chromCodeVector chromCodeList;

int cmpSite(const RBSite *x, const RBSite *y)
{
	return x->pval < y->pval;
}

int cmpChromLocus(const RBSite *x, const RBSite *y)
{
	if (x->chromCode != y->chromCode)
		return (x->chromCode < y->chromCode);

	return (x->start < y->start);
}

void getSiteBhFdr(siteVector &BHhash) {
	// q=p*n/rank log10(q) = log10(p)+log10(n)-log10(rank)
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

void seekRbsSites(FILE *gfp, FILE *faifp, char *outdir, struct parameterInfo *paraInfo, char *treatFile, char *controlFile)
{
	int mapLineNum         = 0;
	int skipSplice         = paraInfo->skipSplice;
	int normalization      = paraInfo->norm;
	double treatTotalNum   = 0;
	double controlTotalNum = 0;
	int keepDup = 1;
	if (paraInfo->PCR) keepDup = 0;

	faidxMap faiHash;
	chromSeqMap chromSeqHash;
	chromVarMap treatVarHash;
	chromVarMap controlVarHash;

	BamReader treatReader;
	BamReader controlReader;
	openBamFile(treatFile, treatReader);

	if (controlFile != NULL) openBamFile(controlFile, controlReader);

	time_t start = 0, end = 0, begin = 0;
	time(&start);
	time(&begin);

	fprintf(stderr, "read genome fai file\n");
	paraInfo->genomeSize = readFai(faifp, faiHash);
	fprintf(stderr, "encode chromosomes\n");
	encodeChrom(faiHash, chromCodeHash, chromCodeList);
	if (paraInfo->transcriptomeSize > paraInfo->genomeSize) paraInfo->transcriptomeSize = paraInfo->genomeSize;
	time(&end);
	fprintf(stderr, "elasped %.2f seconds for reading genome fai file\n", (double)(end - start));

	fprintf(stderr, "fetech fasta sequence from genome file\n");
	fetchAllSeq(gfp, faiHash, chromSeqHash);

	time(&end);
	fprintf(stderr, "elasped %.2f seconds for reading genome fasta file\n", (double)(end - start));
	time(&start);

	fprintf(stderr, "read treatment bam to var list\n");
	treatTotalNum = readBamToVariationMap(chromSeqHash, treatReader, treatVarHash, keepDup,
	                                      paraInfo->maxLocusNum, skipSplice, normalization, paraInfo->minReadLen);
	paraInfo->totalTreatNum = treatTotalNum;
	time(&end);
	fprintf(stderr, "elasped %.2f seconds for reading %.0f reads to variation list in treatment\n", (double)(end - start), treatTotalNum);

	if (controlFile != NULL)
	{
		time(&start);
		fprintf(stderr, "read control bam to sam list\n");
		controlTotalNum = readBamToVariationMap(chromSeqHash, controlReader, controlVarHash, keepDup,
		                                        paraInfo->maxLocusNum, skipSplice, normalization, paraInfo->minReadLen);
		paraInfo->totalCtrlNum = controlTotalNum;
		paraInfo->treatVsCtrlRatio = paraInfo->totalTreatNum / paraInfo->totalCtrlNum;
		time(&end);
		fprintf(stderr, "elasped %.2f seconds for reading %.0f reads to variation list in control\n", (double)(end - start), controlTotalNum);
	}

	fprintf(stderr, "safe free all sequence memories\n");
	safeFreeChromSeq(chromSeqHash);

	time(&start);
	fprintf(stderr, "idenfying peaks and mutations in treatment\n");
	callPeakAndMutation(paraInfo, treatVarHash, controlVarHash, gfp, faiHash, outdir);
	time(&end);
	fprintf(stderr, "elasped %.2f seconds for calling peaks and mutations in treatment\n", (double)(end - start));

	time(&start);
	fprintf(stderr, "output binding sites\n");
	//outputSiteInfo(gfp, faiHash, outdir, paraInfo, totalNum, contigVector);
	time(&end);
	fprintf(stderr, "elasped %.2f seconds for outputing binding sites\n", (double)(end - start));

	fprintf(stderr, "free all memories\n");
	treatReader.Close();
	//freechromVarMap(treatVarHash);

	if (controlFile != NULL) {
		controlReader.Close();
		//freechromVarMap(controlVarHash);
	}
	time(&end);
	fprintf(stderr, "total elasped %.2f seconds for whole workflow\n", (double)(end - begin));
}

void callPeakAndMutation(struct parameterInfo *paraInfo, chromVarMap &treatVarHash,
                         chromVarMap &controlVarHash, FILE *gfp, faidxMap &faiHash, char *outdir)
{
	siteTypeMap treatTypeHash;
	rbsSumTypeMap treatRbsSumHash;
	rbsSumTypeMap ctrlRbsSumHash;
	initiateRbsSumTypeMap(treatRbsSumHash, paraInfo->clipType);
	initiateRbsSumTypeMap(ctrlRbsSumHash, paraInfo->clipType);

	chromVarMap::iterator it;
	faidx *fai = NULL;
	for (it = treatVarHash.begin(); it != treatVarHash.end(); ++it)
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
		varVector treatVarList = it->second;
		varVector ctrlVarList;
		if (controlVarHash.find(it->first) != controlVarHash.end())
		{
			ctrlVarList = controlVarHash[it->first];
		}
		if (treatVarList.size() >= 1)
		{
			int chromLen   = fai->len;
			char* chromSeq = faidxFetchSeq(gfp, fai, 0, chromLen, '+');
			if (paraInfo->verbose) fprintf(stderr, "sort the chromosome: %s with %d reads in treatment\n", chromName, (int)treatVarList.size());
			sort(treatVarList.begin(), treatVarList.end(), compareVar);
			varProfileMap treatPlusVarProHash;
			varProfileMap treatMinusVarProHash;

			varProfileMap ctrlPlusVarProHash;
			varProfileMap ctrlMinusVarProHash;
			// for plus strand
			allocVarProfileMap(treatPlusVarProHash, chromLen, paraInfo->clipType);
			if (paraInfo->verbose) fprintf(stderr, "get peaks and mutations from + strand in treatment from %s\n", chromName);
			varToProfile(paraInfo, treatVarList, treatPlusVarProHash, chromSeq, chromLen, '+'); // plus strand
			if (ctrlVarList.size() > 0)
			{
				if (paraInfo->verbose) fprintf(stderr, "sort the chromosome: %s with %d reads in control\n", chromName, (int)ctrlVarList.size());
				sort(ctrlVarList.begin(), ctrlVarList.end(), compareVar);
				allocVarProfileMap(ctrlPlusVarProHash, chromLen, paraInfo->clipType);
				if (paraInfo->verbose) fprintf(stderr, "get peaks and mutations from + strand in control from %s\n", chromName);
				varToProfile(paraInfo, ctrlVarList, ctrlPlusVarProHash, chromSeq, chromLen, '+'); // plus strand
			}
			if (paraInfo->verbose) fprintf(stderr, "identifying rbs sites from + strand in treatment from %s\n", chromName);
			identifyRbsSites(paraInfo, treatPlusVarProHash, ctrlPlusVarProHash, treatTypeHash, chromName, chromLen, '+', treatRbsSumHash, ctrlRbsSumHash);
			freeVarProfileMap(treatPlusVarProHash);
			if (ctrlPlusVarProHash.size() > 0) freeVarProfileMap(ctrlPlusVarProHash);

			// for minus strand
			allocVarProfileMap(treatMinusVarProHash, chromLen, paraInfo->clipType);
			if (paraInfo->verbose) fprintf(stderr, "get peaks and mutations from - strand in treatment from %s\n", chromName);
			varToProfile(paraInfo, treatVarList, treatMinusVarProHash, chromSeq, chromLen, '-'); // minus strand
			if (ctrlVarList.size() > 0)
			{
				allocVarProfileMap(ctrlMinusVarProHash, chromLen, paraInfo->clipType);
				if (paraInfo->verbose) fprintf(stderr, "get peaks and mutations from - strand in control from %s\n", chromName);
				varToProfile(paraInfo, ctrlVarList, ctrlMinusVarProHash, chromSeq, chromLen, '-'); // minus strand
			}
			if (paraInfo->verbose) fprintf(stderr, "identifying rbs sites from - strand in treatment from %s\n", chromName);
			identifyRbsSites(paraInfo, treatMinusVarProHash, ctrlMinusVarProHash, treatTypeHash, chromName, chromLen, '-', treatRbsSumHash, ctrlRbsSumHash);
			freeVarProfileMap(treatMinusVarProHash);
			if (ctrlMinusVarProHash.size() > 0) freeVarProfileMap(ctrlMinusVarProHash);

			safeFree(chromSeq);
		}
		if (treatVarList.size() >= 1)
		{
			freeVarVector(treatVarList);
		}
		if (ctrlVarList.size() >= 1)
		{
			freeVarVector(ctrlVarList);
		}
	} // for bed hash
	fprintf(stderr, "#calculate genome-wide variations in treatment\n");
	calculateGenomeVariation(paraInfo, treatRbsSumHash);
	fprintf(stderr, "calculate pvalue in treatment\n");
	calculatePval(paraInfo, treatTypeHash, treatRbsSumHash, ctrlRbsSumHash);
	fprintf(stderr, "output all RBS sites\n");
	outputAllTypeRbsSite(paraInfo, treatTypeHash, gfp, faiHash, outdir);
	fprintf(stderr, "free memory of all sites\n");
	freeSiteTypeMap(treatTypeHash);
	freeRbsSumTypeMap(treatRbsSumHash);
	freeRbsSumTypeMap(ctrlRbsSumHash);
}

void initiateRbsSumTypeMap(rbsSumTypeMap &rbsSumHash, int clipType)
{
	int i = 0;
	double initVal = MIN_RATIO;
	for (i = 0; i < SEEKER_TYPE; i++)
	{
		rbsSumType *rbsSum   = NULL;
		rbsSum               = (rbsSumType *)safeMalloc(sizeof(rbsSumType));
		rbsSum->totalReadNum = initVal;
		rbsSum->totalHeight  = initVal;
		rbsSum->totalLen     = 1;
		rbsSumHash[allSeekerType[i]] = rbsSum;
	}
}

void freeRbsSumTypeMap(rbsSumTypeMap &rbsSumHash)
{
	rbsSumTypeMap::iterator it;
	for (it = rbsSumHash.begin(); it != rbsSumHash.end(); ++it)
	{
		rbsSumType *rbsSum = it->second;
		if (rbsSum != NULL) safeFree(rbsSum);
	}
	rbsSumHash.clear();
}

void allocVarProfileMap(varProfileMap &varProHash, int span, int clipType)
{
	int i = 0;
	int j = 0;
	double pseudoProCount = 0;
	for (i = 0; i < SEEKER_TYPE; i++)
	{
		double *profile = NULL;
		profile = (double *)safeMalloc(sizeof(double) * span);
		for (j = 0; j < span; j++)
		{
			profile[j] = pseudoProCount;
		}
		varProHash[allSeekerType[i]] = profile;
	}
}

void freeVarProfileMap(varProfileMap &varProHash)
{
	varProfileMap::iterator it;
	for (it = varProHash.begin(); it != varProHash.end(); ++it)
	{
		double *profile = it->second;
		if (profile != NULL)  safeFree(profile);
	}
	varProHash.clear();
}

double varToProfile(struct parameterInfo *paraInfo, varVector &varList, varProfileMap &varProHash,
                    char *chromSeq, int chromLen, char strand)
{
	int i            = 0;
	int j            = 0;
	int index        = 0;
	double totalNum  = 0;
	int  readLen     = 0;
	int  start       = 0;
	int  end         = chromLen;
	int span         = end - start;

	double *truProfile = varProHash["Truncation"];
	double *endProfile = varProHash["End"];
	for (varVector::iterator vecItr = varList.begin(); vecItr != varList.end(); vecItr++)
	{
		Variation *var = *vecItr;
		double readNum = var->readNum;
		totalNum += readNum;
		if (var->strand != strand) // skip other strand
		{
			continue;
		}

		// transform read profile into chromosome profile
		int seMut = getVariation(paraInfo, varProHash, var, paraInfo->clipType, strand);
		// for start and end
		if (strand == '+')
		{
			int startIdx = var->chromStart - start;
			int endIdx   = var->chromEnd - start - 1;
			if (startIdx >= 0 && startIdx < span)
			{
				if (seMut != 1) truProfile[startIdx] += readNum;
			}
			if (endIdx >= 0 && endIdx < span) {
				if (seMut != 2) endProfile[endIdx]   += readNum;
			}
		}
		else
		{
			int startIdx = var->chromEnd - start - 1;
			int endIdx   = var->chromStart - start;
			if (startIdx >= 0 && startIdx < span)
			{
				if (seMut != 1) truProfile[startIdx] += readNum;
			}
			if (endIdx >= 0 && endIdx < span) {
				if (seMut != 2) endProfile[endIdx]   += readNum;
			}
		}
	}// for end
	return totalNum;
}

int getVariation(struct parameterInfo *paraInfo, varProfileMap &varProHash, Variation *var, int clipType, char strand)
{
	int i = 0;
	int j = 0;
	int seMut = 0;
	int varNum          = var->mutNum;
	uint32_t chromStart = var->chromStart;
	uint32_t chromEnd   = var->chromEnd;
	uint32_t exonStart  = chromStart;
	uint32_t exonEnd    = chromEnd;
	double readNum = var->readNum;
	exonVector exonList;

	double *conProfile = varProHash["Conversion"];
	double *mutProfile = varProHash["Mutation"];
	double *delProfile = varProHash["Deletion"];
	double *insProfile = varProHash["Insertion"];
	double *peaProfile = varProHash["Peak"];

	for (i = 0; i < varNum;)
	{
		uint32_t pos = var->position[i];
		uint8_t  mut = var->mutation[i];
		if (paraInfo->cvs != NULL && paraInfo->cvsCode == mut)
		{
			//conProfile[pos] += readNum;
			if (paraInfo->rmSeMutation)
			{
				if (pos >= chromStart + mutDist && pos <= chromEnd - mutDist - 1) conProfile[pos] += readNum;
			}
			else
			{
				conProfile[pos] += readNum;
			}
		}
		if (mut >= 4 && mut <= 15)
		{
			if ((pos == chromStart && strand == '+') || (pos == chromEnd - 1 && strand == '-'))
			{
				seMut = 1;
			}
			if ((pos == chromEnd - 1 && strand == '+') || (pos == chromStart && strand == '-'))
			{
				seMut = 2;
			}
			if (paraInfo->rmSeMutation)
			{
				if (pos >= chromStart + mutDist && pos <= chromEnd - mutDist - 1) mutProfile[pos] += readNum;
			}
			else
			{
				mutProfile[pos] += readNum;
			}

		}
		if (mut == 16)
		{
			//delProfile[pos] += readNum;
			if (paraInfo->rmSeMutation)
			{
				if (pos >= chromStart + mutDist && pos <= chromEnd - mutDist - 1) delProfile[pos] += readNum;
			}
			else
			{
				delProfile[pos] += readNum;
			}
		}
		if (mut == 17)
		{
			//insProfile[pos] += readNum;
			if (paraInfo->rmSeMutation)
			{
				if (pos >= chromStart + mutDist && pos <= chromEnd - mutDist - 1) insProfile[pos] += readNum;
			}
			else
			{
				insProfile[pos] += readNum;
			}
		}
		if (mut == 18)
		{
			exonEnd = pos;
			exonPair ep = make_pair(exonStart, exonEnd);
			exonList.push_back(ep);
			exonStart = var->position[i + 1];
			i += 1;
		}
		i++;
	}
	exonEnd     = chromEnd;
	exonPair ep = make_pair(exonStart, exonEnd);
	exonList.push_back(ep);
	for (i = 0; i < exonList.size(); i++)
	{
		uint32_t start = exonList[i].first;
		uint32_t end   = exonList[i].second;
		for (j = start; j < end; j++)
		{
			peaProfile[j] += readNum;
		}
	}
	exonList.clear();

	return seMut;
}

void identifyRbsSites(struct parameterInfo *paraInfo, varProfileMap &treatVarProHash, varProfileMap &ctrlVarProHash,
                      siteTypeMap &allSiteTypeHash, char *chrom, int chromLen, char strand,
                      rbsSumTypeMap &treatRbsSumHash, rbsSumTypeMap &ctrlRbsSumHash)
{
	int i = 0;
	int j = 0;
	int k = 0;
	int peakStart = 0;
	int peakEnd   = 0;
	double maxCov = 0;
	int    maxPos = 0;
	int stepLen   = 5000;
	string chromStr(chrom);
	uint16_t chromCode = chromCodeHash[chromStr];
	string peakType    = allSeekerType[0];
	double *treatHgt   = treatVarProHash[peakType];
	double *ctrlHgt    = NULL;
	if (ctrlVarProHash.size() > 0) ctrlHgt = ctrlVarProHash[peakType];
	rbsSumType *tpRbsSum = treatRbsSumHash[peakType];
	rbsSumType *cpRbsSum = ctrlRbsSumHash[peakType];
	siteVector peakSiteList = allSiteTypeHash[peakType];

	int ctrlFlag = 0;
	double ctrlRatio = paraInfo->totalTreatNum / paraInfo->totalCtrlNum;
	if (paraInfo->totalCtrlNum > minCtrlNum) ctrlFlag = 1;

	// for peaks
	for (i = 1; i < chromLen - 1; i++)
	{
		// for peaks
		int start = i;
		int end   = start + 1;
		double tHeight    = treatHgt[i];
		double tHeightRpm = tHeight / paraInfo->totalTreatNum * rpmFactor;
		double cHeight    = 0;
		if (ctrlHgt != NULL) cHeight = ctrlHgt[i];
		double cHeightRpm = 0;
		if (ctrlHgt != NULL) cHeightRpm = cHeight / paraInfo->totalCtrlNum * rpmFactor;
		if (tHeight > 0) {
			tpRbsSum->totalHeight  += tHeight;
			tpRbsSum->totalReadNum += tHeight;
			tpRbsSum->totalLen     += 1;
		}
		if (cHeight > 0) {
			cpRbsSum->totalHeight  += cHeight;
			cpRbsSum->totalReadNum += cHeight;
			cpRbsSum->totalLen     += 1;
		}

		if (tHeight > 0)
		{
			if (peakStart == 0)
			{
				peakStart = start;
				peakEnd   = end;
			}
			else
			{
				peakEnd   = end;
			}
			if (tHeight > maxCov)
			{
				maxCov = tHeight;
				maxPos = i;
			}
		}
		else
		{
			int peakLen = peakEnd - peakStart;
			if (peakStart > 0 && peakLen > 0)
			{
				identifyPeaks(paraInfo, treatHgt, ctrlHgt, chromCode, chromLen, peakStart,
				              peakEnd, strand, maxCov, maxPos, ctrlRatio, peakSiteList);
			}
			peakStart = 0;
			peakEnd   = 0;
			maxCov = 0;
			maxPos = 0;
		}
		if (i == (chromLen - 1) )
		{
			int peakLen = peakEnd - peakStart;
			if (peakStart > 0 && peakLen > 0)
			{
				identifyPeaks(paraInfo, treatHgt, ctrlHgt, chromCode, chromLen, peakStart,
				              peakEnd, strand, maxCov, maxPos, ctrlRatio, peakSiteList);
			}
			peakStart = 0;
			peakEnd   = 0;
			maxCov = 0;
			maxPos = 0;
		}
	}
	allSiteTypeHash[peakType] = peakSiteList;
	// other site types
	for (j = 1; j < SEEKER_TYPE; j++)
	{
		string siteType = allSeekerType[j];
		rbsSumType *tRbsSum = treatRbsSumHash[siteType];
		rbsSumType *cRbsSum = ctrlRbsSumHash[siteType];
		if (treatVarProHash.find(siteType) == treatVarProHash.end()) continue;
		double *treatVarHgt = treatVarProHash[siteType];
		double *ctrlVarHgt  = NULL;
		if (ctrlVarProHash.size() > 0) ctrlVarHgt = ctrlVarProHash[siteType];
		siteVector mutSiteList = allSiteTypeHash[siteType];
		for (i = 1; i < chromLen - 1; i++)
		{
			// for peaks
			int start = i;
			int end   = start + 1;
			double tHeight    = treatHgt[i];
			double tHeightRpm = tHeight / paraInfo->totalTreatNum * rpmFactor;
			double cHeight   = 0;
			if (ctrlHgt != NULL) cHeight = ctrlHgt[i];
			double tvar = treatVarHgt[i];
			double uvar = treatVarHgt[i - 1];
			double dvar = treatVarHgt[i + 1];
			double uHeight = treatHgt[i - 1];
			double dHeight = treatHgt[i + 1];
			if (strand == '-')
			{
				double tmpVar = uvar;
				double tmpHgt = uHeight;
				uvar = dvar;
				dvar = tmpVar;
				uHeight = dHeight;
				dHeight = tmpHgt;
			}
			if (uvar <= 0) uvar = 0;
			if (dvar <= 0) dvar = 0;

			if (uHeight <= 0) uHeight = 1;
			if (dHeight <= 0) dHeight = 1;

			double cvar = 0;
			if (ctrlVarHgt != NULL) cvar = ctrlVarHgt[i];
			if (tvar > 0 && tHeight >= paraInfo->minHeight)
			{
				tRbsSum->totalHeight  += tHeight;
				tRbsSum->totalReadNum += tvar;
				tRbsSum->totalLen += 1;
			}
			if (cvar > 0)
			{
				if (cHeight >= paraInfo->minHeight) {
					cRbsSum->totalHeight  += cHeight;
					cRbsSum->totalReadNum += cvar;
					cRbsSum->totalLen += 1;
				}
			}
			if ( tvar >= 1 && tHeight >= paraInfo->minHeight && tHeightRpm >= paraInfo->minRpm)
			{
				if (uvar > cvar)
				{
					cvar = uvar;
					cHeight = uHeight;
				}
				if (dvar > cvar)
				{
					cvar = dvar;
					cHeight = dHeight;
				}
				RBSite *rbs   = allocRBSite(chromCode, start, end, strand);
				rbs->tVarNum  = tvar;
				rbs->tHeight  = tHeight;
				rbs->cVarNum  = cvar;
				rbs->cHeight  = cHeight;
				mutSiteList.push_back(rbs);
			}
		}// for chrom length
		allSiteTypeHash[siteType] = mutSiteList;
	} // for site type
}

int identifyPeaks(struct parameterInfo *paraInfo, double *treatHgt, double *ctrlHgt,
                  uint16_t chromCode, int chromLen, int peakStart, int peakEnd, char strand,
                  double maxCov, int maxPos, double ctrlRatio, siteVector &siteList)
{
	int i = 0;
	int peakTag = 0;
	int peakLen = peakEnd - peakStart;
	double tSumHeight = 0;
	double cSumHeight = 0;
	int orgPeakStart = peakStart;
	int orgPeakEnd   = peakEnd;
	int treatPeakLen = 0;
	int ctrlPeakLen  = 0;
	orgPeakStart = orgPeakStart - paraInfo->windowLen;
	orgPeakEnd   = orgPeakEnd + paraInfo->windowLen;
	if (orgPeakStart < 0) orgPeakStart = 0;
	if (orgPeakEnd > chromLen) orgPeakEnd = chromLen;
	double ctrlCov = 0;
	if (ctrlHgt != NULL) ctrlCov = ctrlHgt[maxPos];
	if (maxPos - peakExtLen > orgPeakStart)
	{
		peakStart = maxPos - peakExtLen;
	}
	if (maxPos + peakExtLen < orgPeakEnd)
	{
		peakEnd = maxPos + peakExtLen;
	}

	for (i = orgPeakStart; i < orgPeakEnd; i++)
	{
		if (treatHgt[i] > 0)
		{
			tSumHeight += treatHgt[i];
			treatPeakLen += 1;
		}
		if (ctrlHgt != NULL && ctrlHgt[i] > 0)
		{
			cSumHeight += ctrlHgt[i];
			ctrlPeakLen += 1;
		}
	}
	int newPeakLen = peakEnd - peakStart;
	if (newPeakLen > 0)
	{
		RBSite *rbs   = allocRBSite(chromCode, peakStart, peakEnd, strand);
		double tMeanPeak = tSumHeight / (double)treatPeakLen;
		double cMeanPeak = 0;
		if (ctrlPeakLen > 0) cMeanPeak = cSumHeight / (double)ctrlPeakLen * ctrlRatio;
		ctrlCov = ctrlCov * ctrlRatio;
		rbs->tVarNum  = maxCov;
		rbs->tHeight  = maxCov;
		if (tMeanPeak > cMeanPeak) cMeanPeak = tMeanPeak;
		if (ctrlCov > cMeanPeak) cMeanPeak = ctrlCov;
		rbs->cVarNum  = cMeanPeak;
		rbs->cHeight  = cMeanPeak;
		siteList.push_back(rbs);
		peakTag = 1;
	}
	return peakTag;
}

void calculatePval(struct parameterInfo *paraInfo, siteTypeMap &siteTypeHash,
                   rbsSumTypeMap &treatRbsSumHash, rbsSumTypeMap &ctrlRbsSumHash)
{
	int binomialFlag = 0;
	if (paraInfo->totalCtrlNum > minCtrlNum) binomialFlag = 1;
	double ctrlRatio = 1;
	rbsSumType *tRbsSum  = treatRbsSumHash["Peak"];
	rbsSumType *cRbsSum  = ctrlRbsSumHash["Peak"];
	if (paraInfo->transcriptomeSize < tRbsSum->totalLen)
	{
		tRbsSum->totalLen = paraInfo->transcriptomeSize;
	}
	double lambda = tRbsSum->totalHeight / (double)tRbsSum->totalLen;
	if (binomialFlag)
	{
		double ctrlLambda = cRbsSum->totalHeight / (double)cRbsSum->totalLen;
		ctrlRatio = lambda / ctrlLambda;
		fprintf(stderr, "Ratio:%.5f\ttreatLambda:%.5f\tctrlLambda:%.5f\n", ctrlRatio, lambda, ctrlLambda);
		fprintf(stderr, "treatTotalHeigth:%.5f\tctrlTotalHeigth:%.5f\ttreatTotalLen:%ld\tctrlTotalLen:%ld\n", tRbsSum->totalHeight, cRbsSum->totalHeight, tRbsSum->totalLen, cRbsSum->totalLen);
	}
	else
	{
		fprintf(stderr, "treatLambda:%.5f\n", lambda);
	}
	for (siteTypeMap::iterator mapItr = siteTypeHash.begin(); mapItr != siteTypeHash.end(); mapItr++) {
		string siteType     = mapItr->first;
		siteVector siteList = mapItr->second;
		rbsSumType *rbsSum  = treatRbsSumHash[siteType];
		for (siteVector::iterator vecItr = siteList.begin(); vecItr != siteList.end(); vecItr++)
		{
			RBSite *rbs = *vecItr;
			if (siteType == "Peak")
			{
				double ipRead = rbs->tHeight;
				double inRead = lambda;
				if (rbs->cHeight > inRead)
				{
					inRead = rbs->cHeight;
				}
				rbs->cHeight = inRead;
				rbs->cVarNum = inRead;
				if (ipRead < pseudoCount) ipRead = pseudoCount;
				if (inRead < pseudoCount) inRead = pseudoCount;
				rbs->mfold    = ipRead / inRead;
				double logPoissonP = ilogCumulativePoisson(ipRead, inRead);
				rbs->pval = log10Lnp(logPoissonP);
			}
			else
			{
				//double totalNum = rbs->cHeight + rbs->height;
				//double stopNum  = rbs->cReadNum + rbs->readNum;
				//double pval     = loghypergeoD(totalNum, stopNum, rbs->height, rbs->readNum);
				int tHeight = (int) (rbs->tHeight + 0.5);
				int cHeight = (int) (rbs->cHeight + 0.5);
				int tVarNum = (int) (rbs->tVarNum + 0.5);
				int cVarNum = (int) (rbs->cVarNum + 0.5);
				//int totalNum = treatHeight + ctrlHeight;
				//int stopNum  = treatNum + ctrlNum;
				//double pval  = loghypergeoD(totalNum, stopNum, treatHeight, treatNum);
				double tRatio = rbs->tVarNum / rbs->tHeight;
				double cRatio = rbsSum->varRatio;
				//if (binomialFlag)
				//{
				cRatio = 1.0 / (double)tHeight;
				//if (cRatio < rbsSum->varRatio) cRatio = rbsSum->varRatio;
				if (cHeight > 5 ) {
					if (cVarNum > 0)
						cRatio = (double)cVarNum / (double)cHeight;
					else
						cRatio = (double)(cVarNum + 0.1) / (double)(cHeight + 1.0);
				}
				//}
				rbs->mfold  = tRatio / cRatio;
				double pval = 1.0;
				if (binomialFlag)
				{
					int maxInt  = tHeight * 20;
					pval = logbinomialD(tHeight, tVarNum, cRatio, maxInt);
				}
				else
				{
					if (cRatio < rbsSum->varRatio) cRatio = rbsSum->varRatio;
					rbs->mfold  = tRatio / cRatio;
					int maxInt  = tHeight * 20;
					pval = logbinomialD(tHeight, tVarNum, cRatio, maxInt);
					//pval = loghypergeoD(rbsSum->totalHeight, rbsSum->totalReadNum, tHeight, tVarNum);
				}
				rbs->pval = log10Lnp(pval);
			}
		}
	}
}

void calculateGenomeVariation(struct parameterInfo *paraInfo, rbsSumTypeMap &treatRbsSumHash)
{
	int i = 0;
	for (i = 1; i < SEEKER_TYPE; i++)
	{
		string siteType    = allSeekerType[i];
		rbsSumType *rbsSum = treatRbsSumHash[siteType];

		rbsSum->varRatio = rbsSum->totalReadNum / rbsSum->totalHeight;
		if (siteType == "Conversion")
		{
			if (paraInfo->cvs != NULL) fprintf(stderr, "the number of %c->%c %s is %.0f, total number is %.0f, ratio is %.5f\n", paraInfo->cvs[0], paraInfo->cvs[1], siteType.c_str(), rbsSum->totalReadNum, rbsSum->totalHeight, rbsSum->varRatio);
		}
		else
		{
			fprintf(stderr, "the number of %s is %.0f, total number is %.0f, ratio is %.5f\n", siteType.c_str(), rbsSum->totalReadNum, rbsSum->totalHeight, rbsSum->varRatio);
		}
	}
}

void outputAllTypeRbsSite(struct parameterInfo *paraInfo, siteTypeMap &treatSiteTypeHash,
                          FILE *gfp, faidxMap &faiHash, char *outdir)
{
	for (siteTypeMap::iterator mapItr = treatSiteTypeHash.begin(); mapItr != treatSiteTypeHash.end(); mapItr++) {
		string siteType          = mapItr->first;
		siteVector treatSiteList = mapItr->second;
		char *siteTypeStr        =  (char *)siteType.c_str();
		if (siteType == "Peak")
		{
			if (treatSiteList.size() > 0) outputPeakTypeSites(paraInfo, gfp, faiHash, treatSiteList, outdir, siteTypeStr);
		}
		else
		{
			if (treatSiteList.size() > 0) outputVariationTypeSites(paraInfo, gfp, faiHash, treatSiteList, outdir, siteTypeStr);
		}
	}
}

int outputPeakTypeSites(struct parameterInfo * paraInfo, FILE * gfp, faidxMap & faiHash,
                        siteVector & siteArray, char *outdir, char *siteType)
{
	int extendLen = EXTEND_LEN;
	int siteIdx = 1;
	char outfile[512];
	strcpy(outfile, outdir);
	strcat(outfile, "_");
	strcat(outfile, siteType);
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
	fprintf(outfp, "extendSeq\tmotifPos\ttype\tlog10(p-value)\tlog10(q-value)\tvariationNum\theight\theightRpm\tratio\tmfold\tctrlNum\tctrlHeight");
	if (paraInfo->rnafold)
	{
		fprintf(outfp, "\tstructure\tmfe\n");
	}
	else
	{
		fprintf(outfp, "\n");
	}
	fflush(outfp);
	sort(siteArray.begin(), siteArray.end(), cmpChromLocus);
	siteVector mergeSiteArray;
	mergeAllPeakSite(paraInfo, siteArray, mergeSiteArray, '+');
	mergeAllPeakSite(paraInfo, siteArray, mergeSiteArray, '-');
	sort(mergeSiteArray.begin(), mergeSiteArray.end(), cmpSite);
	for (siteVector::iterator vecItr = mergeSiteArray.begin(); vecItr != mergeSiteArray.end(); vecItr++)
	{
		RBSite *rbs = *vecItr;
		outputOnePeakSite(paraInfo, outfp, gfp, faiHash, rbs, siteIdx, siteType);
		siteIdx++;
	}
	fprintf(stderr, "get %d\t%ss\tfrom %lu sites\n", (siteIdx - 1), siteType, siteArray.size());
	freeSiteVector(mergeSiteArray);
	fclose(outfp);
	return siteIdx;
}

void mergeAllPeakSite(struct parameterInfo * paraInfo, siteVector & siteArray, siteVector & mergeSiteArray, char strand)
{
	// output + strand
	RBSite *rbsPeak = NULL;
	int firstTag = 0;
	int maxDist = mergePeakDist; // the maximum distance between two sites
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
			if (rbsPeak->chromCode == rbs->chromCode && overlapLength(rbsPeak->start, rbsPeak->end, rbs->start, rbs->end) >= maxDist)
			{
				int minStart = MIN(rbsPeak->start, rbs->start);
				int maxEnd   = MAX(rbsPeak->end, rbs->end);
				if ((rbsPeak->qval > rbs->qval) || (rbsPeak->qval == rbs->qval && rbsPeak->mfold < rbs->mfold))
				{
					copyRBSiteNoChrom(rbsPeak, rbs);
				}
				rbsPeak->start = minStart;
				rbsPeak->end   = maxEnd;
			}
			else {
				if (rbsPeak->end - rbsPeak->start >= paraInfo->minClusterLen)
				{
					mergeSiteArray.push_back(rbsPeak);
				}
				rbsPeak = NULL;
				rbsPeak = (RBSite *)safeMalloc(sizeof(RBSite));
				copyRBSite(rbsPeak, rbs);
			}
		}
	}// for site vectors
	if (rbsPeak != NULL)
	{
		if (rbsPeak->end - rbsPeak->start >= paraInfo->minClusterLen)
		{
			mergeSiteArray.push_back(rbsPeak);
		}
	}
}

int outputOnePeakSite(struct parameterInfo * paraInfo, FILE * outfp, FILE * gfp, faidxMap & faiHash,
                      RBSite *rbs, int siteIdx, char *siteType)
{
	uint16_t chromCode = rbs->chromCode;
	int chromStart  = rbs->start;
	int chromEnd    = rbs->end;
	char strand     = rbs->strand;
	faidx *fai = NULL;
	string chrom = chromCodeList[chromCode];
	fai = faiHash[chrom];
	int extendStart = chromStart;
	int extendEnd = chromEnd;
	if (extendStart < 0) extendStart = 0;
	if (extendEnd > fai->len) extendEnd = fai->len;
	char *extendSeq   = (char *)faidxFetchSeq(gfp, fai, extendStart, extendEnd, strand);
	int motifPos = searchMotif(paraInfo->motifSeq, extendSeq);
	double heightRpm = rbs->tHeight / paraInfo->totalTreatNum * 1000000;
	double ratio = rbs->tVarNum / rbs->tHeight;
	if (rbs->pval < paraInfo->pval
	        && rbs->qval < paraInfo->qval
	        && rbs->mfold >= paraInfo->mfold)
	{
		fprintf(outfp, "%s\t%d\t%d\trbsSeeker_%s_%d\t%.0f\t%c\t%s\t%d\t", chrom.c_str(), chromStart, chromEnd, siteType, siteIdx, rbs->tHeight, strand, extendSeq, motifPos);
		fprintf(outfp, "%s\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t", siteType, rbs->pval, rbs->qval, rbs->tVarNum, rbs->tHeight, heightRpm, ratio, rbs->mfold);
		fprintf(outfp, "%.5f\t%.5f", rbs->cVarNum, rbs->cHeight);
		if (paraInfo->rnafold)
		{
			char *structure = (char *)safeMalloc(strlen(extendSeq) + 1);
			double mfe = fold(extendSeq, structure);
			fprintf(outfp, "\t%s\t%.5f\n", structure, mfe);
			free_arrays;
			safeFree(structure);
		}
		else
		{
			fprintf(outfp, "\n");
		}
	}
	safeFree(extendSeq);
	siteIdx++;
	fflush(outfp);
	return siteIdx;
}

void copyRBSiteNoChrom(RBSite * tRbs, RBSite * oRbs)
{
	tRbs->start   = oRbs->start;
	tRbs->end     = oRbs->end;
	tRbs->strand  = oRbs->strand;
	tRbs->tVarNum = oRbs->tVarNum;
	tRbs->tHeight = oRbs->tHeight;
	tRbs->cVarNum = oRbs->cVarNum;
	tRbs->cHeight = oRbs->cHeight;
	tRbs->mfold   = oRbs->mfold;
	tRbs->ufold   = oRbs->ufold;
	tRbs->dfold   = oRbs->dfold;
	tRbs->pval    = oRbs->pval;
	tRbs->qval    = oRbs->qval;
}

void copyRBSite(RBSite * tRbs, RBSite * oRbs)
{
	tRbs->chromCode = oRbs->chromCode;
	tRbs->start   = oRbs->start;
	tRbs->end     = oRbs->end;
	tRbs->strand  = oRbs->strand;
	tRbs->tVarNum = oRbs->tVarNum;
	tRbs->tHeight = oRbs->tHeight;
	tRbs->cVarNum = oRbs->cVarNum;
	tRbs->cHeight = oRbs->cHeight;
	tRbs->mfold   = oRbs->mfold;
	tRbs->ufold   = oRbs->ufold;
	tRbs->dfold   = oRbs->dfold;
	tRbs->pval    = oRbs->pval;
	tRbs->qval    = oRbs->qval;
}

void mergeAllVariationSite(struct parameterInfo * paraInfo, siteVector & siteArray, siteVector & mergeSiteArray, char strand)
{
	// output + strand
	RBSite *rbsPeak = NULL;
	int firstTag = 0;
	int maxDist  = mergeVarDist; // the maximum distance between two sites
	int minStart = siteArray[0]->start;
	int maxEnd   = siteArray[0]->end;
	for (siteVector::iterator vecItr = siteArray.begin(); vecItr != siteArray.end(); vecItr++)
	{
		RBSite *rbs = *vecItr;
		if (rbs->strand != strand) continue;
		if (rbs->pval < paraInfo->pval && rbs->qval < paraInfo->qval && rbs->mfold >= paraInfo->mfold)
		{
			if (firstTag == 0)
			{
				firstTag = 1;
				minStart = rbs->start;
				maxEnd   = rbs->end;
				rbsPeak = (RBSite *)safeMalloc(sizeof(RBSite));
				copyRBSite(rbsPeak, rbs);
				continue;
			}
			if (rbsPeak->chromCode == rbs->chromCode && overlapLength(minStart, maxEnd, rbs->start, rbs->end) >= maxDist)
			{
				minStart = MIN(minStart, rbs->start);
				maxEnd   = MAX(maxEnd, rbs->end);
				if ((rbsPeak->qval > rbs->qval) || (rbsPeak->qval == rbs->qval && rbsPeak->mfold < rbs->mfold))
				{
					copyRBSiteNoChrom(rbsPeak, rbs);
				}
			}
			else {
				mergeSiteArray.push_back(rbsPeak);
				//rbsPeak = NULL;
				rbsPeak = (RBSite *)safeMalloc(sizeof(RBSite));
				copyRBSite(rbsPeak, rbs);
				minStart = rbs->start;
				maxEnd   = rbs->end;
			}
		}// for filter
	}// for site vectors
	if (rbsPeak != NULL)
	{
		mergeSiteArray.push_back(rbsPeak);
	}
}

int outputVariationTypeSites(struct parameterInfo * paraInfo, FILE * gfp, faidxMap & faiHash,
                             siteVector & siteArray, char *outdir, char *siteType)
{
	int siteIdx = 0;
	char outfile[512];
	strcpy(outfile, outdir);
	strcat(outfile, "_");
	if (strcmp(siteType, "Conversion") == 0 && paraInfo->cvs != NULL)
	{
		strcat(outfile, paraInfo->cvs);
		//strcat(outfile, "_");
		//strcat(outfile, siteType);
	}
	else
	{
		strcat(outfile, siteType);
	}
	strcat(outfile, ".bed");
	FILE *outfp = NULL;
	int i = 0;
	outfp = (FILE *) fopen(outfile, "w");
	if (outfp == NULL)
	{
		fprintf(stderr, "ERROR: Can't open %s\n", outfile);
		exit(1);
	}
	//fprintf(stderr, "calculate the Benjamini-Hochberg FDR\n");
	getSiteBhFdr(siteArray);
	fprintf(outfp, "#chrom\tchromStart\tchromEnd\tname\tscore\tstrand\t");
	fprintf(outfp, "extendSeq\tmotifPos\ttype\tlog10(p-value)\tlog10(q-value)\tvariationNum\theight\theightRpm\tratio\tmfold\tctrlNum\tctrlHeight");
	if (paraInfo->rnafold)
	{
		fprintf(outfp, "\tstructure\tmfe\n");
	}
	else
	{
		fprintf(outfp, "\n");
	}
	fflush(outfp);
	sort(siteArray.begin(), siteArray.end(), cmpChromLocus);
	siteVector mergeSiteArray;
	mergeAllVariationSite(paraInfo, siteArray, mergeSiteArray, '+');
	mergeAllVariationSite(paraInfo, siteArray, mergeSiteArray, '-');
	sort(mergeSiteArray.begin(), mergeSiteArray.end(), cmpSite);
	for (siteVector::iterator vecItr = mergeSiteArray.begin(); vecItr != mergeSiteArray.end(); vecItr++)
	{
		RBSite *rbs = *vecItr;
		outputOneVarSite(paraInfo, outfp, gfp, faiHash, rbs, siteIdx, siteType);
		siteIdx++;
	}
	if (strcmp(siteType, "Conversion") == 0 && paraInfo->cvs != NULL)
	{
		fprintf(stderr, "get significant %d\t%c->%c mutations\tfrom %lu sites\n", siteIdx, paraInfo->cvs[0], paraInfo->cvs[1], siteArray.size());
	}
	else
	{
		fprintf(stderr, "get significant %d\t%ss\tfrom %lu sites\n", siteIdx, siteType, siteArray.size());
	}
	freeSiteVector(mergeSiteArray);
	fclose(outfp);
	return siteIdx;
}

int outputOneVarSite(struct parameterInfo * paraInfo, FILE * outfp, FILE * gfp, faidxMap & faiHash,
                     RBSite *rbs, int siteIdx, char *siteType)
{
	int extendLen = EXTEND_LEN;
	double ratio = rbs->tVarNum / rbs->tHeight;
	double heightRpm = rbs->tHeight / paraInfo->totalTreatNum * 1000000;
	if (rbs->mfold >= paraInfo->mfold
	        && rbs->pval < paraInfo->pval
	        && rbs->qval < paraInfo->qval
	        && ratio >= paraInfo->minRatio
	        && ratio <= paraInfo->maxRatio
	        && rbs->tVarNum >= paraInfo->minMutNum)
	{
		uint16_t chromCode = rbs->chromCode;
		int chromStart  = rbs->start;
		int chromEnd    = rbs->end;
		char strand     = rbs->strand;
		faidx *fai = NULL;
		string chrom = chromCodeList[chromCode];
		fai = faiHash[chrom];
		int extendStart = chromStart - extendLen;
		int extendEnd = chromEnd + extendLen;
		if (extendStart < 0) extendStart = 0;
		if (extendEnd > fai->len) extendEnd = fai->len;
		char *extendSeq   = (char *)faidxFetchSeq(gfp, fai, extendStart, extendEnd, strand);
		int motifPos = searchMotif(paraInfo->motifSeq, extendSeq);
		fprintf(outfp, "%s\t%d\t%d\trbsSeeker_%s_%d\t%.0f\t%c\t%s\t%d\t", chrom.c_str(), chromStart, chromEnd, siteType, siteIdx, rbs->tHeight, strand, extendSeq, motifPos);
		fprintf(outfp, "%s\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t", siteType, rbs->pval, rbs->qval, rbs->tVarNum, rbs->tHeight, heightRpm, ratio, rbs->mfold);
		fprintf(outfp, "%.5f\t%.5f", rbs->cVarNum, rbs->cHeight);
		if (paraInfo->rnafold)
		{
			char *structure = (char *)safeMalloc(strlen(extendSeq) + 1);
			double mfe = fold(extendSeq, structure);
			fprintf(outfp, "\t%s\t%.5f\n", structure, mfe);
			free_arrays;
			safeFree(structure);
		}
		else
		{
			fprintf(outfp, "\n");
		}
		siteIdx++;
		safeFree(extendSeq);
	}
	return siteIdx;
}


void freeSiteTypeMap(siteTypeMap &siteTypeHash) /*free site map */
{
	for (siteTypeMap::iterator mapItr = siteTypeHash.begin(); mapItr != siteTypeHash.end(); mapItr++) {
		siteVector siteList = mapItr->second;
		freeSiteVector(siteList);
	}
	siteTypeHash.clear();
}

RBSite *allocRBSite(uint16_t chromCode, int start, int end, char strand)
{
	RBSite *rbs    = (RBSite *)safeMalloc(sizeof(RBSite));
	rbs->chromCode = chromCode;
	rbs->start     = start;
	rbs->end       = end;
	rbs->strand    = strand;
	rbs->tVarNum   = 0;
	rbs->tHeight   = 0;
	rbs->mfold     = 0;
	rbs->cHeight   = 0;
	rbs->cVarNum   = 0;
	rbs->pval      = 1.5;
	rbs->qval      = 1.5;
	return rbs;
}

char revComChar(char ch, char strand)
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
	return ch;
}

int encodeNucleotide(char base)
{
	int i = 0;
	base = toupper(base);
	if (base == 'A') return 0;
	else if (base == 'C') return 1;
	else if (base == 'G') return 2;
	else if (base == 'T' || base == 'U') return 3;
	else if (base == 'R') return 4;
	else if (base == 'Y') return 5;
	else if (base == 'K') return 6;
	else if (base == 'M') return 7;
	else if (base == 'S') return 8;
	else if (base == 'W') return 9;
	else if (base == 'B') return 10;
	else if (base == 'D') return 11;
	else if (base == 'H') return 12;
	else if (base == 'V') return 13;
	else if (base == 'N') return 14;

	return i;
}

void decodeNucleotide(char base, nuclVector &nuclList)
{
	int i = 0;
	char baseStr[] = "ACGTN";
	base = toupper(base);
	int pairMatrix[nuclNum][baseNum] =
	{
		/* A C G T N*/
		{1, 0, 0, 0, 0}, /*A*/
		{0, 1, 0, 0, 0}, /*C*/
		{0, 0, 1, 0, 0}, /*G*/
		{0, 0, 0, 1, 0}, /*T*/
		{1, 0, 1, 0, 0}, /*R*/
		{0, 1, 0, 1, 0}, /*Y*/
		{0, 0, 1, 1, 0}, /*K*/
		{1, 1, 0, 0, 0}, /*M*/
		{0, 1, 1, 0, 0}, /*S*/
		{1, 0, 0, 1, 0}, /*W*/
		{0, 1, 1, 1, 0}, /*B*/
		{1, 0, 1, 1, 0}, /*D*/
		{1, 1, 0, 1, 0}, /*H*/
		{1, 1, 1, 0, 0}, /*V*/
		{1, 1, 1, 1, 1}, /*N*/
	};
	int *baseMatrix = pairMatrix[encodeNucleotide(base)];
	for (i = 0; i < baseNum; i++)
	{
		if (baseMatrix[i] == 1)
		{
			char bp = baseStr[i];
			int cd = encodeNucleotide(bp);
			nuclList.push_back(cd);
		}
	}
}

int searchMotif(char *motif, char *seq)
{
	int pos = -1;
	if (motif == NULL) return pos;

	int i = 0;
	int j = 0;
	int k = 0;
	int seqLen   = strlen(seq);
	int motifLen = strlen(motif);
	int midLen   = int(seqLen / 2);

	for (i = 0; i < seqLen - motifLen; i++)
	{
		int motifScore = 0;
		for (j = 0; j < motifLen; j++)
		{
			int seqCode = encodeNucleotide(seq[i + j]);
			nuclVector nuclList;
			decodeNucleotide(motif[j], nuclList);
			for (k = 0; k < nuclList.size(); k++)
			{
				int mtfCode = nuclList[k];
				if (seqCode == mtfCode)
				{
					motifScore++;
					break;
				}
				//fprintf(stderr, "%c\t%d\n", motif[j], mtfCode);
			}
		}
		if (motifScore == motifLen)
		{
			pos = i;
			if (abs(pos - midLen) < 5) // nearby to variation site
				break;
		}
	}
	return pos;
}

void freeSite(RBSite * rbs)
{
	//safeFree(rbs->chrom);
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

// exp(-0.5 ((x-m)/s))**2 )/sqrt(2 pi s**2).
//double getGaussKD(double rh)
//{
//	double denominator = 1.0 / (double)sqrt(2.0 * M_PI);
//	double numerator = exp((-0.5) * pow(rh, 2));
//	double result = numerator * denominator;
//	return result;
//}

//double *computeKDEs(parameterInfo * paraInfo, double *data, int length)
//{
//	int i = 0;
//	double bandwidth = (double)paraInfo->bandwidth;
//	double *kdes  = (double *) safeMalloc(sizeof(double) * length);
//	int extendNum = paraInfo->bandwidth * 4;
//	double *kds   = (double *) safeMalloc(sizeof(double) * extendNum);
//
//	for (i = 0; i < extendNum; i++) kds[i] = getGaussKD(double(i) / bandwidth);
//
//	for (int j = 0; j < length; ++j)
//	{
//		double kde = 0.0;
//		int start = j - extendNum;
//		if (start < 0) start = 0;
//		int end = j + extendNum;
//		if (end > length) end = length;
//		for (i = start; i < end; i++)
//		{
//			int idx = j - i;
//			kde += data[i] * kds[idx];
//		}
//		kdes[j] = kde / bandwidth;
//	}
//	safeFree(kds);
//	return kdes;
//}

