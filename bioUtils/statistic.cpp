#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include <getopt.h>
#include<assert.h>
#include<math.h>

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
#include "statistic.h"

double log10Val(double pval)
{
	double returnVal = 0;
	if (pval == 0)
	{
		returnVal = -325.00;
	}
	else {
		returnVal = log(pval) / log(10);
	}
	return returnVal;
}

double log10Lnp(double lnP)
{
	double returnVal = 0;
	returnVal = lnP / log(10);
	if (returnVal > 0) returnVal = 0;
	return returnVal;
}

double log10Pval(double pval)
{
	double returnVal = 0;
	if (pval == 0)
	{
		returnVal = -325.00;
	}
	else {
		returnVal = log(pval) / log(10);
	}
	if (returnVal > 0) returnVal = 0;
	return returnVal;
}

double simpson(int overlapRegionLen, int queryRegionLen, int sampleRegionLen) {
	double 	simpson = 0;
	int minVal = queryRegionLen;
	if (minVal > sampleRegionLen) minVal = sampleRegionLen;
	simpson = overlapRegionLen / (double)minVal;
	return simpson;
}

double jaccard(int overlapRegionLen, int queryRegionLen, int sampleRegionLen)
{
	double jaccardVal = 0;
	jaccardVal = overlapRegionLen / (double)(queryRegionLen + sampleRegionLen - overlapRegionLen);
	return jaccardVal;
}
