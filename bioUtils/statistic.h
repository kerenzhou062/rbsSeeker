/* statistic head file */

#ifndef STATISTIC_H
#define STATISTIC_H

double log10Val(double pval);

double log10Lnp(double lnP);

double log10Pval(double pval);

double simpson(int overlapRegionLen, int queryRegionLen, int sampleRegionLen);

double jaccard(int overlapRegionLen, int queryRegionLen, int sampleRegionLen);

#endif /* End STATISTIC_H */