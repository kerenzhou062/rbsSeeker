/********************************************************************
 * rbsSeeker: identify RBP binding sites at single-base resolution
 * Author: jianhua yang
 * Email: yangjh7@mail.sysu.edu.cn
 * Copyright: School of Life Sciences, Sun Yat-sen University
 * $ 2019/12/09
 ********************************************************************/
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include<getopt.h>
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
#include "statistic.h"
#include "rbsSeeker.h"

char version[] = "rbsSeeker version 0.1";
void usage(void);

int main(int argc, char *argv[])
{
  char *outdir        = NULL;
  char *prefix        = NULL;
  char *strLine       = NULL;
  char *faFile        = NULL;
  char *faiFile       = NULL;
  FILE *genomefp      = NULL;
  FILE *faifp         = NULL;
  char *treatFile     = NULL;
  char *controlFile   = NULL;
  const char *defpre  = "rbsSeeker";
  const char *defout  = "rbsSeekerResults";
  int showVersion     = 0;
  int showHelp        = 0;
  int i               = 0;
  int c               = 0;
  char createDir[255];
  char outputDir[255];

  parameterInfo *paraInfo = (parameterInfo *)safeMalloc(sizeof(parameterInfo));

  /* parse commmand line parameters */

  if (argc == 1)
  {
    usage();
  }

  const char *shortOptions = "vhVRUekNKo:n:i:a:c:H:d:t:p:T:r:q:l:L:f:F:b:m:P:r:M:u:s:S:B:C:D:w:";

  const struct option longOptions[] =
  {
    { "verbose" , no_argument , NULL, 'v' },
    { "help" , no_argument , NULL, 'h' },
    { "version" , no_argument , NULL, 'V' },
    { "PCR" , no_argument , NULL, 'R' },
    { "full" , no_argument , NULL, 'U' },
    { "norm" , no_argument , NULL, 'N' },
    { "skip" , no_argument , NULL, 'k' },
    { "rm" , no_argument , NULL, 'e' },
    { "rnafold" , no_argument , NULL, 'K' },
    { "outdir" , required_argument , NULL, 'o' },
    { "prefix" , required_argument , NULL, 'P' },
    { "min-read-num" , required_argument, NULL, 'n' },
    { "min-read-len" , required_argument, NULL, 'i' },
    { "max-read-len" , required_argument, NULL, 'a' },
    { "min-peak-len" , required_argument, NULL, 'c' },
    { "max-peak-len" , required_argument, NULL, 'C' },
    { "min-height" , required_argument, NULL, 'H' },
    { "min-var" , required_argument, NULL, 'd' },
    { "max-overlap-len" , required_argument, NULL, 'l' },
    { "max-locus-num" , required_argument, NULL, 'L' },
    { "transcriptome" , required_argument, NULL, 't' },
    { "rpm" , required_argument, NULL, 'r' },
    { "cvs" , required_argument, NULL, 'T' },
    { "pval" , required_argument, NULL, 'p' },
    { "qval" , required_argument, NULL, 'q' },
    { "mfold" , required_argument, NULL, 'm' },
    { "pfold" , required_argument, NULL, 'O' },
    { "fa" , required_argument, NULL, 'f' },
    { "fai" , required_argument, NULL, 'F' },
    { "treat" , required_argument, NULL, 'b' },
    { "control" , required_argument, NULL, 'B' },
    { "motif" , required_argument, NULL, 'M' },
    { "brc-len" , required_argument, NULL, 'u' },
    { "min-ratio" , required_argument, NULL, 's' },
    { "max-ratio" , required_argument, NULL, 'S' },
    { "drop" , required_argument, NULL, 'D' },
    { "window" , required_argument, NULL, 'w' },
    {NULL, 0, NULL, 0} ,  /* Required at end of array. */
  };

  //initial parameters
  paraInfo->maxReadDist   = 0; // the reads were clusted when the distance is less than 0 nt
  paraInfo->minClusterLen = 10;
  paraInfo->maxClusterLen = 1000000;
  paraInfo->minReadNum = 1;
  paraInfo->minReadLen = 18;
  paraInfo->maxReadLen = 1000000;
  paraInfo->minT2cNum  = 1;
  paraInfo->minHeight  = 5;
  paraInfo->minMutNum  = 1;
  paraInfo->verbose = 0;
  paraInfo->PCR = 0;
  paraInfo->bam = 0;
  paraInfo->genomeSize = 0;
  paraInfo->pval = 0.05;
  paraInfo->qval = 0.05;
  paraInfo->cvs  = NULL;
  paraInfo->maxLocusNum = 5;
  paraInfo->mfold = 2;
  paraInfo->pfold = 2;
  paraInfo->fullLength = 0;
  paraInfo->transcriptomeSize = 129600000; // refgene sizes in human 2018/07/15
  paraInfo->minRpm       = 1;
  paraInfo->barcodeLen   = 0;
  paraInfo->motifSeq     = NULL;
  paraInfo->rmSeMutation = 0;
  paraInfo->minRatio = 0.0001;
  paraInfo->maxRatio = 1.0;
  paraInfo->clipType = 0;
  paraInfo->totalTreatNum    = 1;
  paraInfo->totalCtrlNum     = 1;
  paraInfo->treatVsCtrlRatio = 1;
  paraInfo->cvsCode    = 0;
  paraInfo->norm       = 0;
  paraInfo->dropRatio  = 0.25;
  paraInfo->skipSplice = 0;
  paraInfo->windowLen  = 500;
  paraInfo->rnafold    = 0;

  while ((c = getopt_long(argc, argv, shortOptions, longOptions, NULL)) >= 0)
  {
    switch (c)
    {
    case 'v':
      paraInfo->verbose = 1;
      break;
    case 'h':
      showHelp = 1;
      break;
    case 'V':
      showVersion = 1;
      break;
    case 'R':
      paraInfo->PCR = 1;
      break;
    case 'e':
      paraInfo->rmSeMutation = 1;
      break;
    case 'N':
      paraInfo->norm = 1;
      break;
    case 'k':
      paraInfo->skipSplice = 1;
      break;
    case 'K':
      paraInfo->rnafold = 1;
      break;
    case 'P':
      prefix    = optarg;
      break;
    case 'o':
      outdir  = optarg;
      break;
    case 'f':
      faFile  = optarg;
      break;
    case 'F':
      faiFile = optarg;
      break;
    case 'b':
      treatFile   = optarg;
      break;
    case 'B':
      controlFile = optarg;
      break;
    case 'T':
      paraInfo->cvs  = strClone(optarg);
      break;
    case 'M':
      paraInfo->motifSeq = optarg;
      break;
    case 'n':
      paraInfo->minReadNum = atof(optarg);
      break;
    case 'd':
      paraInfo->minMutNum = atof(optarg);
      break;
    case 'r':
      paraInfo->minRpm = atof(optarg);
      break;
    case 'l':
      paraInfo->maxReadDist = atoi(optarg);
      break;
    case 'L':
      paraInfo->maxLocusNum = atoi(optarg);
      break;
    case 'i':
      paraInfo->minReadLen = atoi(optarg);
      break;
    case 'a':
      paraInfo->maxReadLen = atoi(optarg);
      break;
    case 't':
      paraInfo->transcriptomeSize = atof(optarg);
      break;
    case 'c':
      paraInfo->minClusterLen = atoi(optarg);
      break;
    case 'C':
      paraInfo->maxClusterLen = atoi(optarg);
      break;
    case 'H':
      paraInfo->minHeight = atof(optarg);
      break;
    case 'p':
      paraInfo->pval = atof(optarg);
      break;
    case 'm':
      paraInfo->mfold = atof(optarg);
      break;
    case 'O':
      paraInfo->pfold = atof(optarg);
      break;
    case 'q':
      paraInfo->qval = atof(optarg);
      break;
    case 'u':
      paraInfo->barcodeLen = atoi(optarg);
      break;
    case 's':
      paraInfo->minRatio = atof(optarg);
      break;
    case 'S':
      paraInfo->maxRatio = atof(optarg);
      break;
    case 'D':
      paraInfo->dropRatio = atof(optarg);
      break;
    case 'w':
      paraInfo->windowLen = atoi(optarg);
      break;
    case '?':
      showHelp = 1;
      break;
    default:
      usage();
    }
  }
  // help for version
  if (showVersion)
  {
    fprintf(stderr, "%s\n", version);
    exit(1);
  }

  if (showHelp)
  {
    usage();
    exit(1);
  }

  if (paraInfo->cvs != NULL)
  {
    if (strlen(paraInfo->cvs) != 2)
    {
      fprintf(stderr, "the cvs: %s string must be two characters. e.g. TC\n", paraInfo->cvs);
      exit(1);
    }
  }
  if (treatFile == NULL)
  {
    fprintf(stderr, "ERROR: please set the option: --treat <mapped alignments, BAM format>\n");
    usage();
  }

  if (faFile != NULL)
  {
    genomefp = (FILE *) fopen(faFile, "r");
    if (genomefp == NULL)
    {
      fprintf(stderr, "ERROR: Can't open genome file: %s\n", faFile);
      usage();
    }
  }
  else
  {
    fprintf(stderr, "ERROR: please set the option: --fa <genome file>\n");
    usage();
  }

  if (faiFile != NULL)
  {
    faifp = (FILE *) fopen(faiFile, "r");
    if (faifp == NULL)
    {
      fprintf(stderr, "ERROR: Can't open genome file: %s\n", faiFile);
      usage();
    }
  }
  else
  {
    fprintf(stderr, "ERROR: please set the option: --fai <fai file>\n");
    usage();
  }

  strcpy(createDir, "mkdir -p ");
  if (outdir != NULL)
  {
    strcpy(outputDir, outdir);
    strcat(createDir, outdir);
  }
  else
  {
    strcpy(outputDir, defout);
    strcat(createDir, defout);
  }

  strcat(outputDir, "/");
  if (prefix != NULL)
  {
    strcat(outputDir, prefix);
  }
  else
  {
    strcat(outputDir, defpre);
  }

  if (paraInfo->windowLen < 0) paraInfo->windowLen = 0;
  if (paraInfo->minMutNum < 1) paraInfo->minMutNum = 1;
  if (paraInfo->maxClusterLen < 50 ) paraInfo->maxClusterLen = 50;
  if (paraInfo->maxReadLen < 50 ) paraInfo->maxReadLen = 50;
  if (paraInfo->minHeight < 2) paraInfo->minHeight = 2;

  if (paraInfo->verbose)  fprintf(stderr, "#create dir \"%s\"\n", outputDir);
  system(createDir);
  if (paraInfo->pfold < 2.0) fprintf(stderr, "Warnning: we suggest the pfold option must be large than 2.0\n");
  if (paraInfo->mfold < 2.0) fprintf(stderr, "Warnning: we suggest the mfold option must be large than 2.0\n");
  if (paraInfo->pval > 0.05) fprintf(stderr, "Warnning: we suggest the pval option must be less than 0.05\n");
  if (paraInfo->qval > 0.05) fprintf(stderr, "Warnning: we suggest the qval option must be less than 0.05\n");

  double orgPval = paraInfo->pval;
  double orgQval = paraInfo->qval;
  paraInfo->pval = log10Val(paraInfo->pval);
  paraInfo->qval = log10Val(paraInfo->qval);

  if (paraInfo->verbose)  fprintf(stderr, "#log10(p-val:%.5f) cutoff is %g\tlog10(q-val:%.5f) cutoff is %g\n", orgPval, paraInfo->pval, orgQval, paraInfo->qval);
  if (paraInfo->cvs != NULL) paraInfo->cvsCode = encodeVariation(paraInfo->cvs[0], paraInfo->cvs[1]);
  seekRbsSites(genomefp, faifp, outputDir, paraInfo, treatFile, controlFile);
  fprintf(stderr, "program end\n");
  freeParameters(paraInfo);
  fclose(faifp);
  fclose(genomefp);

  return 0;
}

void usage(void)
{
  fprintf(stderr, "%s", "Usage: rbsSeeker [options] --fa <genome file> --fai <genome fai> --treat <mapped alignments>\n\
[options]\n\
-v/--verbose               : verbose information\n\
-V/--version               : rbsSeeker version\n\
-h/--help                  : help informations\n\
-R/--PCR                   : remove pcr duplictions[default is not removed]\n\
-e/--rm                    : remove the muations in start or end sites[default is not removed]\n\
-N/--norm                  : normalized the reads to locus[default is not normalized]\n\
-k/--skip                  : skip reads spanning the intron[default is not skip]\n\
-K/--rnafold               : RNAfold for your sequences[default is not]\n\
--fa <string>              : genome file with FASTA format\n\
--fai <string>             : genome fai file with FAI format\n\
--treat <string>           : alignments treatment file with BAM format\n\
--control <string>         : alignments control file with BAM format\n\
-o/--outdir <string>       : output dir\n\
-P/--prefix <string>       : prefix for output files\n\
-t/--transcriptome <int>   : transcriptome size[e.g. in human, default=129600000]\n\
-T/--cvs <string>          : conversion string[e.g. TC in PAR-CLIP, CT in miCLIP]\n\
-c/--min-peak-len <int>    : minimum length for a peak [default>=10]\n\
-C/--max-peak-len <int>    : maximum length for a peak [default<=1000000]\n\
-i/--min-read-len <int>    : minimum read length [default>=18]\n\
-a/--max-read-len <int>    : maximum read length [default<=1000000]\n\
-n/--min-read-num <double> : minimum number of reads for calling a peak [Default=1]\n\
-L/--max-locus-num <int>   : maximum locus number of reads for mapping to genome [Default=5]\n\
-H/--min-height <double>   : minimum read height for calling a site [Default=5]\n\
-r/--rpm <double>          : minimum rpm height for calling a site [Default=1]\n\
-d/--min-var <double>      : minimum read number for calling a variational site [Default=1]\n\
-p/--pval <double>         : minimum p value for calling a site [Default<=0.05]\n\
-q/--qval <double>         : minimum q value for calling a site [Default<=0.05]\n\
-M/--motif<string>         : search motif in the identified sequences  [default=NULL]\n\
-s/--min-ratio<double>     : minimum ratio for variation [default>=0]\n\
-S/--max-ratio<double>     : maximum ratio for variation [default<=1.0]\n\
-m/--mfold<double>         : minimum fold-change for variation[default>=2]\n\
-w/--window<int>           : window length for calculating the lamda of each peak [default=500]\n\
");
  exit(1);
}
