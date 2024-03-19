#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include<assert.h>
#include<math.h>
#include<stdint.h>
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
#include "varFile.h"

void getOneVariation(char *chrom, char *chromSeq, char *seq, char *cigar, Variation *var)
{
	int cigarLen         = strlen(cigar);
	uint32_t chromStart  = var->chromStart;
	uint32_t chromEnd    = var->chromEnd;
	char strand          = var->strand;
	variationVector varList;
	var->mutNum     = 0;
	var->mutation   = NULL;
	var->position   = NULL;
	int start      = 0;
	int end        = 0;
	int i          = 0;
	int j          = 0;
	int readLen    = 0;
	int alignLen   = 0;
	int profileLen = 0;

	uint32_t mutPos  = 0;
	uint8_t  mutCode = 0;
	varPair vp;

	for (i = 0; i < cigarLen; i++)
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
					mutPos  = chromStart + profileLen + j;
					mutCode = encodeVariation('D', 'D');
					vp      = make_pair(mutPos, mutCode);
					varList.push_back(vp);
				}
				alignLen   += len;
				profileLen += len;
				break;
			case 'I' :
				mutPos  = chromStart + profileLen;
				mutCode = encodeVariation('I', 'I');
				vp      = make_pair(mutPos, mutCode);
				varList.push_back(vp);
				alignLen += len;
				readLen  += len;
				break;
			case 'M' :
			case '=' :
				for (j = 0; j < len; j++)
				{
					mutPos       = chromStart + profileLen + j;
					char chrBase = chromSeq[mutPos];
					char seqBase = seq[readLen + j];
					chrBase      = revComBase(chrBase, strand);
					seqBase      = revComBase(seqBase, strand);
					if ( chrBase != seqBase && chrBase != 'N' && seqBase != 'N')
					{
						mutCode = encodeVariation(chrBase, seqBase);
						vp      = make_pair(mutPos, mutCode);
						varList.push_back(vp);
					}
				}
				alignLen   += len;
				readLen    += len;
				profileLen += len;
				break;
			case 'N' :
				// intron start
				mutPos  = chromStart + profileLen;
				mutCode = encodeVariation('N', 'N');
				vp      = make_pair(mutPos, mutCode);
				varList.push_back(vp);
				profileLen += len; // plus intron length
				// intron end
				mutPos  = chromStart + profileLen;
				mutCode = encodeVariation('N', 'N');
				vp      = make_pair(mutPos, mutCode);
				varList.push_back(vp);
				break;
			case 'S' :
				readLen += len;
				break;
			}
			cigar[i] = c;
			start = i + 1;
		}
	}
	int seqMutNum = varList.size();

	if (seqMutNum > 0)
	{
		var->mutNum   = seqMutNum;
		var->position = (uint32_t *) safeMalloc(sizeof(uint32_t) * seqMutNum);
		var->mutation = (uint8_t *)  safeMalloc(sizeof(uint8_t) * seqMutNum);
		for (i = 0; i < seqMutNum; i++)
		{
			var->position[i] = varList[i].first;
			var->mutation[i] = varList[i].second;
		}
	}
	varList.clear();
}

void freeChromVarMap(chromVarMap &varHash) /*free var map */
{
  for (chromVarMap::iterator mapItr = varHash.begin(); mapItr != varHash.end(); mapItr++) {
    varVector varList = mapItr->second;
    freeVarVector(varList);
  }
  varHash.clear();
}

void freeVarVector(varVector &varList)
{
  for (varVector::iterator vecItr = varList.begin(); vecItr != varList.end(); vecItr++)
  {
    Variation *var = *vecItr;
    freeVariationItem(var);
  }
  varList.clear();
}

void freeVariationItem(Variation *var)
{
	safeFree(var->chrom);
	if (var->position != NULL) safeFree(var->position);
	if (var->mutation != NULL) safeFree(var->mutation);
	safeFree(var);
}

char revComBase(char ch, char strand)
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

uint8_t encodeVariation(char a, char b)
{
	uint8_t code = 0;
	a = toupper(a);
	b = toupper(b);

	// matches
	if (a == 'A' && b == 'A')
	{
		code = 0;
	}
	else if (a == 'C' && b == 'C')
	{
		code = 1;
	}
	else if (a == 'G' && b == 'G')
	{
		code = 2;
	}
	else if (a == 'T' && b == 'T')
	{
		code = 3;
	}
	// mismatches-A
	else if (a == 'A' && b == 'C')
	{
		code = 4;
	}
	else if (a == 'A' && b == 'G')
	{
		code = 5;
	}
	else if (a == 'A' && b == 'T')
	{
		code = 6;
	}
	// mismatches-C
	else if (a == 'C' && b == 'A')
	{
		code = 7;
	}
	else if (a == 'C' && b == 'G')
	{
		code = 8;
	}
	else if (a == 'C' && b == 'T')
	{
		code = 9;
	}
	// mismatches-G
	else if (a == 'G' && b == 'A')
	{
		code = 10;
	}
	else if (a == 'G' && b == 'C')
	{
		code = 11;
	}
	else if (a == 'G' && b == 'T')
	{
		code = 12;
	}
	// mismatches-T
	else if (a == 'T' && b == 'A')
	{
		code = 13;
	}
	else if (a == 'T' && b == 'C')
	{
		code = 14;
	}
	else if (a == 'T' && b == 'G')
	{
		code = 15;
	}
	// deletions
	else if (a == 'D' && b == 'D')
	{
		code = 16;
	}
	// insertions
	else if (a == 'I' && b == 'I')
	{
		code = 17;
	}
	// introns
	else if (a == 'N' && b == 'N')
	{
		code = 18;
	}
	// start/truncation
	//else if (a == 'S' && b == 'S')
	//{
	//	code = 19;
	//}
	return code;
}

bool compareVar(const Variation *oVar, const Variation *tVar)
{
  return (oVar->chromStart < tVar->chromStart);
}
