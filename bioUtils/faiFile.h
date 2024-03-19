#ifndef FAI_HEAD_H
#define FAI_HEAD_H

struct faidxInfo {
  int lineLen;
  int lineBlen;
  long len;
  long offset;
};

typedef struct faidxInfo faidx;

typedef map<string, faidx *> faidxMap;

typedef map<string, char *> chromSeqMap;

typedef map<string, uint16_t> chromCodeMap;

typedef vector<string> chromCodeVector;

long readFai(FILE *fp, faidxMap &faiHash);

void freeFaiList(faidxMap &fai);

char *faidxFetchSeq(FILE *gfp, const faidx *fai, int start, int end, char strand);

void fetchAllSeq(FILE *gfp, faidxMap &faidxHash, chromSeqMap &chromSeqHash);

void safeFreeChromSeq(chromSeqMap &chromSeqHash);

uint16_t encodeChrom(faidxMap &faidxHash, chromCodeMap &chromCodeHash, chromCodeVector &chromCodeList);

#endif /* End FAI_HEAD_H */
