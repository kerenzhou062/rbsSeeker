# rbsSeeker
A new and unified software with Poisson and Hypergeometric modeling, to identify significant RBP-binding events from various kinds of CLIP-seq data.
- [rbsSeeker](#rbsSeeker)
- [installation](#installation)
- [usage](#usage)

# installation
* It's very easy to install rbsSeeker on a linux server.
```

git clone https://github.com/kerenzhou062/rbsSeeker.git

cd ./rbsSeeker

sh install.sh
```

# usage
```shell
Usage: rbsSeeker [options] --fa <genome file> --fai <genome fai> --bam <mapped alignments>
[options]
  -v, --verbose                   : verbose information
  -V, --version                   : rbsSeeker version
  -h, --help                      : help informations
  -R, --PCR                       : remove pcr duplictions[default is not removed]
  -e, --rm                        : remove the muations in start or end sites[default is not removed]
  --fa <string>                   : genome file with FASTA format
  --fai <string>                  : genome fai file with FAI format
  --bam <string>                  : alignments file with BAM format
  -o, --outdir <string>           : output dir
  -P, --prefix <string>           : prefix for output files
  -s, --transcriptome <int>       : transcriptome size[e.g. in human, default=129600000]
  -T, --cvs <string>              : conversion string[e.g. TC in PAR-CLIP, CT in miCLIP]
  -c, --min-peak-len <int>        : minimum length for a peak [default>=10]
  -i, --min-read-len <int>        : minimum read length [default>=10]
  -a, --max-read-len <int>        : maximum read length [default<=5000000]
  -n, --min-read-num <double>     : minimum number of reads for calling a peak [Default=1]
  -L, --max-locus-num <int>       : maximum locus number of reads for mapping to genome [Default=20]
  -H, --min-height <double>       : minimum read height for calling a site [Default=5]
  -r, --rpm <double>              : minimum rpm height for calling a site [Default=0.01]
  -d, --min-var <double>          : minimum read number for calling a variational site [Default=1]
  -p, --pval <double>             : minimum p value for calling a site [Default<=0.05]
  -q, --qval <double>             : minimum q value for calling a site [Default<=0.05]
  --primer<string>                : primer sequene for removing the mispriming [default=NULL]
  -u, --brc-len<int>              : barcode length. extend the barcode length for mispriming[default=0]
  -s, --min-ratio<double>         : minimum ratio for variation [default>=0]
  -S, --max-ratio<double>         : maximum ratio for variation [default<=1.0]
```
