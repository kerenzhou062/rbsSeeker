# rbsSeeker
Here, we develop RNA-Binding Site Seeker (rbsSeeker), a new and unified software with Poisson and Hypergeometric modeling, to identify significant RBP-binding events from various kinds of CLIP-seq data, such as HITS-CLIP, PAR-CLIP, iCLIP and eCLIP data. rbsSeeker is also suitable for miCLIP data analysis.

rbsSeeker can identify highly convinced peaks and individual cross-linking sites, of which includes deletions, truncations and mutations.

- [rbsSeeker](#rbsSeeker)
- [installation](#installation)
- [input](#input)
- [output](#output)
- [usage](#usage)

# installation
* It's very easy to install rbsSeeker on a linux server with following commands:
```bash
# suppose your path to install software is /username/software

cd /username/software

git clone https://github.com/kerenzhou062/rbsSeeker.git

cd ./rbsSeeker

sh install.sh

export PATH=$PATH:/username/software/rbsSeeker/bin
```

# input
* A genome sequences in fasta format is required by `--fa` argument, which chromsome's name must be the same with input `bam` file.

* An fai index file is required by `--fai` argument, which can be generated by using `samtools`:

```bash
# suppose your input genome fasta file is hg38.fa

samtools faidx hg38.fa
```

* A `bam` file contains sequence alignment data is required by `--bam` argument, which can be generated by alignment tools, such as `STAR`, `HISAT2` and `Tophat`.
Pleae note that rbsSeeker can only accept bam from single-end reads as input for peak/site calling.

# output
rbsSeeker may have 8 following output files in bed format (0-base) depends on the input arguments and your dataset.
Output file | Description
-----------|----------
`peaks` | 
`peaks height` | 
`mutation sites` | 
`specific mutation sites` | 
`truncation sites` | 
`insertion sites` | 
`deletion sites` | 
`end sites` | 


Here's column descriptions in the outputs:

Column name | Description
-----------|----------
`chrom` | chromosome name
`chromStart` | start genomic coordiate of specific event (e.g. peak, deletion, truncation or mutation) (0-base)
`chromEnd` | end genomic coordiate of specific event
`name` | uniq event id
`score` | RPM of height
`strand` | genomic sense (+) or antisense (-) strand of specific event
`extendSeq` | sequence extended ± 10bp from the individual sites or the peak center
`motifPos` | start position of input motif in `extendSeq` column. Starts with `0`. `-1` means there is no motif found in `extendSeq`.
`type` | indicate the type of peak or site
`log10(p-value)` | log10 of p-value.
`log10(q-value)` | log10 of q-value.
`readNum` | read number.
`height` | maximum peak or site height.
`heightRpm` | maximum peak or site height in Reads Per Million (RPM).
`mfold` | fold enrichment ( maximum height / average coverage).
`ratio` | ratio of specific event ( read number / total read number).

# usage
The available options of rbsSeeker are as follow:

```shell
Usage: rbsSeeker [options] --fa <genome file> --fai <genome fai> --bam <mapped alignments>
[options]
  -v, --verbose                   : verbose information
  -V, --version                   : rbsSeeker version
  -h, --help                      : help informations
  -R, --PCR                       : remove pcr duplictions [default is not removed]
  -e, --rm                        : remove the muations in start or end sites [default is not removed]
  --fa <string>                   : genome file with FASTA format
  --fai <string>                  : genome fai file with FAI format
  --bam <string>                  : alignments file with BAM format
  -o, --outdir <string>           : output dir
  -P, --prefix <string>           : prefix for output files
  -s, --transcriptome <int>       : transcriptome size [e.g. in human, default=129600000]
  -T, --cvs <string>              : conversion string [e.g. TC in PAR-CLIP, CT in miCLIP]
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
  -u, --brc-len<int>              : barcode length. extend the barcode length for mispriming [default=0]
  -s, --min-ratio<double>         : minimum ratio for variation [default>=0]
  -S, --max-ratio<double>         : maximum ratio for variation [default<=1.0]
```


