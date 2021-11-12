# rbsSeeker
Here, we describe the RNA-Binding Site Seeker (rbsSeeker), a new and unified software with Poisson and Hypergeometric modeling, to identify significant RBP-binding sites (`peaks` and `individual cross-linking sites`) from various kinds of CLIP-seq data, such as `HITS-CLIP`, `PAR-CLIP`, `iCLIP` and `eCLIP` data. rbsSeeker is also suitable for `miCLIP` data analysis.

rbsSeeker can identify highly convinced peaks and individual cross-linking sites, of which includes deletions, truncations and mutations.

- [rbsSeeker](#rbsSeeker)
- [Installation](#Installation)
- [Input](#Input)
- [Output](#Output)
- [Example](#Example)
- [Basic usage](#Basic-usage)
- [Acknowledgements](#Acknowledgements)
- [Contact](#Contact)

# Installation
* It's very easy to install rbsSeeker on a linux server with following commands:
```bash
# suppose your path to install software is /username/software

cd /username/software

git clone https://github.com/kerenzhou062/rbsSeeker.git

cd ./rbsSeeker

sh install.sh

export PATH=$PATH:/username/software/rbsSeeker/bin
```

# Input
* A genome sequences in fasta format is required by `--fa` argument, which chromsome's name must be the same with input `bam` file.

* An fai index file is required by `--fai` argument, which can be generated by using `samtools`:

    ```bash
    # suppose your input genome fasta file is hg38.fa
    
    samtools faidx hg38.fa
    ```

* A `bam` file contains sequence alignment data is required by `--bam` argument, which can be generated by alignment tools, such as `STAR`, `HISAT2` and `Tophat`.

* Pleae note that rbsSeeker can only accept bam from `single-end reads` as input for peak/site calling.

# Output
rbsSeeker may have 8 following output files in `bed format (0-base)` depends on the input arguments and your dataset.
| Output file           | Description
| -----------           |----------
| `*_Peak.bed`          | identified peak-calling results
| `*_PeakHeight.bed`    | identified peak-calling results
| `*_Mutation.bed`      | identified mutation sites
| `*_CT.bed`            | identified `C-to-T` mutation sites when set `--cvs CT`
| `*_Truncation.bed`    | identified truncation sites
| `*_Deletion.bed`      | identified deletion sites
| `*_Insertion.bed`     | identified insertion sites
| `*_End.bed`           | identified ending sites


Here's the description of columns in the outputs:

| Column name      | Description
| -----------      |----------
| `chrom`          | chromosome name
| `chromStart`     | start genomic coordiate of specific event (e.g. peak, deletion, truncation or mutation) (0-base)
| `chromEnd`       | end genomic coordiate of specific event
| `name`           | uniq event id
| `score`          | RPM of peak/site height
| `strand`         | genomic sense (+) or antisense (-) strand of specific event
| `extendSeq`      | sequence extended ± 10bp from the individual sites or the peak center
| `motifPos`       | start position of input motif in `extendSeq` column. Starts with `0`. `-1` means there is no motif found in `extendSeq`.
| `type`           | indicate the type of peak or site
| `log10(p-value)` | log10 of p-value.
| `log10(q-value)` | log10 of q-value.
| `readNum`        | read number of peak/site.
| `height`         | maximum peak/site height.
| `heightRpm`      | maximum peak/site height in Reads Per Million (RPM).
| `mfold`          | fold enrichment ( maximum height / average coverage).
| `ratio`          | ratio of specific event ( read number / total read number).

# Example
Here is an example that shows how to use `rbsSeeker` to identify N6-methyladenosine (`m6A`) sites at `single-base resolution` from `miCLIP` data.

* Supposed the raw reads from miCLIP data were properly processed (e.g. adapters trimmed, PCR duplicates removed, bacodes removed) and aligned to the propper genome (e.g. hg38). So then the reads alignment results (`miCLIP.sorted.bam`) were used for downstream analysis.

    * Run `rbsSeeker` on `miCLIP.sorted.bam`
        ```bash
        # rbsSeeker m6A sites calling, -t 129600000 is for human transcriptome
        # this step is usually finished within 20 minutes, it depends on the sizes of your datasets
        rbsSeeker -T CT -L 20 -t 129600000 -n 1 -H 3 -d 1 -p 0.05 -q 0.1 \
          -o ./output -P miCLIP --fa hg38.fa --fai hg38.fa.fai --bam miCLIP.sorted.bam > miCLIP.rbsSeeker.log
        ```

* Output files from `rbsSeeker` results<BR>
    * miCLIP_CT.bed
    * miCLIP_End.bed
    * miCLIP_Insertion.bed
    * miCLIP_Mutation.bed
    * miCLIP_Peak.bed
    * miCLIP_PeakHeight.bed
    * miCLIP_Truncation.bed

* Identify potential m6A sites (`DRACH motif`) (supposed the miCLIP was performed with [Abcam antibody](https://www.nature.com/articles/nmeth.3453/figures/1))

    * Substitutions at m6A site from `miCLIP_CT.bed`
        ```bash
        #substitutions at m6A site
        awk 'BEGIN{FS="\t";OFS="\t";}
        {
          if (FNR >1) {
            if ($8 == 8) {
              seq = substr($7, 9, 5);
              if (seq ~ /[AGT][AG]AC[ACT]/) {
                $4="mut|"seq"|"FNR;
                print $1, $2, $3, $4, $5, $6;
              }
            }
          }
        }' miCLIP_Mutation.bed > miCLIP.mut.bed
        ```

    * C->T mucations at +1 position of m6A site from `miCLIP_Mutation.bed`
        ```bash
        #C->T mucations at +1 position of m6A site
        awk 'BEGIN{FS="\t";OFS="\t";}
        {
          if (FNR >1) {
            if ($8 == 7) {
              if ($6 == "+") {
                $2 = $2 - 1;
                $3 = $2 + 1;
              }else{
                $3 = $3 + 1;
                $2 = $3 - 1;
              }
              seq = substr($7, 8, 5);
              if (seq ~ /[AGT][AG]AC[ACT]/) {
                $4="CT|"seq"|"FNR;
                print $1, $2, $3, $4, $5, $6;
              }
            }
          }
        }' miCLIP_CT.bed > miCLIP.CT.bed
        ```

    * Truncations at +2 position of m6A site from `miCLIP_Truncation.bed`
        ```bash
        #truncations at +2 position of m6A site
        awk 'BEGIN{FS="\t";OFS="\t";}
        {
          if (FNR >1) {
            if ($8 == 6) {
              if ($6 == "+") {
                $2 = $2 - 2;
                $3 = $2 + 1;
              }else{
                $3 = $3 + 2;
                $2 = $3 - 1;
              }
              seq = substr($7, 7, 5);
              if (seq ~ /[AGT][AG]AC[ACT]/) {
                $4="trunc|"seq"|"FNR;
                print $1, $2, $3, $4, $5, $6;
              }
            }
          }
        }' miCLIP_Truncation.bed > miCLIP.trunc.bed
        ```

* Pool identified sites together and get the final result `miCLIP.merge.bed` ([`bedtools`](https://bedtools.readthedocs.io/en/latest/content/bedtools-suite.html) is required)
    ```bash
    cat miCLIP.mut.bed miCLIP.CT.bed miCLIP.trunc.bed | awk 'BEGIN{FS="\t";OFS="\t";}
      {
        split($4,arr,"|");
        seq=arr[2];
        key=$1"\t"$2"\t"$3"\t"seq"\t"$6;
        readNumArr[key] += $5;
      }
      END{
        for (key in readNumArr) {
          split(key, arr,"\t");
          print arr[1], arr[2], arr[3], arr[4], readNumArr[key], arr[5];
        }
      }' | sort -k1,1 -k2,2n | bedtools intersect -a stdin \
      -b miCLIP.mut.bed miCLIP.CT.bed miCLIP.trunc.bed \
      -names mut CT trunc -s -wa -wb | \
      awk 'BEGIN{FS="\t";OFS="\t";}
      {
        key=$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6;
        if(key in hashArrA){
          hashArrA[key] += 1;
          hashArrB[key] = hashArrB[key]","$7;
        }else{
          hashArrA[key] = 1;
          hashArrB[key] = $7;
        }
      }
      END{
        for (key in hashArrA){
          print key, hashArrB[key], hashArrA[key];
        }
      }' | sort -k1,1 -k2,2n | \
      awk 'BEGIN{FS="\t";OFS="\t";}{$4=$4"|"FNR;print}'> miCLIP.merge.bed
    ```

# Basic Usage
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

# Acknowledgements
Thanks a lot to everyone who contributed to the public codes and libraries (e.g. BamTools) used by rbsSeeker.

# Contact
* Jian-Hua Yang <yangjh7@mail.sysu.edu.cn>, RNA Information Center, School of Life Sciences, Sun Yat-Sen University<BR>
* Keren Zhou <kzhou@coh.org>, Department of Systems Biology, Beckman Research Institute of City of Hope, Monrovia, CA, USA<BR>
