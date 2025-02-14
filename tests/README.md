# run rbsSeeker on testing CLIP-seq datasets
Here, we describe how to run rbsSeeker on the testing datasets

- [download](#Downloading)
- [run rbsSeeker](#Run)
- [Output](#Output)

# Downloading
Download the testing datasets to your local folder.

    ```bash
    # suppose your path to run rbsSeeker is ./rbsSeeker_test
    mkdir ./rbsSeeker_test
    cd ./rbsSeeker_test
    wget 'https://rnasysu.com/encori/software_test/rbsSeeker.test.bam'
    wget 'https://rnasysu.com/encori/software_test/rbsSeeker.test.bam.bai'
    ```

# Run
* It's very easy to run rbsSeeker on testing dataset. Supposed `rbsSeeker` has been added to you `PATH` environment and the `hg38.fa` and `hg38.fa.fai` have been generated.
    ```bash
    rbsSeeker --fa hg38.fa --fai hg38.fa.fai --bam rbsSeeker.test.bam --prefix rbsSeeker --outdir ./
    ```

# Output
rbsSeeker may have 8 following output files in `bed format (0-base)` depends on the input arguments and your dataset.
| Output file               | Description
| -----------               |----------
| `*_Peak.bed`              | identified peak-calling results
| `*_Mutation.bed`          | identified mutation sites
| `*_TC.bed` or `*_GA.bed`  | identified `T-to-C` mutation sites when set `--cvs TC`
| `*_Truncation.bed`        | identified truncation sites
| `*_Deletion.bed`          | identified deletion sites
| `*_Insertion.bed`         | identified insertion sites
| `*_End.bed`               | identified ending sites


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

# Acknowledgements
Thanks a lot to everyone who contributed to the public codes and libraries (e.g. BamTools) used by rbsSeeker.

# Contact
* Jian-Hua Yang <yangjh7@mail.sysu.edu.cn>, RNA Information Center, School of Life Sciences, Sun Yat-Sen University<BR>
* Keren Zhou <kzhou@coh.org>, Department of Systems Biology, Beckman Research Institute of City of Hope, Monrovia, CA, USA<BR>
