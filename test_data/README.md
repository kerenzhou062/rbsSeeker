# run rbsSeeker on testing CLIP-seq datasets
Here, we describe how to run rbsSeeker on the testing datasets

- [download](#Downloading)
- [run rbsSeeker](#Run)
- [Output](#Output)

# Downloading
* Download the testing datasets to your local folder.

    ```bash
    # suppose your path to run rbsSeeker is ./rbsSeeker_test
    mkdir ./rbsSeeker_test
    cd ./rbsSeeker_test
    wget 'https://rnasysu.com/encori/software_test/rbsSeeker/IP/ENCSR987FTF-rep1.STAR.Aligned.sortedByCoord.out.bam'
    wget 'https://rnasysu.com/encori/software_test/rbsSeeker/input/ENCSR799EKA-rep1.STAR.Aligned.sortedByCoord.out.bam'
    ```

# Run
* It's very easy to run rbsSeeker on testing dataset. Supposed `rbsSeeker` has been added to you `PATH` environment and the `hg38.fa` and `hg38.fa.fai` have been generated.
    ```bash
    rbsSeeker -p 0.05 -q 0.05 --rm -S 0.99 -d 5 -N -L 2 --skip --rnafold --fa hg38.fa --fai hg38.fa.fai --treat ./ENCSR987FTF-rep1.STAR.Aligned.sortedByCoord.out.bam --control ./ENCSR799EKA-rep1.STAR.Aligned.sortedByCoord.out.bam --outdir ./ --prefix ENCSR987FTF-rep1
    ```

# Output
* The results should be the same with the following result files.

    ```bash
    # suppose your path to run rbsSeeker is ./rbsSeeker_test
    mkdir validate_results
    cd validate_results
    wget 'https://rnasysu.com/encori/software_test/rbsSeeker/rbsSeeker_output/ENCSR987FTF-repl_Insertion.bed'
    wget 'https://rnasysu.com/encori/software_test/rbsSeeker/rbsSeeker_output/ENCSR987FTF-repl_Peak.bed'
    wget 'https://rnasysu.com/encori/software_test/rbsSeeker/rbsSeeker_output/ENCSR987FTF-repl_Truncation.bed'
    wget 'https://rnasysu.com/encori/software_test/rbsSeeker/rbsSeeker_output/ENCSR987FTF-repl_Mutation.bed'
    wget 'https://rnasysu.com/encori/software_test/rbsSeeker/rbsSeeker_output/ENCSR987FTF-repl_End.bed'
    wget 'https://rnasysu.com/encori/software_test/rbsSeeker/rbsSeeker_output/ENCSR987FTF-repl_Deletion.bed'
    ```

Here's the description of columns in the outputs:

| Column name      | Description
| -----------      |----------
| `chrom`          | chromosome name
| `chromStart`     | Start genomic coordiate of the event (e.g. peak, deletion, truncation or mutation) (0-base)
| `chromEnd`       | End genomic coordinate of the event
| `name`           | Unique event ID
| `score`          | RPM of peak/site height
| `strand`         | Genomic sense (+) or antisense (-) strand of the event
| `extendSeq`      | Sequence extended ±10 bp from the site or peak center
| `motifPos`       | Start position of the motif in the `extendSeq` column (0-based). `-1` indicates no motif found.
| `type`           | Type of peak or site
| `log10(p-value)` | Log10-transformed p-value
| `log10(q-value)` | Log10-transformed q-value
| `readNum`        | This column represents the number of reads that support a specific binding event.
| `height`         | This value represents the total number of reads covering a specific nucleotide position, regardless of whether they contain a variation event or not.
| `heightRpm`      | This is the `height` normalized to Reads Per Million (RPM) mapped reads.
| `mfold`          | This column represents the fold-enrichment of a binding signal compared to a background model.
| `ratio`          | this is the proportion of reads covering a site that also carry the specific variation.

# Acknowledgements
Thanks a lot to everyone who contributed to the public codes and libraries (e.g. BamTools) used by rbsSeeker.

# Contact
* Jian-Hua Yang <yangjh7@mail.sysu.edu.cn>, RNA Information Center, School of Life Sciences, Sun Yat-Sen University<BR>
* Keren Zhou <kzhou@stjude.org>, Department of Pathology, St. Jude Children’s Research Hospital, Memphis, TN, USA<BR>
