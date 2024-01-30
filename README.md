# Curvibacter Transcriptomics

This repository stores all scripts necessary to perform a transcriptome analysis with reads from *Curvibacter* sp. AEP1-3.
RNA extraction and sequencing was performed similar to [Pietschke et. al. 2017](https://www.pnas.org/doi/10.1073/pnas.1706879114).

Trimming and mapping was performed with customized Python 3 scripts, which are located in the scripts directory of this repository.
The mapping procedure was performed on the HPC-System Hilbert of the Heinrich-Heine-University Düsseldorf.

## Quality Control and Trimming of RAW reads

RAW reads have been examined for quality using the FASTQC software (Version: v0.11.9). Following the quality assessment, RAW reads were trimmed using the Trimmomatic software (Version: 0.39). For detailed program settings, please refer to the script: trimming_curvibacter.py.

## Deduplication of trimmed reads

Trimmed reads were subsequently employed as input for a deduplication procedure conducted by BBMap (Version: 38.96). The deduplication process was executed using the dedupe.sh script. Following deduplication, the reads were reformatted using the reformat.sh script. For detailed program settings, kindly consult the script: deduplication_curvibacter.py.
Deduplication was performed due to the warning raised by FASTQC and to remove PCR-artefacts. The trimmed and deduplicated reads were then used as input for a second quality check with FASTQC.

## Kallisto read mapping procedure

The trimmed reads, as well as the trimmed and deduplicated reads, were mapped against the RefSeq genome (GCF_002163715.1) of *Curvibacter* sp. AEP1-3. The mapping procedure was executed using the Kallisto software (Version: 0.45.0).
Kallisto utilizes an innovative approach known as pseudoalignment to rapidly determine the compatibility of RNA-seq reads with potential transcripts, rather than aligning each read to a reference genome. This approach significantly speeds up the quantification process while maintaining accuracy.

Before mapping with Kallisto we need to generate an index file, this was done by using following command:

`kallisto index -i curvibacter_kallisto_index curvibacter_cds_genome.fna`

For detailed settings of the mapping procedure with Kallisto, please refer to the script: kallisto_procedure_curvibacter.py and/or kallisto_procedure_curvibacter_deduplicated_reads.py.

## STAR mapping procedure and RSEM read abundance measurements

In addition to the mapping procedure performed with Kallisto, trimmed and deduplicated reads were aligned using the STAR aligner (Version: 2.7.3a), and the resulting alignments were quantified using RSEM (Version: 1.2.31). The mapping and read abundance estimations were conducted using the rsem-calculate-expression command within the RSEM package.

The STAR-RSEM approach refers to the integration of two bioinformatics tools, STAR (Spliced Transcripts Alignment to a Reference) and RSEM (RNA-Seq by Expectation-Maximization), for the analysis of RNA-Seq data.

Here's a breakdown of the approach:

- STAR (Spliced Transcripts Alignment to a Reference): STAR is a fast and highly accurate RNA-seq read aligner. It is used to align RNA-seq reads to a reference genome, allowing the identification of potential splice junctions and accurate mapping of reads spanning these junctions. STAR is known for its efficiency in handling large-scale RNA-seq datasets.
- RSEM (RNA-Seq by Expectation-Maximization): RSEM is a software package designed for quantifying gene and isoform abundances from RNA-seq data. It uses the expectation-maximization (EM) algorithm to estimate expression levels, and it accounts for the uncertainty in read assignments to different isoforms and genes.

The rsem-calculate-expression command is a part of the RSEM package. This specific command is used to perform the quantification step, taking aligned reads (often generated by STAR) as input and estimating transcript or gene expression levels. It calculates the expected read counts and fragments per kilobase of transcript per million mapped reads (FPKM) values.

### References
- [Pietschke et. al. 2017](https://www.pnas.org/doi/10.1073/pnas.1706879114)
- [Kallisto](https://www.nature.com/articles/nbt.3519)

# TODO's
- [ ] add deseq2 scripts and description for differential gene expression analysis
- [ ] add clusterProfiler scripts and description for pathway analysis with DEGs
- [ ] add References for each software
