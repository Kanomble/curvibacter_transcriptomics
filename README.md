# Curvibacter Transcriptomics

## Sample Description

<div style="text-align: justify;white-space: pre-line;">
This repository stores all scripts necessary to perform a transcriptome analysis with reads from <i>Curvibacter</i> sp. AEP1-3.
RNA extraction and sequencing was performed similar to <a href="https://www.pnas.org/doi/10.1073/pnas.1706879114"> Pietschke et. al. 2017 </a>.


Trimming and mapping was performed with customized Python 3 scripts, which are located in the scripts directory of this
repository.
The mapping procedure was performed on the HPC-System Hilbert of the Heinrich-Heine-University Düsseldorf.

Differential gene expression analysis was performed with DESeq2 and custom R scripts on a local desktop computer. For
the analysis of the metaorganism, DESeq2 was used to compare normalized read counts from <i>Curvibacter</i> cells in
liquid R2A media with those from <i>Curvibacter</i> cells of the host organism <i>Hydra vulgaris</i>.
Samples of free living <i>Curvibacter</i> in R2A are labelled with a "G" as prefix, whereas samples with reads obtained
from the metaorganism are labelled with "Hydra" as prefix. Both prefixes are followed by a number indicating the
biological replicate.
For the analysis of mono-colonizing <i>Curvibacter</i> cells on <i>Hydra vulgaris</i>, please refer to the sample
description in BioProject <a href="https://www.ncbi.nlm.nih.gov/bioproject/PRJNA887579">PRJNA887579</a>
</div>

## From RAW-Reads to Read Counts

### Quality Control and Trimming of RAW reads

<div style="text-align: justify">
RAW reads have been examined for quality using the FastQC software (Version: v0.11.9). Following the quality assessment, RAW reads were trimmed using the Trimmomatic software (Version: 0.39). For detailed program settings, please refer to the script: trimming_curvibacter.py.
</div>

### Deduplication of trimmed reads

<div style="text-align: justify;white-space: pre-line;">
Trimmed reads were subsequently employed as input for a deduplication procedure conducted by BBMap (Version: 38.96). The deduplication process was executed using the dedupe.sh script. Following deduplication, the reads were reformatted using the reformat.sh script. For detailed program settings, refer to the script: deduplication_curvibacter.py.
Deduplication was performed due to the warning raised by FastQC and to remove PCR-artefacts. The trimmed and deduplicated reads were then used as input for a second quality check with FastQC.
</div>

### Kallisto mapping procedure and read abundance measurements

<div style="text-align: justify;white-space: pre-line;">
The trimmed reads, as well as the trimmed and deduplicated reads, were mapped against the RefSeq genome (GCF_002163715.1) of <i>Curvibacter</i> sp. AEP1-3. The mapping procedure was executed using the Kallisto software (Version: 0.45.0).
Kallisto utilizes an innovative approach known as pseudoalignment to rapidly determine the compatibility of RNA-seq reads with potential transcripts, rather than aligning each read to a reference genome. This approach significantly speeds up the quantification process while maintaining accuracy.

Before mapping with Kallisto we need to generate an index file, this was done by using following command:
</div>

`kallisto index -i curvibacter_kallisto_index curvibacter_cds_genome.fna`

<div style="text-align: justify;white-space: pre-line;">
For detailed settings of the mapping procedure with Kallisto, please refer to the script: kallisto_procedure_curvibacter.py and/or kallisto_procedure_curvibacter_deduplicated_reads.py.
</div>

### STAR mapping procedure and RSEM read abundance measurements

<div style="text-align: justify;white-space: pre-line;">
In addition to the mapping procedure performed with Kallisto, trimmed and deduplicated reads were aligned using the STAR aligner (Version: 2.7.3a), and the resulting alignments were quantified using RSEM (Version: 1.2.31). The mapping and read abundance estimations were conducted using the rsem-calculate-expression command within the RSEM package.

The STAR-RSEM approach refers to the integration of two bioinformatics tools, STAR (Spliced Transcripts Alignment to a
Reference) and RSEM (RNA-Seq by Expectation-Maximization), for the analysis of RNA-Seq data.

Here's a breakdown of the approach:
</div>

- STAR (Spliced Transcripts Alignment to a Reference): STAR is a fast and highly accurate RNA-seq read aligner. It is
  used to align RNA-seq reads to a reference genome, allowing the identification of potential splice junctions and
  accurate mapping of reads spanning these junctions. STAR is known for its efficiency in handling large-scale RNA-seq
  datasets.
- RSEM (RNA-Seq by Expectation-Maximization): RSEM is a software package designed for quantifying gene and isoform
  abundances from RNA-seq data. It uses the expectation-maximization (EM) algorithm to estimate expression levels, and
  it accounts for the uncertainty in read assignments to different isoforms and genes.

<div style="text-align: justify;white-space: pre-line;">

The rsem-calculate-expression command is a part of the RSEM package. This specific command is used to perform the
quantification step, taking aligned reads (often generated by STAR) as input and estimating transcript or gene
expression levels. It calculates the expected read counts and fragments per kilobase of transcript per million mapped
reads (FPKM) values.
</div>

## Differential gene expression analysis with deseq2

<div style="text-align: justify;white-space: pre-line;">
To estimate differentially expressed genes between <i>Curvibacter</i> cells in liquid R2A media and <i>Curvibacter</i> cells living on the host <i>Hydra vulgaris</i>, DESeq2 was used.
The samples designated as HydraX consist of reads derived from the meta-organism <i>Hydra vulgaris</i>. This implies that within the initial raw read samples, genes from other bacterial colonizers are also present. Mapping those raw reads against the <i>Curvibacter</i> genome results in low read abundances.

However, DESeq2 performs a normalization procedure to account for differences in sequencing depth between samples. This
normalization ensures that the comparisons between samples are not biased by these differences. In addition, DESeq2
estimates the dispersion of counts. Dispersion measures the degree of variability in expression for each gene across
samples. Accurate estimation of dispersion is important for determining statistical significance.
DESeq2 uses a negative binomial distribution model to account for the inherent over-dispersion observed in count data.
Unlike a simple Poisson distribution, which assumes that the variance equals the mean, the negative binomial
distribution allows for greater flexibility in modeling the variance. Once normalization and dispersion estimation are
performed, DESeq2 conducts a statistical test (Wald test) to assess whether the observed differences in gene expression
between conditions are statistically significant. This test accounts for both biological variability and technical
variability in the data.
</div>

The DESeq2 analysis is implemented within two R scripts. One script was used for the Kallisto abundance files (abundance.tsv), the other for the RSEM genes.results files.
Due to the fact that *Curvibacter* is a bacterial species, isoforms.results files of RSEM have not been used in this analysis step.

## Gene Enrichment Analysis with AnnotationForge and clusterProfiler

With the <i>Curvibacter</i> proteome as input, BlastKOALA was used to identify KeggOntology (KO) identifier for each gene. The GenBank file of <i>Curvibacter</i> was used to extract all GeneOntology (GO) terms for each protein sequence.
Subsequently, the R package [AnnotationForge](https://bioconductor.org/packages/release/bioc/html/AnnotationForge.html) was used to construct an OrgDb package for <i>Curvibacter</i>. This OrgDb package was used as input for the enrichGO function of the [clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html) R package.

Differentially expressed genes have been filtered for significance (adjusted p-value "padj" <= 0.05). Downregulated genes with an log2FoldChange value of <= -1 and upregulated genes with an log2FoldChange value of >= 1 have been used as input for the enrichKEGG and enrichGO functions of the clusterProfiler R package.

## Gene Enrichment Analysis with GOATOOLS.

For the GO-term analysis with GOATOOLS, an association file with gene_id to GO identifier mapping was prepared using the <i>Curvibacter</i> genome GFF file obtained from NCBIs [FTP-Server](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/163/715/GCF_002163715.1_ASM216371v1/). The analysis was conducted on significantly up- and downregulated genes filtered by the adjusted p-value and the log2FoldChange parameters of the differential expression analysis tables. A custom python script was then used to plot the enriched gene sets.

# References

- [Pietschke et. al. 2017](https://www.pnas.org/doi/10.1073/pnas.1706879114) website
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) website
- [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) website
- [BBMap](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbmap-guide/) website
- [RSEM](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-323) paper
- [STAR](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-323) paper
- [Kallisto](https://www.nature.com/articles/nbt.3519) paper
- [DESeq2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8) paper
- [AnnotationForge](https://bioconductor.org/packages/release/bioc/html/AnnotationForge.html) website
- [clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html) website
- [GOATOOLS](https://www.nature.com/articles/s41598-018-28948-z)
