[![Travis-ci tests:](https://travis-ci.org/onuryukselen/rnaseq.svg?branch=master)](https://travis-ci.org/onuryukselen/rnaseq)

RNA-seq pipeline includes Quality Control, rRNA filtering, Genome Alignment using HISAT2, STAR and Tophat2, and estimating gene and isoform expression levels by RSEM and featureCounts.  
  
Steps:
  1. For Quality Control, we use FastQC to create qc outputs. There are optional read quality filtering (trimmomatic), read quality trimming (trimmomatic), adapter removal (cutadapt) processes available.
  2. Bowtie2/Bowtie/STAR is used to count or filter out common RNAs reads (eg. rRNA, miRNA, tRNA, piRNA etc.). 
  3. RSEM is used to align RNA-Seq reads to a reference transcripts and estimates gene and isoform expression levels.
  4. HISAT2, STAR and Tophat2 is used to align RNA-Seq reads to a genome. Optional estimation of gene and isoform expression levels could be done by featureCount.
  5. Genome-wide Bam analysis is done by RseQC, Picard.
  6. Optionally you can create Integrative Genomics Viewer (IGV)  and Genome Browser Files (TDF and Bigwig, respectively)

Program Versions:
  - FastQC v0.11.8
  - Star v2.6.1
  - Hisat2 v2.1.0
  - Picard v2.18.27
  - Rseqc v2.6.2
  - Samtools v1.3
  - Subread v1.6.4
  - Multiqc v1.7
  - Tophat v2.1.1
  - RSEM v1.3.1
  - Bowtie2 v2.3.5
  - Bowtie v1.2.2
  - Trimmomatic v0.39
  - Igvtools v2.5.3
  - Bedtools v2.27.1
  - Fastx_toolkit v0.0.14
  - Ucsc-wigToBigWig v366
  - Pdfbox-App v2.0.0
 


# Local installation:

To start using the dolphinnext/rnaseq pipeline, follow the steps in the following link:[Local installation](https://github.com/onuryukselen/rnaseq/blob/master/docs/local.md)





