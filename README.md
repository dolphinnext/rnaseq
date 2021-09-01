[![Travis-ci tests:](https://travis-ci.org/dolphinnext/rnaseq.svg?branch=master)](https://travis-ci.org/dolphinnext/rnaseq) [![DOI:10.1101/689539](https://zenodo.org/badge/DOI/10.1101/689539.svg)](https://doi.org/10.1101/689539)

RNA-seq pipeline includes Quality Control, rRNA filtering, Genome Alignment using HISAT2, STAR and Tophat2, and estimating gene and isoform expression levels by RSEM and featureCounts. 
Alternatively, Kallisto could be used for quantifying abundances of transcripts based on pseudoalignments, without the need for alignment.
  
##### Steps:
  1. For Quality Control, we use FastQC to create qc outputs. There are optional read quality filtering (trimmomatic), read quality trimming (trimmomatic), adapter removal (cutadapt) processes available.
  2. Bowtie2/Bowtie/STAR is used to count or filter out common RNAs (eg. rRNA, miRNA, tRNA, piRNA etc.). 
  3. RSEM is used to align RNA-Seq reads to a reference transcripts and estimates gene and isoform expression levels. Alternatively, Kallisto used for quantifying abundances of transcripts based on pseudoalignments.
  4. HISAT2, STAR and Tophat2 are used to align RNA-Seq reads to a genome. Optionally, counting reads to genomic features such as genes, exons, promoters and genomic bins could be done by featureCounts.
  5. Genome-wide Bam analysis is done by RseQC, Picard.
  6. Optionally you can create Integrative Genomics Viewer (IGV)  and Genome Browser Files (TDF and Bigwig, respectively)

##### Pipeline Container:
  * Docker: dolphinnext/rnaseq:3.0

##### Program Versions:
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
  - Kallisto v0.46.0
  - Megadepth v1.1.0
 

##### Run through DolphinNext User Interface:

To start using the dolphinnext/rnaseq pipeline please go to [*DolphinNext Web page*](https://dolphinnext.umassmed.edu/index.php?np=1&id=755) and click run button.

##### Run through Command Line:

To install and start using the dolphinnext/rnaseq pipeline by using command line, please follow these steps: [*Installation*](https://github.com/dolphinnext/rnaseq/blob/master/docs/local.md).