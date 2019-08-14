RNA-seq pipeline includes Quality Control, rRNA filtering, Genome Alignment using HISAT2, STAR and Tophat2, and estimating gene and isoform expression levels by RSEM and featureCounts.  
  
Steps:
  1. For Quality Control, we use FastQC to create qc outputs. There are optional read quality filtering (trimmomatic), read quality trimming (trimmomatic), adapter removal (cutadapt) processes available.
  2. Bowtie2/Bowtie/STAR is used to count or filter out common RNAs reads (eg. rRNA, miRNA, tRNA, piRNA etc.). 
  3. RSEM is used to align RNA-Seq reads to a reference transcripts and estimates gene and isoform expression levels.
  4. HISAT2, STAR and Tophat2 is used to align RNA-Seq reads to a genome. Optional estimation of gene and isoform expression levels could be done by featureCount.
  5. Genome-wide Bam analysis is done by RseQC, Picard.
  6. Optionally you can create Integrative Genomics Viewer (IGV)  and Genome Browser Files (TDF and Bigwig, respectively)

Program Versions:
  - FastQC v0.11.7
  - Picard-tools v1.131
  - Bedtools v2.25.0
  - Samtools v1.2
  - FASTX Toolkit v0.0.13
  - WigToBigWig v4
  - Pdfbox-App v2.0.0
  - IGV v2.3.32
  - Trimmomatic v0.32
