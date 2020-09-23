# dolphinnext/rnaseq: Usage

## Introduction
Nextflow handles job submissions on SLURM or other environments, and supervises running the jobs. Thus the Nextflow process must run until the pipeline is finished. We recommend that you put the process running in the background through `screen` / `tmux` or similar tool. Alternatively you can run nextflow within a cluster job submitted your job scheduler.

It is recommended to limit the Nextflow Java virtual machines memory. We recommend adding the following line to your environment (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

## Running the pipeline
The typical command for running the pipeline is as follows:

```bash
nextflow run dolphinnext/rnaseq -profile docker --DOWNDIR /path/to/save/genome-data --reads '*_R{1,2}.fastq.gz' --mate 'pair' --genome_build mouse_mm10_refseq
```

If you're running for the first time, you need to enable `--run_checkAndBuild` paramater as follows:

```bash
nextflow run dolphinnext/rnaseq -profile docker --DOWNDIR /path/to/save/genome-data --reads '*_R{1,2}.fastq.gz' --mate 'pair' --genome_build mouse_mm10_refseq --run_checkAndBuild 'yes'
```


This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
results         # Finished results 
.nextflow_log   # Log file from Nextflow
```

### Updating the pipeline
When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version. In order to download latest version of the pipeline you need to run following command:

```bash
nextflow pull dolpinnext/rnaseq
```

## Main arguments

### `-profile`
Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments. Note that multiple profiles can be loaded, for example: `-profile docker, test` - the order of arguments is important!

If `-profile` is not specified at all the pipeline will be run locally and expects all software to be installed and available on the `PATH`.

* `docker`
  * A generic configuration profile to be used with [Docker](http://docker.com/)
  * Pulls software from Dockerhub: [`dolphinnext/rnaseq`](http://hub.docker.com/r/dolphinnext/rnaseq/)
* `singularity`
  * A generic configuration profile to be used with [Singularity](http://singularity.lbl.gov/)
  * Pulls software from DockerHub
* `test`
  * A profile with a complete configuration for automated testing
  ```bash
  ## First, download sample fastq files into `inputs` folder with the following command:
  mkdir -p inputs && cd inputs && wget https://galaxyweb.umassmed.edu/pub/dnext_data/test/reads/control_rep1.1.fastq.gz https://galaxyweb.umassmed.edu/pub/dnext_data/test/reads/exper_rep1.1.fastq.gz  && cd ..
  ## Start testing pipeline:
  nextflow run dolphinnext/rnaseq -profile docker,test 
  ## In the test profile, --reads parameter assinged as: 'inputs/*.fastq.gz'
  ```

### `--reads`
Use this to specify the location of your input FastQ files. For example:

```bash
--reads 'path/to/data/sample_*_{1,2}.fastq' --mate 'pair'
--reads 'path/to/data/sample_*.fastq' --mate 'single'
```

Please note the following requirements:

1. The path must be enclosed in quotes
2. The path must have at least one `*` wildcard character
3. When using the pipeline with paired end data, the path must use `{1,2}` notation to specify read pairs.


### `--mate`
Two options (single or pair) available for `--mate` parameter. If you have single-end data, you need to specify as 'single' and for paired-end data, you need to specify as 'pair'. For example:

```bash
--reads 'path/to/data/sample_*_{1,2}.fastq' --mate 'pair'
--reads 'path/to/data/sample_*.fastq' --mate 'single'
```

It is not possible to run a mixture of single-end and paired-end files in one run.


## Reference genomes

### `--genome_build` 
There are 8 different species supported in the UMMS-Biocore references. To run the pipeline, you must specify which to use with the `--genome_build` flag.

List of genomes that are supported are:

* Human
  * `--genome_build human_hg19_refseq`
  * `--genome_build human_hg38_gencode_v28`
  * `--human_hg38_gencode_v34`
* Mouse
  * `--genome_build mouse_mm10_refseq`
  * `--genome_build mouse_mm10_gencode_m25`
* Rat
  * `--genome_build rat_rn6_refseq`
  * `--genome_build rat_rn6_ensembl_v86`
* Zebrafish
  * `--genome_build zebrafish_GRCz11_ensembl_v95`
  * `--genome_build zebrafish_GRCz11_refseq`
  * `--genome_build zebrafish_GRCz11_v4.3.2`
* C. elegans
  * `--genome_build c_elegans_ce11_ensembl_ws245`
* S. cerevisiae
  * `--genome_build s_cerevisiae_sacCer3_refseq` 
* S. pombe
  * `--genome_build s_pombe_ASM294v2_ensembl_v31` 
* D. melanogaster
  * `--genome_build d_melanogaster_dm6_refseq`

Note: For new genome requests, please send e-mail to UMMS-Biocore(biocore@umassmed.edu).

### `--DOWNDIR` `--run_checkAndBuild`
If your indexes are not build before, you can enable `--run_checkAndBuild` by assinging it's value to 'yes' which will check genome files in `--DOWNDIR` and download into that directory. Afterwards it will start building indexes based on the selected parameters in the pipeline. 


### `--star_index`, `--bowtie_index`, `--bowtie2_index`, `--hisat2_index`, `--rsem_ref_using_bowtie_index`, `--rsem_ref_using_bowtie2_index`, `--rsem_ref_using_star_index`, `--genome`, `--gtf`, `--bed`, `--genome_sizes`, `--commondb`
If you prefer, you can specify the full path to your reference genome and disable `--run_checkAndBuild` option.

```bash
--genome '[path to Fasta reference]' \
--genome_sizes '[path to genome_sizes file]' \
--gtf '[path to GTF file]' \
--bed '[path to bed12 file]' \
--commondb '[path to commondb directory when Bowtie/Bowtie2 indexes found for common RNA's (eg. rRNA, miRNA, tRNA, etc.)] \

--star_index '[path to STAR index]' \
--bowtie_index '[path to Bowtie index]' \
--bowtie2_index '[path to Bowtie index]' \
--hisat2_index '[path to HISAT2 index]' \
--kallisto_index '[path to Kallisto index]' \
--rsem_ref_using_bowtie_index '[path to RSEM reference build with bowtie index]' \
--rsem_ref_using_bowtie2_index '[path to RSEM reference build with bowtie2 index]' \
--rsem_ref_using_star_index    '[path to RSEM reference build with STAR index]' \

```

## Alignment tool
By default, the pipeline uses [STAR](https://github.com/alexdobin/STAR) to align the raw FastQ reads to the reference genome. STAR is fast and common, but requires a lot of memory to run, typically around 38GB for the Human hg19 reference genome.

If you prefer, you can use [HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml) or [TOPHAT2](https://ccb.jhu.edu/software/tophat/) as the alignment tool. Both tools developed by the same group and have a much smaller memory footprint compared to STAR. However, HISAT2 is recommended by the same group.  

You can choose multiple aligner to compare their results by enabling/disabling following parameters:
```bash
To enable HISAT2  : `--run_HISAT2 yes`
To disable HISAT2 : `--run_HISAT2 no`
To enable STAR    : `--run_STAR yes`
To disable STAR   : `--run_STAR no`
To enable Tophat2 : `--run_Tophat yes`
To disable Tophat2: `--run_Tophat no`
```

## Transcripts Quantification
By default, RSEM is used to align RNA-Seq reads to a reference transcripts and estimates gene and isoform expression levels. 
Besides, you can enable featureCounts after running HISAT2 and/or STAR and/or Tophat2 as well.
Alternatively, Kallisto could be used for quantifying abundances of transcripts based on pseudoalignments, without the need for alignment.

You can choose multiple tools to compare their results by enabling/disabling following parameters:
```bash
To enable RSEM  : `--run_RSEM yes`
To disable RSEM : `--run_RSEM no`
To enable featureCounts after running HISAT   : `--run_HISAT2 yes --run_FeatureCounts_after_Hisat2 yes`
To disable featureCounts after running HISAT  : `--run_FeatureCounts_after_Hisat2 no`
To enable featureCounts after running STAR    : `--run_STAR yes --run_FeatureCounts_after_STAR yes`
To disable featureCounts after running STAR   : `--run_FeatureCounts_after_STAR no`
To enable featureCounts after running Tophat2 : `--run_Tophat yes --run_FeatureCounts_after_Tophat2 yes`
To disable featureCounts after running Tophat2: `--run_FeatureCounts_after_Tophat2 no`
To enable Kalliso  : `--run_Kalliso yes`
To disable Kalliso : `--run_Kalliso no`
```

## Feature Counts Parameters
You can change feature count parameters by assigning new parameters to following options:

```bash
# FeatureCounts after Hisat2 Parameters:
--BAM_Analysis_Hisat2_featureCounts_Prep.run_name [array @default:["gene_id","transcript_id"] ]
# Prefix for run output

--BAM_Analysis_Hisat2_featureCounts_Prep.run_parameters =  [array @default:["-g gene_id -s 0 -Q 20 -T 2 -B -d 50 -D 1000 -C --fracOverlap 0 --minOverlap 1","-g transcript_id -s 0 -Q 20 -T 2 -B -d 50 -D 1000 -C --fracOverlap 0 --minOverlap 1"]]  
# -s Indicate strand-specific read counting: 0 (unstranded, default), 1 (stranded) and 2 (reversely stranded) 
# -Q The minimum mapping quality score 
# -T Number of threads 
# -B requireBothEndsMapped 
# -C countChimericFragments 
# −−fracOverlap Minimum fraction of overlapping bases 
# −−minOverlap Minimum number of overlapping bases

# FeatureCounts after Tophat2 Parameters:
--BAM_Analysis_Tophat2_featureCounts_Prep.run_name [array @default:["gene_id","transcript_id"] ]
--BAM_Analysis_Tophat2_featureCounts_Prep.run_parameters =  [array @default:["-g gene_id -s 0 -Q 20 -T 2 -B -d 50 -D 1000 -C --fracOverlap 0 --minOverlap 1","-g transcript_id -s 0 -Q 20 -T 2 -B -d 50 -D 1000 -C --fracOverlap 0 --minOverlap 1"]]  

# FeatureCounts after STAR Parameters:
--BAM_Analysis_STAR_featureCounts_Prep.run_name [array @default:["gene_id","transcript_id"] ]
--BAM_Analysis_STAR_featureCounts_Prep.run_parameters =  [array @default:["-g gene_id -s 0 -Q 20 -T 2 -B -d 50 -D 1000 -C --fracOverlap 0 --minOverlap 1","-g transcript_id -s 0 -Q 20 -T 2 -B -d 50 -D 1000 -C --fracOverlap 0 --minOverlap 1"]]  
```

## Kalliso Parameters
You can change Kalliso parameters by assigning new parameters to following options:

```bash
params.Kallisto_module_kallisto_quant.single_or_paired_end_reads = [@options:"single","pair" @default:"pair"] 
# When single is selected, the fragment_length and standard_deviation parameters are REQUIRED.

params.Kallisto_module_kallisto_quant.fragment_length = [int @default:200] 
# Estimated average fragment length. Required when using single reads. Typical value: 200.

params.Kallisto_module_kallisto_quant.standard_deviation = [int @default:30] 
# Estimated standard deviation of fragment length. Required when using single reads. Typical value: 30.

params.Kallisto_module_kallisto_quant.kallisto_parameters = [string @default:"--threads 4"] 
# Custom Kallisto parameters.

params.Kallisto_module_kallisto_quant.genomebam = [@options:"true","false" @default:"true"] 
# If true is selected, Kallisto will project pseudoalignments to genome sorted BAM file."
```

## Adapter Removal
If specific Adapter Removal is required, you can enable trimmomatic and enter the adapter sequence. 

```bash
To enable adapter_removal: 
--run_Adapter_Removal "yes"

--Adapter_Trimmer_Quality_Module_Adapter_Removal.phred = [@options:33,64 @default:33]
# Specifies the fastq quality encoding. Default is 33 which is now almost universally used, and 64 which is used in some older Illumina data

--Adapter_Trimmer_Quality_Module_Adapter_Removal.Adapter_Sequence [string]
# You can enter a single sequence or multiple sequences in different lines. Reverse sequences will not be removed.

--Adapter_Trimmer_Quality_Module_Adapter_Removal.min_length [int @default:10]
# Specifies the minimum length of reads to be kept

--Adapter_Trimmer_Quality_Module_Adapter_Removal.seed_mismatches [int @default:1]
# Specifies the maximum mismatch count which will still allow a full match to be performed

--Adapter_Trimmer_Quality_Module_Adapter_Removal.palindrome_clip_threshold [int @default:30]
# Specifies how accurate the match between the two -adapter ligated- reads must be for PE palindrome read alignment

--Adapter_Trimmer_Quality_Module_Adapter_Removal.simple_clip_threshold [int @default:5]
# Specifies how accurate the match between any adapter etc. sequence must be against a read.

--Adapter_Trimmer_Quality_Module_Adapter_Removal.discard_non_clipped [@options:"yes","no" @default:"yes"]
# Discard_non_clipped sequences (keep only sequences which contained the adapter)
```

## Trimmer
Optianally, you can trim your reads by defining trimming lenghts as shown at below: 

```bash

--run_Trimmer [@options:"yes","no" @default:"no"]
# Enables Trimmer by setting this parameter as "yes"

--Adapter_Trimmer_Quality_Module_Trimmer.phred = [@options:33,64 @default:33]
# Specifies the fastq quality encoding. Default is 33 which is now almost universally used, and 64 which is used in some older Illumina data

For Single End Reads  : 
--Adapter_Trimmer_Quality_Module_Trimmer.single_or_paired_end_reads "single"
--Adapter_Trimmer_Quality_Module_Trimmer.trim_length_5prime [int]
--Adapter_Trimmer_Quality_Module_Trimmer.trim_length_3prime [int]

For Paired End Reads  : 
--Adapter_Trimmer_Quality_Module_Trimmer.single_or_paired_end_reads "pair"
--Adapter_Trimmer_Quality_Module_Trimmer.trim_length_5prime_R1 [int]
--Adapter_Trimmer_Quality_Module_Trimmer.trim_length_3prime_R1 [int]
--Adapter_Trimmer_Quality_Module_Trimmer.trim_length_5prime_R2 [int]
--Adapter_Trimmer_Quality_Module_Trimmer.trim_length_3prime_R2 [int]
```

## Quality Filtering
Optianally, you can trim your reads based on their quality. Trimmomatic works on both paired-end and single ended data. Alternatively fastx option (fastx_toolkit) could be used for single reads. 

```bash
To use Trimmomatic  : 
--run_Quality_Filtering "yes"
--Adapter_Trimmer_Quality_Module_Quality_Filtering.tool "trimmomatic"

--Adapter_Trimmer_Quality_Module_Quality_Filtering.phred = [@options:33,64 @default:33]
# Specifies the fastq quality encoding. Default is 33 which is now almost universally used, and 64 which is used in some older Illumina data

--Adapter_Trimmer_Quality_Module_Quality_Filtering.window_size [int @default:10]
# Performs a sliding window trimming approach. It starts scanning at the 5' end and clips the read once the average quality within the window falls below a threshold (=required_quality).

--Adapter_Trimmer_Quality_Module_Quality_Filtering.required_quality_for_window_trimming [int @default:15]
# Specifies the average quality required for window trimming approach

--Adapter_Trimmer_Quality_Module_Quality_Filtering.leading [int @default:5]
# Cut bases off the start of a read, if below a threshold quality

--Adapter_Trimmer_Quality_Module_Quality_Filtering.trailing [int @default:5]
# Cut bases off the end of a read, if below a threshold quality

--Adapter_Trimmer_Quality_Module_Quality_Filtering.minlen [int @default:36]
# Specifies the minimum length of reads to be kept
```

```bash
To use fastx_toolkit  : 
--run_Quality_Filtering "yes"
--Adapter_Trimmer_Quality_Module_Quality_Filtering.tool "fastx"
--Adapter_Trimmer_Quality_Module_Quality_Filtering.minQuality [int @default:20]
# Minimum quality score to keep reads

--Adapter_Trimmer_Quality_Module_Quality_Filtering.minPercent [int @default:100]
# Minimum percent of bases that must have entered minQuality
```

## Sequential Mapping
Optianally,Bowtie2/Bowtie/STAR is used to count or filter out common RNAs reads (eg. rRNA, miRNA, tRNA, piRNA etc.). You need to specify mapping set by entering following paramters in array format.

```bash
--run_Sequential_Mapping "yes"
--Sequential_Mapping_Module_Sequential_Mapping.remove_duplicates [@options:"yes","no" @default:"no"] 
# Duplicates (both PCR and optical) will be removed from alignment file (bam) and separate count table will be created for comparison

--Sequential_Mapping_Module_Sequential_Mapping._select_sequence [array @options:"rRNA","ercc","miRNA","tRNA","piRNA","snRNA","rmsk", "custom"]   
# Sequence Set for Mapping. eg. ["rRNA", "rmsk", "custom"]

--Sequential_Mapping_Module_Sequential_Mapping.index_directory [array]
# If custom sequence is defined please enter index directory of custom sequence(full path), otherwise you need to enter empty string. The index directory must include the full path and the name of the index file must only be the prefix of the fasta or index file. Index files and fasta files also need to have the same prefix.For STAR alignment, gtf file which has the same prefix, must be found in same directory. eg. ["", "", "/share/custom_seq_dir"]

--Sequential_Mapping_Module_Sequential_Mapping.name_of_the_index_file [array]  
# If custom sequence is defined please enter name of the index or fasta file (prefix), otherwise you need to enter selected sequence as string. eg. ["rRNA", "rmsk", "custom_seq_prefix"]

--Sequential_Mapping_Module_Sequential_Mapping._aligner =  [array @options:"bowtie","bowtie2" @default:"bowtie2"] 
# Aligner set for mapping: eg. ["bowtie", "bowtie2", "bowtie2"]

--Sequential_Mapping_Module_Sequential_Mapping.aligner_Parameters [array]
# Aligner parameters." eg. ["--threads 1","-N 1","-N 1"]

--params.Sequential_Mapping_Module_Sequential_Mapping.description [array] 
# Description of index file (please don't use comma or quotes in this field). eg. ["rRNA", "rmsk", "custom_seq_explanation"]

--Sequential_Mapping_Module_Sequential_Mapping.filter_Out =  "[array @options:"Yes","No" @default:"Yes"] 
# Select whether or not you want the reads mapped to this index filtered out of your total reads.

```

## TDF Conversion for IGV Genome Browser
Optionally, you can convert bam files to TDF for IGV Genome Browser visualization by using IGVtools.
```bash
--run_IGV_TDF_Conversion "yes"
## For RSEM BAM output
--BAM_Analysis_RSEM_IGV_BAM2TDF_converter.igv_extention_factor [int @default:0]
# The read or feature is extended by the specified distance in bp prior to counting. This option is useful for chip-seq and rna-seq applications. The value is generally set to the average fragment length of the library minus the average read length.
--BAM_Analysis_RSEM_IGV_BAM2TDF_converter.igv_window_size [int @default:5]
# The window size over which coverage is averaged.

## For Tophat2 BAM output
--BAM_Analysis_Tophat2_IGV_BAM2TDF_converter.igv_extention_factor [int @default:0]
--BAM_Analysis_Tophat2_IGV_BAM2TDF_converter.igv_window_size [int @default:5]

## For STAR BAM output
--BAM_Analysis_STAR_IGV_BAM2TDF_converter.igv_extention_factor [int @default:0]
--BAM_Analysis_STAR_IGV_BAM2TDF_converter.igv_window_size [int @default:5]

## For HISAT2 BAM output
--BAM_Analysis_Hisat2_IGV_BAM2TDF_converter.igv_extention_factor [int @default:0]
--BAM_Analysis_Hisat2_IGV_BAM2TDF_converter.igv_window_size [int @default:5]

```

## BigWig Conversion for UCSC Genome Browser
Optionally, you can convert bam files to bigwig files for UCSC Genome Browser visualization.

```bash
--run_BigWig_Conversion "yes"
```

## BigWig Conversion for UCSC Genome Browser
Optionally, you can convert bam files to bigwig files for UCSC Genome Browser visualization.

```bash
--run_BigWig_Conversion "yes"
```

## RSeQC Analysis
Optionally, you can enable RSeQC to calculate how mapped reads were distributed over genome feature (like CDS exon, 5’UTR exon, 3’ UTR exon, Intron, Intergenic regions).

```bash
--run_RSeQC "yes"
```

## Picard Analysis
Optionally, you can enable Picard to calculate multiple metrics such as CollectAlignmentSummaryMetrics, CollectInsertSizeMetrics, QualityScoreDistribution, MeanQualityByCycle, CollectBaseDistributionByCycle, CollectGcBiasMetrics, RnaSeqMetrics, CollectSequencingArtifactMetrics, and CollectQualityYieldMetrics. 

```bash
--run_Picard_CollectMultipleMetrics "yes"
```


## Other command line parameters

### `--outdir`
The output directory where the results will be saved.

### `-resume`
Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`
Specify the path to a specific config file (this is a core NextFlow command). 

Note - you can use this to override pipeline defaults.

## Stand-alone scripts
The `dolphinnext/tools` repository contains some scripts used by the pipeline which may also be run manually:

* `gtf2bed`
  * Script used to generate the BED12 reference files used by RSeQC. Takes a `.gtf` file as input