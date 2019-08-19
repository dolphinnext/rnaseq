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
  * Pulls software from Dockerhub: [`nfcore/rnaseq`](http://hub.docker.com/r/nfcore/rnaseq/)
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
There are 5 different species supported in the UMMS-Biocore references. To run the pipeline, you must specify which to use with the `--genome_build` flag.

List of genomes that are supported are:

* Human
  * `--genome_build human_hg19_refseq`
  * `--genome_build human_hg38_genecode_v28`
* Mouse
  * `--genome_build mouse_mm10_refseq`
* Rat
  * `--genome_build rat_rn6_refseq`
  * `--genome_build rat_rn6_ensembl_v86`
* Zebrafish
  * `--genome_build 'zebrafish_danRer10'`


### `--DOWNDIR` `--run_checkAndBuild`
If your indexes are not build before, you can enable `--run_checkAndBuild` by assinging it's value to 'yes' which will check genome files in `--DOWNDIR` and download into that directory. Afterwards it will start building indexes based on the selected parameters in the pipeline. 


### `--star_index`, `--bowtie_index`, `--bowtie2_index`, `--hisat2_index`, `--rsem_ref_using_bowtie_index`, `--rsem_ref_using_bowtie2_index`, `--rsem_ref_using_star_index`, `--genome`, `--gtf`, `--bed`, `--genome_sizes`, `--commondb`
If you prefer, you can specify the full path to your reference genome and without disable `--run_checkAndBuild` option.

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
--rsem_ref_using_bowtie_index '[path to RSEM reference build with bowtie index]' \
--rsem_ref_using_bowtie2_index '[path to RSEM reference build with bowtie2 index]' \
--rsem_ref_using_star_index    '[path to RSEM reference build with STAR index]' \

```



## Alignment tool
By default, the pipeline uses [STAR](https://github.com/alexdobin/STAR) to align the raw FastQ reads to the reference genome. STAR is fast and common, but requires a lot of memory to run, typically around 38GB for the Human hg19 reference genome.

If you prefer, you can use [HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml) or TOPHAT2 as the alignment tool instead. Both tools developed by the same group, and HISAT2 has a much smaller memory footprint.

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
By default, RSEM is used to align RNA-Seq reads to a reference transcripts and estimates gene and isoform expression levels. But you can enable featureCounts after running HISAT2 and/or STAR and/or Tophat2 as well.

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
```

## Adapter Removal
If specific Adapter Removal is required, you can enable trimmomatic and enter the adapter sequence. 

```bash
To enable adapter_removal  : `--run_Adapter_Removal yes`
```
#### `--Adapter_Trimmer_Quality_Module_Adapter_Removal.Adapter_Sequence [string]`
You can enter a single sequence or multiple sequences in different lines. Reverse sequences will not be removed.

#### `--Adapter_Trimmer_Quality_Module_Adapter_Removal.min_length [int]`
Specifies the minimum length of reads to be kept

#### `--Adapter_Trimmer_Quality_Module_Adapter_Removal.seed_mismatches [int]`
Specifies the maximum mismatch count which will still allow a full match to be performed

#### `--Adapter_Trimmer_Quality_Module_Adapter_Removal.palindrome_clip_threshold [int (default=30)]`
Specifies how accurate the match between the two -adapter ligated- reads must be for PE palindrome read alignment

#### `--Adapter_Trimmer_Quality_Module_Adapter_Removal.simple_clip_threshold [int (default=5)]`
Specifies how accurate the match between any adapter etc. sequence must be against a read.

#### `--Adapter_Trimmer_Quality_Module_Adapter_Removal.discard_non_clipped [@dropdown @options:"yes","no" @default="yes"]`
Discard_non_clipped sequences (keep only sequences which contained the adapter)




// Process Parameters for Adapter_Trimmer_Quality_Module_Trimmer:
params.Adapter_Trimmer_Quality_Module_Trimmer.single_or_paired_end_reads =  ""  //* @dropdown @options:"single","pair" 
params.Adapter_Trimmer_Quality_Module_Trimmer.trim_length_5prime =  0  //* @input @description:"Trimming length from 5' end"  
params.Adapter_Trimmer_Quality_Module_Trimmer.trim_length_3prime =  0  //* @input @description:"Trimming length from 3' end"  
params.Adapter_Trimmer_Quality_Module_Trimmer.trim_length_5prime_R1 =  0  //* @input @description:"Trimming length from 5' end of R1 reads"  
params.Adapter_Trimmer_Quality_Module_Trimmer.trim_length_3prime_R1 =  0  //* @input @description:"Trimming length from 3' end of R1 reads"  
params.Adapter_Trimmer_Quality_Module_Trimmer.trim_length_5prime_R2 =  0  //* @input @description:"Trimming length from 5' end of R2 reads"  
params.Adapter_Trimmer_Quality_Module_Trimmer.trim_length_3prime_R2 =  0  //* @input @description:"Trimming length from 3' end of R2 reads" 

// Process Parameters for Adapter_Trimmer_Quality_Module_Quality_Filtering:
params.Adapter_Trimmer_Quality_Module_Quality_Filtering.tool =  "trimmomatic"  //* @dropdown @options:"trimmomatic","fastx" @description:"Choose quality removal tool to be used. Note:fastx option (fastx_toolkit fastq_quality_filter) is not suitable for paired reads." 
params.Adapter_Trimmer_Quality_Module_Quality_Filtering.window_size =  10  //* @input @description:"Performs a sliding window trimming approach. It starts scanning at the 5' end and clips the read once the average quality within the window falls below a threshold (=required_quality)."  
params.Adapter_Trimmer_Quality_Module_Quality_Filtering.required_quality_for_window_trimming =  15  //* @input @description:"specifies the average quality required for window trimming approach" 
params.Adapter_Trimmer_Quality_Module_Quality_Filtering.leading =  5  //* @input @description:"Cut bases off the start of a read, if below a threshold quality" 
params.Adapter_Trimmer_Quality_Module_Quality_Filtering.trailing =  5  //* @input @description:"Cut bases off the end of a read, if below a threshold quality"  
params.Adapter_Trimmer_Quality_Module_Quality_Filtering.minlen =  36  //* @input @description:"Specifies the minimum length of reads to be kept"  
params.Adapter_Trimmer_Quality_Module_Quality_Filtering.minQuality =  20  //* @input @description:"Minimum quality score to keep reads"
params.Adapter_Trimmer_Quality_Module_Quality_Filtering.minPercent =  100  //* @input @description:"Minimum percent of bases that must have entered minQuality"

// Process Parameters for Sequential_Mapping_Module_Sequential_Mapping:
params.Sequential_Mapping_Module_Sequential_Mapping.remove_duplicates =  "no"  //* @dropdown @description:"Duplicates (both PCR and optical) will be removed from alignment file (bam) and separate count table will be created for comparison" @title:"General Mapping Options" @options:{"yes","no"}
params.Sequential_Mapping_Module_Sequential_Mapping.remove_duplicates_based_on_UMI_after_mapping =  "no"  //* @dropdown @description:"UMI extract process should have executed before this step. Read headers should have UMI tags which are separated with underscore.(eg. NS5HGY:2:11_GTATAACCTT)" @options:{"yes","no"}
params.Sequential_Mapping_Module_Sequential_Mapping._select_sequence =  ""  //* @dropdown @description:"Select sequence for mapping" @title:"Sequence Set for Mapping" @options:{"rRNA","ercc","miRNA","tRNA","piRNA","snRNA","rmsk","genome",custom"},{_nucleicAcidType="dna","ercc","rmsk","genome","custom"}
params.Sequential_Mapping_Module_Sequential_Mapping.index_directory =  ""  //* @input  @description:"index directory of sequence(full path)" @tooltip:"The index directory must include the full path and the name of the index file must only be the prefix of the fasta or index file. Index files and Fasta files also need to have the same prefix.For STAR alignment, gtf file which has the same prefix, must be found in same directory" 
params.Sequential_Mapping_Module_Sequential_Mapping.name_of_the_index_file =  ""  //* @input  @autofill:{_select_sequence=("rRNA","ercc","miRNA","tRNA","piRNA","snRNA","rmsk","genome"), _select_sequence},{_select_sequence="custom", " "} @description:"Name of the index or fasta file (prefix)" @tooltip:"The index directory must include the full path and the name of the index file must only be the prefix of the fasta or index file. Index files and Fasta files also need to have the same prefix.For STAR alignment, gtf file which has the same prefix, must be found in same directory" 
params.Sequential_Mapping_Module_Sequential_Mapping._aligner =  "bowtie2"  //* @dropdown @description:"Select aligner tool"  @options:{_select_sequence=("rRNA","ercc","miRNA","tRNA","piRNA","snRNA","rmsk"),"bowtie","bowtie2"},{_select_sequence=("genome","custom"),"bowtie","bowtie2","STAR"}
params.Sequential_Mapping_Module_Sequential_Mapping.aligner_Parameters =  ""  //* @input @description:"Aligner parameters." @autofill:{_aligner="bowtie", "--threads 1"},{_aligner="bowtie2", "-N 1"},{_aligner="STAR", "--runThreadN 1"} 
params.Sequential_Mapping_Module_Sequential_Mapping.description =  ""  //* @input @autofill:{_select_sequence=("rRNA","ercc","miRNA","tRNA","piRNA","snRNA","rmsk","genome"), _select_sequence},{_select_sequence="custom", " "} @description:"Description of index file (please don't use comma or quotes in this field" 
params.Sequential_Mapping_Module_Sequential_Mapping.filter_Out =  "Yes"  //* @dropdown @dropdown @options:"Yes","No" @description:"Select whether or not you want the reads mapped to this index filtered out of your total reads." 




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

