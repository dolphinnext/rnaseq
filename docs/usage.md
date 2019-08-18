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
To enable HISAT2  : `--run_HISAT2 yes`
To disable HISAT2 : `--run_HISAT2 no`
To enable STAR    : `--run_STAR yes`
To disable STAR   : `--run_STAR no`
To enable Tophat2 : `--run_Tophat yes`
To disable Tophat2: `--run_Tophat no`

## Transcripts Quantification
By default, RSEM is used to align RNA-Seq reads to a reference transcripts and estimates gene and isoform expression levels. But you can enable featureCounts after running HISAT2 and/or STAR and/or Tophat2 as well.

You can choose multiple tools to compare their results by enabling/disabling following parameters:
To enable RSEM  : `--run_RSEM yes`
To disable RSEM : `--run_RSEM no`
To enable featureCounts after running HISAT   : `--run_HISAT2 yes --run_FeatureCounts_after_Hisat2 yes`
To disable featureCounts after running HISAT  : `--run_FeatureCounts_after_Hisat2 no`
To enable featureCounts after running STAR    : `--run_STAR yes --run_FeatureCounts_after_STAR yes`
To disable featureCounts after running STAR   : `--run_FeatureCounts_after_STAR no`
To enable featureCounts after running Tophat2 : `--run_Tophat yes --run_FeatureCounts_after_Tophat2 yes`
To disable featureCounts after running Tophat2: `--run_FeatureCounts_after_Tophat2 no`


## Adapter Trimming
If specific additional trimming is required (for example, from additional tags),
you can use any of the following command line parameters. These affect the command
used to launch TrimGalore!

### `--clip_r1 [int]`
Instructs Trim Galore to remove bp from the 5' end of read 1 (or single-end reads).

### `--clip_r2 [int]`
Instructs Trim Galore to remove bp from the 5' end of read 2 (paired-end reads only).

### `--three_prime_clip_r1 [int]`
Instructs Trim Galore to remove bp from the 3' end of read 1 _AFTER_ adapter/quality trimming has been performed.

### `--three_prime_clip_r2 [int]`
Instructs Trim Galore to re move bp from the 3' end of read 2 _AFTER_ adapter/quality trimming has been performed.



## Skipping QC steps
The pipeline contains a large number of quality control steps. Sometimes, it may not be desirable to run all of them if time and compute resources are limited.
The following options make this easy:

* `--skip_qc` -                Skip **all QC steps**, apart from MultiQC
* `--skip_fastqc` -            Skip FastQC
* `--skip_rseqc` -             Skip RSeQC
* `--skip_genebody_coverage` - Skip calculating the genebody coverage
* `--skip_preseq` -            Skip Preseq
* `--skip_dupradar` -          Skip dupRadar (and Picard MarkDups)
* `--skip_edger` -             Skip edgeR MDS plot and heatmap
* `--skip_multiqc` -           Skip MultiQC

## Job resources
### Automatic resubmission
Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

### Custom resource requests
Wherever process-specific requirements are set in the pipeline, the default value can be changed by creating a custom config file. See the files hosted at [`nf-core/configs`](https://github.com/nf-core/configs/tree/master/conf) for examples.

If you are likely to be running `nf-core` pipelines regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter (see definition below). You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

If you have any questions or issues please send us a message on [`Slack`](https://nf-core-invite.herokuapp.com/).


## Other command line parameters

### `--outdir`
The output directory where the results will be saved.


### `-resume`
Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

**NB:** Single hyphen (core Nextflow option)

### `-c`
Specify the path to a specific config file (this is a core NextFlow command).

**NB:** Single hyphen (core Nextflow option)

Note - you can use this to override pipeline defaults.

### `--custom_config_version`
Provide git commit id for custom Institutional configs hosted at `nf-core/configs`. This was implemented for reproducibility purposes. Default is set to `master`.

```bash
## Download and use config file with following git commid id
--custom_config_version d52db660777c4bf36546ddb188ec530c3ada1b96
```

### `--custom_config_base`
If you're running offline, nextflow will not be able to fetch the institutional config files
from the internet. If you don't need them, then this is not a problem. If you do need them,
you should download the files from the repo and tell nextflow where to find them with the
`custom_config_base` option. For example:

```bash
## Download and unzip the config files
cd /path/to/my/configs
wget https://github.com/nf-core/configs/archive/master.zip
unzip master.zip

## Run the pipeline
cd /path/to/my/data
nextflow run /path/to/pipeline/ --custom_config_base /path/to/my/configs/configs-master/
```

> Note that the nf-core/tools helper package has a `download` command to download all required pipeline
> files + singularity containers + institutional configs in one go for you, to make this process easier.

### `--max_memory`
Use to set a top-limit for the default memory requirement for each process.
Should be a string in the format integer-unit. eg. `--max_memory '8.GB'`

### `--max_time`
Use to set a top-limit for the default time requirement for each process.
Should be a string in the format integer-unit. eg. `--max_time '2.h'`

### `--max_cpus`
Use to set a top-limit for the default CPU requirement for each process.
Should be a string in the format integer-unit. eg. `--max_cpus 1`

### `--hisatBuildMemory`
Required amount of memory in GB to build HISAT2 index with splice sites.
The HiSAT2 index build can proceed with or without exon / splice junction information.
To work with this, a very large amount of memory is required.
If this memory is not available, the index build will proceed without splicing information.
The `--hisatBuildMemory` option changes this threshold. By default it is `200GB` - if your system
`--max_memory` is set to `128GB` but your genome is small enough to build using this, then you can
allow the exon build to proceed by supplying `--hisatBuildMemory 100GB`

### `--subsampFilesizeThreshold`
This parameter defines the threshold in BAM file size (in bytes) at which data subsampling is used prior to the RSeQC `gene_body_coverage` step. This step is done to speed up and reduce compute resources for the gene body coverage analysis .
For very large files this means, that the BAM file will be subsampled to compute the `gene_body_coverage`, for small files there will not be a subsampling step.
By default this parameter is set to `10000000000` - ten gigabytes.

### `--sampleLevel`
Used to turn of the edgeR MDS and heatmap. Set automatically when running on fewer than 3 samples.

### `--plaintext_email`
Set to receive plain-text e-mails instead of HTML formatted.

### `--monochrome_logs`
Set to disable colourful command line output and live life in monochrome.

### `--multiqc_config`
Specify a path to a custom MultiQC configuration file.

## Stand-alone scripts
The `bin` directory contains some scripts used by the pipeline which may also be run manually:

* `gtf2bed`
  * Script used to generate the BED12 reference files used by RSeQC. Takes a `.gtf` file as input
* `dupRadar.r`
  * dupRadar script used in the _dupRadar_ pipeline process.
* `edgeR_heatmap_MDS.r`
  * edgeR script used in the _Sample Correlation_ process
