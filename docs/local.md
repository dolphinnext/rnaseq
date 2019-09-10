# dolphinnext/rnaseq: Local Configuration
<!-- Install Atom plugin markdown-toc-auto for this ToC -->
<!-- TOC START min:2 max:3 link:true asterisk:true -->
* [Install NextFlow](#install-nextflow)
* [Install the pipeline](#install-the-pipeline)
  * [Automatic](#automatic)
  * [Docker](#docker)
  * [Singularity](#singularity)
* [Run the pipeline](usage.md)
<!-- TOC END -->

## Install NextFlow
Nextflow runs on most POSIX systems (Linux, Mac OSX etc). It can be installed by running the following commands:

```bash
# Make sure that Java v8+ is installed:
java -version

# Install Nextflow
curl -fsSL get.nextflow.io | bash

# Add Nextflow binary to your bin PATH or any accessible path in your environment:
chmod 755 nextflow
mv nextflow ~/bin/
# OR system-wide installation:
# sudo mv nextflow /usr/local/bin
```

See [nextflow.io](https://www.nextflow.io/) for further instructions on how to install and configure Nextflow.

## Install the pipeline

### Automatic
This pipeline itself needs no installation - NextFlow will automatically fetch it from GitHub if `dolphinnext/rnaseq` is specified as the pipeline name.

### Docker
First, install docker on your system: [Docker Installation Instructions](https://docs.docker.com/engine/installation/)

Then, running the pipeline with the option `-profile docker` tells Nextflow to enable Docker for this run. An image containing all of the software requirements will be automatically fetched and used from dockerhub ([https://hub.docker.com/r/dolphinnext/rnaseq](https://hub.docker.com/r/dolphinnext/rnaseq)).

```
nextflow run dolphinnext/rnaseq -profile docker --DOWNDIR /path/to/save/genome-data --reads '*_R{1,2}.fastq.gz' --genome_build mouse_mm10_refseq
```

### Singularity
If you're not able to use Docker then [Singularity](http://singularity.lbl.gov/) is a great alternative.
The process is very similar: running the pipeline with the option `-profile singularity` tells Nextflow to enable singularity for this run. An docker image will be automatically converted into singularity image and used in the pipeline.

```
nextflow run dolphinnext/rnaseq -profile singularity --DOWNDIR /path/to/save/genome-data --reads '*_R{1,2}.fastq.gz' --genome_build mouse_mm10_refseq
```