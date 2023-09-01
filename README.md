# OHRBID
One Health Research into Bacterial Infectious Diseases Lab, University of Glasgow

## Bactocap Nextflow Workflow

The bactocap nextflow workflow is for the analysis and SNP calling of Illumina short-read sequences. It is a translation of Matej Medvecky's bactocap python pipeline and the following scripts (https://github.com/matejmedvecky/anthraxdiversityscripts).

## Quick Start

A version of nextflow that supports DSL2 is required. It is recommended to run the pipeline with `NXF_VER=22.10.4`, as the pipeline has been tested using this version.
E.g. to download
```
export NXF_VER="22.10.4"
curl -fsSL https://get.nextflow.io | bash
```
To download this workflow
```
git clone https://github.com/tlforde/OHRBID.git
```

The workflow is designed to run with either docker `-profile docker`, singularity `-profile singularity` or conda `-profile conda`. An additional profile has also been added for the CLIMB-BIG-DATA notebook servers `-profile climbnotebook` (this has been configured to be compatible with the CLIMB-BIG-DATA notebook global nextflow config; the workflow will run on the Kubernetes infrastructure).
The parameters for the workflow are set in the `nextflow.config`.
The container images are pulled from the CLIMB-BIG-DATA quay.io registry and a singularity cache directory is set in the `nextflow.config`.

E.g. to run the workflow with singularity:
```
nextflow run main.nf -profile singularity
```

## Parameters

The following params must be set by the user in the `nextflow.config`

* **inputDir**<br /> 
Directory containing input fastq 

* **pattern**<br />
Glob pattern to match fastq files in inputDir

* **outputDir**<br /> 
Directory for the output results

* **refDir**<br /> 
The path to supplementary data for the workflow. All databases and supplementary files that the workflow is dependent on should be present in this directory,
e.g. bowtie index, reference fasta, genbank file, custom txt file for discarded regions

* **noMonomorphic**<br /> 
Exclude monomorphic sites. Values: "yes" or "no"

* **discRegFileCustom**<br /> 
File name for custom txt file for discarded regions, set to null if no file

* **bowtieIndexName**<br /> 
Name of bowtie index, EXCLUDING .bt2 file extension 

* **refGenomeName**<br /> 
Name of reference genome

* **genBank**<br /> 
Name of genbank file
  
There are several other params whose value can optionally be changed, described and defined in the `nextflow.config`

## Additional Information

### Containers ###
This workflow uses the CLIMB-BIG-DATA docker container images which are hosted on quay.io (https://quay.io/organization/climb-big-data).
The Dockerfiles used to build the docker images are hosted here: https://github.com/MRC-CLIMB/containers

The tag of the container images refers to the version of the software. E.g.
https://github.com/tlforde/OHRBID/blob/83437f1ef2da3d7d3e7f531fd72df500ae2b4b6f/nextflow.config#L85-L87
uses version 0.39 of trimmomatic.

The workflow uses a custom pybactocap:1.0.0 container for the processes that run python scripts. If you wish to update the python packages available in this container, please raise an issue in the MRC-CLIMB/containers GitHub repo or put in a pull request modifying https://github.com/MRC-CLIMB/containers/blob/main/pybactocap/Dockerfile.pybactocap-1.0.0

### Conda Recipes ###
The conda recipes for the workflow can be found here: https://github.com/tlforde/OHRBID/tree/main/conda
The software versions are pinned in the recipes, please note that the versions of software differ slightly from the containerised software.
E.g. The latest version of bam-readcount available through conda is 0.8, whereas the bam-readcount container uses the most recent version 1.0.1

### Ignore Errors ###
By default, when a nextflow process fails it will cause the workflow to exit. `errorStrategy 'ignore'` can be added to the process declaration to stop the workflow exiting on the failure of the process (https://www.nextflow.io/docs/latest/process.html#errorstrategy).

### Work Directory ###
Execution of the workflow takes place in the `work` directory. This work directory will grow over time, so it should be deleted as required using `rm`. Alternatively `cleanup = true` can be set in `nextflow.config` (https://www.nextflow.io/docs/latest/config.html#miscellaneous), but this will prevent the use of the resume feature (https://www.nextflow.io/docs/latest/getstarted.html#modify-and-resume)

### Publish Directory ###
In a nextflow process, `publishDir` is used to define which files are published to the outputDir.
E.g.
https://github.com/tlforde/OHRBID/blob/f77e0804da23b674e069e2ffe5f9526874c33588/modules/bactocapModules.nf#L246

### Debug Mode ###
For modules that are dependent on python scripts the following directive has been set
https://github.com/tlforde/OHRBID/blob/f1ec0503ab0d02b33f609b3e6e5b3832a87bdfa3/modules/bactocapModules.nf#L185
This ensures that error messages from the scripts are printed to the stdout

### Executors ###

By default, the pipeline will just run on the local machine. To run on a cluster, modifications will have to be made to the `nextflow.config` to add in the executor. E.g. for a SLURM cluster add `process.executor = 'slurm'`. For more information on executor options see the Nextflow docs: https://www.nextflow.io/docs/latest/executor.html
