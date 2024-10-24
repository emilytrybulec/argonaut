# emilytrybulec/argonaut: Usage

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

To get started running Argonaut, four main files are recommended: [samplesheet.csv](#Samplesheet-input), [params.yaml](#Parameters), [my_config](#Configurations), and [nextflow.sh](#Running-the-pipeline)

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use the input parameter to specify its location in params.yaml. The samplesheet has to be a comma-separated file (.csv) with 3 columns, and a header row. More detail about creating the samplesheets for your specific type of input is located [below](#Using-long-and-short-reads).

### Using long and short reads

The `sample` identifiers are important for naming throughout the pipeline, and we recommend that specific sample names are used, indicating read type. The pipeline requires concatenated raw reads before performing any downstream analysis. For best results, please provide the entire path to the reads.
Below is an example samplesheet containing all read types:

`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2,single_end,read_type
maca_jans_ont,SRR11191910.fastq.gz,,TRUE,ont
maca_jans_ill,SRR11191912_1.fastq.gz,SRR11191912_2.fastq.gz,FALSE,ill
maca_jans_pb,SRR11191909.fastq.gz,,TRUE,pb
```
!!! PLEASE ADD "ont", "pb", AND/OR "ill" TO YOUR SAMPLES NAMES !!!   
Flye and canu will automatically detect your read type based on your samplesheet sample names. Please ensure that read type (ont, pb, or ill) is also indicated in the last column.

|   `sample`    | Custom sample name.  
|   `fastq_1`   | Full path to FastQ file for ONT long reads or Illumina short reads 1. File must have the extension ".fastq" or ".fastq.qz".  
|   `fastq_2`   | Full path to FastQ file for Illumina short reads 2. File must have the extension ".fastq" or ".fastq.gz".  
| `single end`  | True/false indicating whether reads are single end (usually long reads) or paired end (usually short reads).  
| `read_type`   | One of three options indicating the sequencing technology used to generate reads; choice of "ont", "pb", or "ill".  

The [example samplesheet](../assets/samplesheet.csv) is available in a more readable format.


## Parameters
We recommend specifying parameters in a params file. 

The params.yaml file feeds the pipeline all of your input paths! Create a params.yaml file and enter the path to your long and/or short samplesheet(s). You will continue to build on this params file to provide contaminant databases, an out directory name, and more. An example of a params.yaml is shown here:


An example `params.yaml` contains:

```yaml
input                  :  "./samplesheet.csv"
outdir                 :  "test_outdir"
fastq                  :  "./SRR11191910.fastq.gz"
centrifuge_db          :  "./f+b+a+v/"
busco_lineage          :  "metazoa_odb10"
summary_txt            :  "./sequencing_summary.txt"
kraken_db              :  "./kraken_db"
rcf_db                 :  "./recentrifuge/taxdump/"
model                  :  "r1041_e82_400bps_sup_g615"
ragtag_reference       :  "./assembly.fasta"
<...>
```
For best results, please provide full paths. Paths have been truncated for readability.  

For detailed information about parameters, please refer to the [config](../nextflow.config) file.

### General tips for your params file:

* The fastq path should point to one of your read files. If more than one read input is available, input the path to any one of the fastqs here. This parameter is mostly for pipeline initiation and not very important for processing.
* When providing a centrifuge database, please ensure that the path points to a DIRECTORY, not a file.
* Please ensure that your BUSCO lineage aligns with your organism type. (e.g. for japanese walnut tree: "embryophyta_odb10")
* For more information about specific parameters, please refer to [nextflow.config](https://github.com/emilytrybulec/argonaut/blob/main/nextflow.config). It contains explanations for each possible parameter that can be set in the params.yaml file (look for the groups of "=null" labels)

[Here](https://github.com/emilytrybulec/argonaut/blob/main/params.yaml) is a full params.yaml example for a test run.

Not all parameters are required, and the default settings can be modified for individualized use. If you would like to change any settings dictating which assemblers run, whether short reads are available, or options like length filtering and scaffolding, please create a config file using the directions [below](#Configurations).

### Inputting your own assembly  
For users who would like to input an existing genome assembly into the pipeline for downstream processing, please add a line to your params file specifying the path to your *.fasta or *.fa file like so:
```yaml
existing_assembly        :  "./assembly.fasta"
```
Additionally modify your my_config file to indicate that ex_assembly is "true" like so:
```config
params{
  ex_assembly           = true
}
```
Learn more about your my_config file directly below.

## Configurations  

Your my_config file acts as the master switch for controlling the pipeline options. Create a my_config file and populate it with your preferences that differ from the default settings found in [nextflow.config](../nextflow.config). The pipeline accepts fastq files for both short and long reads. ONT and PacBio HiFi reads are considered long read and Illumina reads are considered short read. Please indicate whether your input consists of both or one of the read types in the your my_config file. Specify the full path to your config file with '-c' when [running the nextflow command](#Running-the-pipeline).  

An example `my_config` contains:

```config
params{
  shortread           = false
  min_readlength      = 1000
  kmer_num            = 19
  canu                = true
  ragtag_scaffold     = false
}
```

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run emilytrybulec/argonaut \
  -r main \
  -params-file params.yaml \
  -c my_config \
  -profile singularity,xanadu \
```

This will launch the pipeline with the `xanadu` configuration profile, which will allocate resources correctly for Xanadu users running Argonaut at the University of Connecticut. See [below](#Core-Nextflow-arguments) for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir in params.yaml)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull emilytrybulec/argonaut
```

### Reproducibility

It is a good idea to specify a be aware of versions when running the pipeline on your data. Please refer to the collated_versions.yml file for software versions.

To further assist in reproducbility, you can use share and re-use [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

## Core Nextflow arguments

> **NB:** These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer enviroment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)


### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases you may wish to change which container or conda environment a step of the pipeline uses for a particular tool. By default nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

## Azure Resource Requests

To be used with the `azurebatch` profile by specifying the `-profile azurebatch`.
We recommend providing a compute `params.vm_type` of `Standard_D16_v3` VMs by default but these options can be changed if required.

Note that the choice of VM size depends on your quota and the overall workload during the analysis.
For a thorough list, please refer the [Azure Sizes for virtual machines in Azure](https://docs.microsoft.com/en-us/azure/virtual-machines/sizes).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
