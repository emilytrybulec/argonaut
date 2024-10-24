# emilytrybulec/argonaut: Output

## Introduction

This document describes the output produced by the pipeline. All of the programs being run will have their own folder in the output directory, nested within a folder describing which step each program is used for.

The directories created will depend on which options are selected in the configuration and which programs are run. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

* LONG READ QC  
  * [Centrifuge](#centrifuge) - Contaminant detection for long reads
    * [Recentrifuge](#recentrifuge) - Contaminant filtering visualization for short and long reads
    * [Contaminant filter](#seqkit) - Filtering contaminants detected by Centrifuge from long reads for downstream processing  
  * Genome Size Estimation - Genome size estimation using kmerfreq and gce  
    * [Kmerfreq](#kmerfreq)  
    * [GCE](#gce)  
  * [Length Filter](#lengthfilter) - Optional length filtering with bioawk based on minimum read length parameter  
  * [Nanoplot](#nanoplot) - Long read quality checking  
  
* SHORT READ QC  
  * [Fastp](#fastp) - Adapter trimming for short reads  
  * [Fastqc](#fastqc) - Short read quality checking  
  * [Genomescope2](#genomescope2) - Short read ploidy and genome size estimatation  
  * [Jellyfish](#jellyfish) - Short read kmer distribution  
  * [Kraken2](#kraken2) - Contaminant detection for short reads  
    * [Recentrifuge](#recentrifuge) - Contaminant filtering visualization for short and long reads  
  
* ASSEMBLY
  * hybrid
    * [MaSuRCA](#masurca) - Masurca hybrid assembly
  * long read
    * [Canu](#canu) - Canu assembly (for high coverage data) 
    * [Flye](#flye) - Flye assembly
    * [Hifiasm](#hifiasm) - Hifiasm asssembly (for PacBio Hifi data)
  * short read
    * [MaSuRCA](#masurca) - Masurca short read assembly
    * [Redundans](#masurca) - Redundans assembly
  
* POLISH  
  * [Medaka](#medaka) - Long read polishing assemblies
  * [Racon](#racon) - Long read polishing assemblies
  * [POLCA](#polca) - Short (or long) read polishing assemblies
  
* PURGE  
  * [Align](#purge) - Alignment of raw reads to assemblies  
  * [Histogram](#purge) - Histogram of read-depth vs count  
  * [Purge](#purge) - Purged assembly using purge haplotigs to reduce duplication  
  
* ASSEMBLY QC  
  * [Busco](#busco) - Assembly quality checking for completeness  
  * [Merqury](#merqury) - Assembly quality checking for accuracy  
    * [Meryl](#meryl) - Building a database for merqury quality checking  
  * [Minimap2](#minimap2) - Assemblies with aligned reads  
  * [PycoQC](#pycoqc) - Assembly quality checking with sequencing summary  
  * [Quast](#quast) - Assembly quality checking for contiguity  
  * [Samtools](#samtools) - Indexed Minimap2 alignments  
   
* OTHER  
  * [Gunzip](#gunzip/gzip) - Converting files from .gz to unzipped (no .gz)  
  * [Gzip](#gunzip/gzip) - Converting files from unzipped to .gz  
  
* OUTPUT
  * [Genome Size Estimation](#extract) - Genome size estimation from gce OR manually inputted in params.yaml
  * [Long Read Coverage](#coverage) - Coverage calculated from long reads and genome size estimation
  * [Short Read Coverage](#coverage) - Coverage calculated from short reads and genome size estimation
  * [Output](#output) - Final quality stats (quast, busco, merqury) of all assemblies produced
  
* [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution  



### busco

<details markdown="1">
<summary>Output files</summary>

- `busco/`
  - `canu_*/`: Canu assembly BUSCO output directory
  - `flye_*/`: Flye assembly BUSCO output directory
  - `masurca_*/`: MaSuRCA assembly BUSCO output directory
  - `polca_*/`: POLCA assembly BUSCO output directory

</details>

[BUSCO](https://busco.ezlab.org/) gives completeness quality metrics about your assemblies. It provides information about whether your assembly matches with known protein sequences for similar species. Please be sure to specify the correct BUSCO lineage for your species. For further reading and documentation see the [BUSCO user guide](https://busco.ezlab.org/busco_userguide.html).

### canu

### centrifuge

<details markdown="1">
<summary>Output files</summary>

- `centrifuge/`  
  - `centrifuge/`: Directory containting contaminant detection data  
    -`filtered_*.report.txt/`  
    -`filtered_*.results.txt/`  
    -`filtered_*.unmapped.fastq.gz/`: Centrifuge filtered reads (not used for downstream processing, as filtering with seqkit is more accurate)
  - `filtered/`: Directory containting decontaminated reads  
  - `kreport/`: Directory containing a kraken style text report about the contaminants detected
  - `recentrifuge/`: Directory containting user friendly html visualization of contaminants detected

</details>

[Centrifuge](https://ccb.jhu.edu/software/centrifuge/) is a contaminant detection tool that generates a report summarizing which of your reads contain known bacterial, fungal, viral, etc. sequences. These reads are flagged and filtered from the reads used downstream, as they are likely contaminants. Please provide a centrifuge database containing contaminant sequences in params.yaml for the pipeline to run smoothly. For more information, please see the [Centrifuge manual](https://github.com/DaehwanKimLab/centrifuge/blob/master/MANUAL).

[Recentrifuge](https://github.com/khyox/recentrifuge) is an interactive visualization tool for contaminant detection. An html report is generated that summarizes the results of centrifuge or kraken. For more information, see the [Recentrifuge wiki](https://github.com/khyox/recentrifuge/wiki).

### coverage

<details markdown="1">
<summary>Output files</summary>

- `*coverage/`
  - `estimatedCoverage_after.txt`: Estimated coverage after read quality filtering
  - `estimatedCoverage_before.txt`: Estimated coverage before read quality filtering
  - `totalBases_after.txt`: Number of total bases after quality filtering
  - `totalBases_before.txt`: Number of total bases before quality filtering

</details>

The coverage estimation module uses linux commands to extract the total number of bases before and after quality checking and divide total bases over the estimated genome size.

### extract

<details markdown="1">
<summary>Output files</summary>

- `extract/`
  - `finalSize.txt`: Genome size estimation in megabases or gigabases
  - `standardSize.txt`: Genome size estimation in standard numerical format

</details>

The extract module uses awk and numfmt to isolate the genome size estimate for downstream use in assemblers like flye and canu, as well as quality checking with merqury.

### fastp

<details markdown="1">
<summary>Output files</summary>

- `fastp/`
  - `*1.fail.fastq.gz`
  - `*1.fastp.fastq.gz`
  - `*2.fail.fastq.gz`
  - `*2.fastp.fastq.gz`
  - `*fastp.html`
  - `*fastp.json`
  - `*fastp.log`
  - `*merged.fastq.gz`

</details>

[Fastp](https://github.com/OpenGene/fastp) does short read adapter trimming and provides quality metrics about your short reads. For further reading and documentation see the [Fastp usage documentation](https://github.com/OpenGene/fastp#simple-usage).

### fastqc

<details markdown="1">
<summary>Output files</summary>

- `fastqc/`
  - `2_adaper_trim/`
    - `*1.fastqc.html`
    - `*1.fastqc.zip`
    - `*2.fastqc.html`
    - `*2.fastqc.zip`
  - `3_decontam/`
    - `*1.fastqc.html`
    - `*1.fastqc.zip`
    - `*2.fastqc.html`
    - `*2.fastqc.zip`
  - `*1.fastqc.html`
  - `*1.fastqc.zip`
  - `*2.fastqc.html`
  - `*2.fastqc.zip`

</details>

[Fastqc](https://github.com/s-andrews/FastQC) provides quality metrics about your short reads. It is run on the raw short reads, after adapter trimming with fastp, and after contaminant filtering with kraken2. For further reading and documentation see the [Fastqc help page](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

### flye

<details markdown="1">
<summary>Output files</summary>

- `flye/`
  - `assembly_graph.gfa`
  - `assembly_graph.gv`
  - `flye_*.assembly_info.txt`
  - `flye_*.fasta`
  - `flye_*.flye.log`
  - `flye_*.params.json`

</details>

[Flye](https://github.com/fenderglass/Flye) is a long read assembler. The fasta file is used downstream. For further reading and documentation see the [Flye usage documentation](https://github.com/fenderglass/Flye/blob/flye/docs/USAGE.md).

### gce

<details markdown="1">
<summary>Output files</summary>

- `gce/`
  - `gce2.log`

</details>

[GCE](https://github.com/fanagislab/GCE) is a genomic character estimator. It uses the kmerfreq output to estimate genome size. 

### genomescope2

<details markdown="1">
<summary>Output files</summary>

- `genomescope2/`
  - `*linear_plot.png`
  - `*log_plot.png`
  - `*model.txt`
  - `*summary.txt`
  - `*transformed_linear_plot.png`
  - `*transformed_log_plot.png`

</details>

[Genomescope2](https://github.com/tbenavi1/genomescope2.0) gives quality metrics about your short reads. It provides information about ploidy, repetitive content, genome size, heterozygosity, and sequencing coverage. This program will only run if short reads are provided.

### gunzip/gzip

<details markdown="1">
<summary>Output files</summary>

- `gunzip/`
  - `*.fastq`
- `gzip/`
  - `*.fastq.gz`

</details>

[Gzip](https://www.gzip.org/) is a file compression tool that may be used in the pipeline to convert files from fastq to fastq.gz and vice versa. For further reading and documentation see the [Gzip documentation](https://www.gnu.org/software/gzip/manual/gzip.html).

### hifiasm

### jellyfish

<details markdown="1">
<summary>Output files</summary>

- `jellyfish/`
  - `*_mer_out.histo`
  - `*_mer_out.jf`


</details>

[Jellyfish](https://github.com/gmarcais/Jellyfish) performs kmer counting to generate histograms used in Genomescope2 for genome size estimation. This program will only run if short reads are provided. For further reading and documentation see the [Jellyfish documentation](https://github.com/gmarcais/Jellyfish/blob/master/doc/Readme.md).

### kmerfreq

<details markdown="1">
<summary>Output files</summary>

- `kmerfreq/`
  - `kmernum.txt`
  - `read_files.lib.kmer.freq.stat.2column`


</details>

[Kmerfreq](https://github.com/fanagislab/kmerfreq) performs kmer counting for long reads to generate tables used in GCE for genome size estimation.

### kraken2

<details markdown="1">
<summary>Output files</summary>

- `kraken2/`
  - `*.classified_1.fastq.gz`
  - `*.classified_2.fastq.gz`
  - `*.kraken2.classifiedreads.txt`
  - `*.kraken2.report.txt`
  - `*.unclassified_1.fastq.gz`
  - `*.unclassified_2.fastq.gz`
- `recentrifuge/`: Directory containting user friendly html visualization of contaminants detected

</details>

[Kraken2](https://github.com/DerrickWood/kraken2/tree/master) is used to detect contaminant sequences in short reads. Reads flagged as contaminants are filtered from the reads used downstream. Please provide a kraken2 database containing contaminant sequences in params.yaml for the pipeline to run smoothly. For more information, please see the [Kraken2 user manual](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown).

[Recentrifuge](https://github.com/khyox/recentrifuge) is an interactive visualization tool for contaminant detection. An html report is generated that summarizes the results of centrifuge or kraken. For more information, see the [Recentrifuge wiki](https://github.com/khyox/recentrifuge/wiki).

### masurca

<details markdown="1">
<summary>Output files</summary>

- `masurca/`
  - `masurca_*.fasta`


</details>

[Masurca](https://github.com/fenderglass/Flye) is a hybrid assembler which will only run if both long and short reads are provided. The fasta file is used downstream.

### medaka

<details markdown="1">
<summary>Output files</summary>

- `medaka/`
  - `*.fasta`: flye assembly polished with long reads


</details>

[Medaka](https://github.com/nanoporetech/medaka) is a tool that polishes assemblies based on consensus sequences to reduce errors. Please provide the model of flow cell that was used for sequencing in params.yaml for best results. The resulting fasta file is used downstream.

### merqury

<details markdown="1">
<summary>Output files</summary>

- `merqury/`
  - `best_kmer_num.txt`: estimated best kmer number from best_k.sh
  - `*.qv`: stats including merqury quality value

</details>

[Merqury](https://github.com/marbl/merqury) gives accuracy quality metrics about your assemblies. It provides information about whether your assembly matches with kmers present in your long reads (or short reads if available) generated with Meryl. For further reading about interpreting output, see the [Merqury guide](https://github.com/marbl/merqury/wiki/2.-Overall-k-mer-evaluation).

### meryl

<details markdown="1">
<summary>Output files</summary>

- `meryl/`
  - `kmer_db.filtered.meryl`: directory containing meryl indeces and data

</details>

[Meryl](https://github.com/marbl/meryl) creates a database containing raw read kmers for use in Merqury. Meryl will run on the long reads unless short reads are provided.

### minimap2

<details markdown="1">
<summary>Output files</summary>

- `minimap2/`
  - `*.mmi`: mapped alignment
  - `*.bam`: reads aligned to the assembly

</details>

[Minimap2](https://github.com/lh3/minimap2/tree/master) is a sequence alignment program that will align your long reads to your assemblies. For a detailed description of usage, please refer to the [Minimap2 Manual](https://lh3.github.io/minimap2/minimap2.html)

### nanoplot

<details markdown="1">
<summary>Output files</summary>

- `nanoplot/`
  - `*/`: output directory before contaminant filtering with centrifuge
  - `*_filtered/`: output directory after contaminant filtering with centrifuge

</details>

[Nanoplot](https://github.com/wdecoster/NanoPlot) is a quality checking program for ONT long reads that provides visualizations of sequencing data and general statistics. Nanoplot is run on the reads before and after filtering centrifuge contaminants.

### output

<details markdown="1">
<summary>Output files</summary>

- `output/`
  - `*assemblyStats.txt`: assembly statistics from quast, busco, and merqury

</details>

The output module pulls quality checking statistics from assemblies at each step of the pipeline for simpler interpretation of results. 

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.

### polca

<details markdown="1">
<summary>Output files</summary>

- `polca/`
  - `polca_*.fasta`: flye assembly polished with short reads

</details>

[POLCA](https://github.com/alekseyzimin/masurca#polca) a polishing tool built into MaSuRCA. This tool can be used separate of the hybrid assembler and will run in the pipeline if short reads are available, polishing the flye assembly with the Illumina reads. The resulting fasta file is used downstream.

### purge

<details markdown="1">
<summary>Output files</summary>

- `purge/`
  - `align/`: directory containing read alignments to assemblies
  - `histogram/`: directory containing purge histogram
  - `purged/`: directory containing purged assemblies

</details>

[Purge Haplotigs](https://bitbucket.org/mroachawri/purge_haplotigs/src/master/) is a tool for reducing duplication in assemblies. Purge will initially only produce a histogram. Examine the histogram for bimodal peaks and determine the low, mid, and high cutoffs. Specify the cutoffs in my_config file, along with the path to the gencov file in the histogram directory. Once the parameters are updated, resume the pipeline to allow purge to finish running. The resulting fasta files in the purged directory are used downstream. For more information about determining cutoffs, see the [Purge usage tutorial](https://bitbucket.org/mroachawri/purge_haplotigs/wiki/Tutorial).

### pycoqc

<details markdown="1">
<summary>Output files</summary>

- `pycoqc/`
  - `flye_*/`: Flye assembly PycoQC output directory
  - `polca_*/`: Flye assembly polished with POLCA PycoQC output directory

</details>

[PycoQC](https://github.com/a-slide/pycoQC) generates interactive quality checking plots for long reads, in relation to assemblies in this case. The "Bases alignment rate breakdown" figure provides information about how many of the reads are present in the assembly and how many have been clipped. Please be sure to specify the path to your sequencing summary text file in params.yaml for the pipeline to run smoothly.

### quast

<details markdown="1">
<summary>Output files</summary>

- `quast/`
  - `canu_*/`: Canu assembly quast output directory
  - `flye_*/`: Flye assembly quast output directory
  - `masurca_*/`: MaSuRCA assembly quast output directory
  - `polca_*/`: POLCA assembly quast output directory
  - `purge_*/`: Purged assembly quast output directory

</details>

[Quast](https://github.com/a-slide/pycoQC) is an quality checking tool that assesses contiguity. A report is generated that summarizes the statistics for each assembly. For more information, see the [Quast user manual](https://quast.sourceforge.net/docs/manual.html).

### racon

### recentrifuge

<details markdown="1">
<summary>Output files</summary>

- `recentrifuge/`
  - `centrifuge/`: Visualization of centrifuge (long read) detected contaminants
  - `kraken/`: Visualization of kraken (short read) detected contaminants


</details>

[Recentrifuge](https://github.com/khyox/recentrifuge) is an interactive visualization tool for contaminant detection. An html report is generated that summarizes the results of centrifuge or kraken. For more information, see the [Recentrifuge wiki](https://github.com/khyox/recentrifuge/wiki).

### samtools

<details markdown="1">
<summary>Output files</summary>

- `samtools/`
  - `*bam.bai`: reads aligned and indexed to the assemblies

</details>

[Samtools](https://github.com/samtools/samtools) is a tool for indexing aligned files. Output from Minimap2 is indexed with samtools.

### seqkit

<details markdown="1">
<summary>Output files</summary>

- `filtered/`
  - `*_filtered.fastq`: Filtered long reads after Centrifuge contaminant detection

</details>

[Seqkit](https://github.com/shenwei356/seqkit) is a tool for file manipulation and is used in the pipeline to filter flagged contaminants from long reads. For further reading and documentation see the [seqkit user guide](https://bioinf.shenwei.me/seqkit/usage/).
