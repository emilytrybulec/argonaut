[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Nextflow Tower](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Nextflow%20Tower-%234256e7)](https://tower.nf/launch?pipeline=https://github.com/nf-core/genomeassembly)

## Introduction


**Argonaut** performs **a**utomated **r**eads to **g**enome **o**perations for de **n**ovo **a**ssemblies; it is a bioinformatics pipeline that performs genome assembly on long and short read data. A fastq file and input information is fed to the pipeline, resulting in final assemblies with completeness, contiguity, and correctnesss quality checking at each step. The pipeline accepts short reads, long reads, or both. 

<img align="right" height="300" src="https://github.com/emilytrybulec/genomeassembly/assets/114685119/9b900dab-44cb-479e-9362-0c0d9dc00ae0">

## Table of Contents
- [Pipeline Summary](#Pipeline-Summary)
- [Quick Start](#Quick-Start)
- [Output Overview](#Pipeline-Output)
- [Credits](#Credits)
- [Contributions & Support](#Contributions-and-Support)
- [Citations](#Citations)
   
## Pipeline Summary

Illumina Short Read 
1. Read QC, Adaptor Trimming, Contaminant Filtering([`FastQC v0.11.9`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [`FastP v0.23.4`](https://github.com/OpenGene/fastp), [`GenomeScope2 v2.0`](http://qb.cshl.edu/genomescope/),[`Jellyfish v2.2.6`](https://github.com/gmarcais/Jellyfish),[`Kraken2 v2.1.2`](https://ccb.jhu.edu/software/kraken2/), [`Recentrifuge v1.9.1`](https://github.com/khyox/recentrifuge),)
    
PacBio HiFi Long Read (CCS format)
1. Read QC, Adaptor Trimming, Contaminant Filtering([`Nanoplot v1.41.0`](https://github.com/wdecoster/NanoPlot),[`CutAdapt v3.4`](https://cutadapt.readthedocs.io/en/stable/),[`GenomeScope2 v2.0`](http://qb.cshl.edu/genomescope/),[`Jellyfish v2.2.6`](https://github.com/gmarcais/Jellyfish),[`Kraken2 v2.1.2`](https://ccb.jhu.edu/software/kraken2/), [`Recentrifuge v1.9.1`](https://github.com/khyox/recentrifuge))
2. Length Filtering (optional)([`Seqkit v2.4.0`](https://bioinf.shenwei.me/seqkit/usage/#seq), [`Nanoplot v1.41.0`](https://github.com/wdecoster/NanoPlot))

  
ONT Long Read
1. Read QC and Contaminant Filtering([`Nanoplot v1.41.0`](https://github.com/wdecoster/NanoPlot),[`KmerFreq`](https://github.com/fanagislab/kmerfreq), [`GCE`](https://github.com/fanagislab/GCE), [`Centrifuge v1.0.4`](https://ccb.jhu.edu/software/centrifuge/), [`Recentrifuge v1.9.1`](https://github.com/khyox/recentrifuge))
2. Length Filtering (optional)([`Seqkit v2.4.0`](https://bioinf.shenwei.me/seqkit/usage/#seq), [`Nanoplot v1.41.0`](https://github.com/wdecoster/NanoPlot))
   

All reads are used for the following steps:  
    <img align="right" width="600" alt="Argonaut Hybrid Workflow" src="https://github.com/emilytrybulec/argonaut/assets/114685119/a746d16c-5197-486b-b5cd-78c3ab54d840">  

3. Assembly 
- [`Flye v2.9`](https://github.com/fenderglass/Flye) 
- [`Canu v2.2`](https://github.com/marbl/canu)
- [`Verkko v2.2`](https://github.com/marbl/verkko)
- [`Hifiasm v0.19.8`](https://github.com/chhylp123/hifiasm) 
- [`MaSuRCA v4.1.0`](https://github.com/alekseyzimin/masurca) 
- [`Redundans v2.01`](https://github.com/Gabaldonlab/redundans) 

4. Polish  
- [`Medaka v1.8.0`](https://github.com/nanoporetech/medaka)
- [`Racon v1.4.20`](https://github.com/isovic/racon)
- [`POLCA v4.1.0`](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007981)
- [`Pilon v1.23.0`](https://github.com/broadinstitute/pilon)
  
5. Purge 
- [`PurgeHaplotigs v1.1.2`](https://bitbucket.org/mroachawri/purge_haplotigs/src/master/)
- [`Redundans v2.01`](https://github.com/Gabaldonlab/redundans)
  
6. Scaffolding 
- [`RagTag v2.1.0`](https://github.com/malonge/RagTag)
  
7. Quality Checking
- [`BUSCO v5.7.1`](https://busco.ezlab.org/) or [`Compleasm v0.2.6`](https://github.com/huangnengCSU/compleasm)
- [`Quast v5.2.0`](https://quast.sourceforge.net/)
- [`Merqury v1.3`](https://github.com/marbl/merqury)
- [`Minimap2 v2.24`](https://github.com/lh3/minimap2) or [`Winnowmap v2.03`](https://github.com/marbl/Winnowmap)

8. Assembly Visualization
- ([`Blobtools v4.3.9`](https://blobtoolkit.genomehubs.org/blobtools2/))  

To the right is a figure detailing the major workflow steps involved in hybrid assembly.
  
If you indicate that you would like for long read polishers to be run, the pipeline will default to using PacBio HiFi reads, and using ONT if no PacBio HiFi is available. If short reads are also available, they will automatically be used to polish the assemblies after long read polishing (or assembly if long read polishing is off).   
  
Purge Haplotigs is the first step of manual curation, as it produces a histogram that needs to be analyzed for -l, -m, -h flags. The pipeline will stop at the purge step if purge is activated in your configuration and wait for manual input of parameters according to the histogram of your assembly, which can be found in your out directory.

## Quick Start
**Installation**
> Only Nextflow and Singularity need to be installed to run Argonaut. Users that would like to run [Centrifuge](https://ccb.jhu.edu/software/centrifuge/manual.shtml#database-download-and-index-building) and/or [Kraken2](https://ccb.jhu.edu/software/kraken2/index.shtml?t=downloads) will need to provide a database. There are similar restrictions for running [Recentrifuge](https://github.com/khyox/recentrifuge/wiki/Installation#3-getting-the-databases) and Blobtools with [Blast](https://www.nlm.nih.gov/ncbi/workshops/2023-08_BLAST_evol/databases.html) and [NCBI](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/) taxdump. Follow the links provided for database download directions. Xanadu users running Argonaut at the University of Connecticut may use the database paths provided in the example [params.yaml](https://github.com/emilytrybulec/argonaut/blob/main/params.yaml)

**Samplesheets**  

To get started setting up your run, prepare a samplesheet with your input data as follows:

`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2,single_end,read_type
chr3_gibbon_pb,/core/projects/EBP/conservation/gen_assembly_pipeline/hoolock_chrm_3/chr3_pb.fastq.gz,,TRUE,pb
chr3_gibbon_ont,/core/projects/EBP/conservation/gen_assembly_pipeline/hoolock_chrm_3/chr3_ont.fastq.gz,,TRUE,ont
chr3_gibbon_ill,/core/projects/EBP/conservation/gen_assembly_pipeline/hoolock_chrm_3/chr3_ill_R1.paired.fastq.gz,/core/projects/EBP/conservation/gen_assembly_pipeline/hoolock_chrm_3/chr3_ill_R2.paired.fastq.gz,FALSE,ill
```

!!! PLEASE ADD "ont", "pb", AND/OR "ill" TO YOUR SAMPLE NAMES AND FILE NAMES!!! Failure to do so may result in assemblers not recognizing your read type and/or outputs being overwritten.

The sample name inputted in your samplesheet will serve as the prefix for your output files. Please indicate which kind of read is being inputted in the sample name, as well as the read type column. 



After you have your samplesheet, create a params.yaml file to specify the paths to your samplesheet, contaminant databases, etc. Most likely, a config file will also need to be made to modify the default settings of the pipeline. Please look through the [nextflow.config](nextflow.config) file to browse the defaults and specify which you would like to change in your my_config file. More information is located in the [usage](docs/usage.md) section.

Now, you can run the pipeline using:


```bash
nextflow run emilytrybulec/argonaut \
  -r main \
  -params-file params.yaml \
  -c my_config \
  -profile singularity,xanadu \
```

## Pipeline output
All of the output from the programs run in the pipeline pipeline will be located in the out directory specified in params.yaml. The pipeline produces the following labeled directories depending on configurations:

```
├── 01 READ QC
│   ├── centrifuge
│   ├── fastp
│   ├── fastqc
│   ├── genome size est
│   │   ├── genomescope2
│   │   ├── jellyfish
│   │   ├── ont gce 
│   │   ├── ont kmerfreq
│   ├── kraken2
│   ├── nanoplot
│   ├── pacbio cutadapt
├── 02 ASSEMBLY
│   ├── hybrid
│   ├── long read
│   ├── short read
├── 03 POLISH
│   ├── hybrid
│   │   ├── polca
│   ├── long read
│   │   ├── medaka
│   │   ├── racon
├── 04 PURGE
│   ├── align
│   ├── histogram
│   ├── purge haplotigs
│   ├── short read redundans
├── 05 SCAFFOLD
├── ASSEMBLY QC
│   ├── busco
│   ├── bwamem2
│   ├── merqury
│   ├── minimap2
│   ├── quast
│   ├── samtools
├── OUTPUT
│   ├── blobtools visualization
│   ├── coverage
│   ├── genome size estimation
│   ├── *assemblyStats.txt
├── PIPELINE INFO
    └── execution_trace_*.txt
```
Some output files have labels such as "dc", indicating that the reads have been decontaminated, or  "lf", indicating that reads have been length filtered.  
  
Information about interpreting output is located in the [output](docs/output.md) section.

## Credits
emilytrybulec/genomeassembly was originally written by Emily Trybulec.

I thank the following people for their extensive assistance in the development of this pipeline:

University of Connecticut:  
<img align="right" height="210" src="https://github.com/emilytrybulec/argonaut/assets/114685119/7a3fd47c-0fbf-443c-a121-9fd8a3da9ba3">

* Biodiversity and Conservation Genomics Center  
     * Jill Wegrzyn  
     * Cynthia Webster  
     * Anthony He  
     * Laurel Humphrey  
     * Keertana Chagari  
     * Amanda Mueller  
     * Cristopher Guzman  
     * Harshita Akella
  
<img align="right" height="140" src="https://github.com/emilytrybulec/argonaut/assets/114685119/91c25e9f-f70b-481f-8aab-55d2d529eca4">

* Rachel O'Neill Lab  
     * Rachel O’Neill  
     * Michelle Neitzey  
     * Nicole Pauloski  
     * Vel Johnston
  
<img align="right" height="150" src="https://github.com/emilytrybulec/argonaut/assets/114685119/161c0c34-4f05-496d-9436-2d087ba5ccd1">  

* Computational Biology Core  
     * Noah Reid  
     * Gabe Barrett  

* nf-core Community  

* Zbigniew Trybulec

 
## Contributions and Support

Development of this pipeline was funded by the University of Connecticut Office of Undergraduate Research through the Summer Undergraduate Research Fund (SURF) Grant.

<img align="right" height="110" src="https://github.com/emilytrybulec/argonaut/assets/114685119/925bab1e-ba82-44d3-b640-3d7cf6c2028f">  
<img align="right" height="110" src="https://github.com/emilytrybulec/argonaut/assets/114685119/3276ed4c-e272-4ff5-93ba-4a7802bd78aa">  

The Biodiversity and Conservation Genomics Center is a part of the [Earth Biogenome Project](https://www.earthbiogenome.org/), working towards capturing the genetic diversity of life on Earth.

## Citations
Argonaut is currently unpublished. For now, please use the GitHub URL when referencing.

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
