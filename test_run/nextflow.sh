#!/bin/bash
#SBATCH --job-name=ehg_nextflow
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mem=10G
#SBATCH --mail-user=emily.trybulec@uconn.edu

##no xanadu profile, running on mantis general partitions

module load nextflow

nextflow pull emilytrybulec/argonaut
nextflow run emilytrybulec/argonaut -params-file params.yaml -c my_config,mantis_config -profile singularity -r main -resume
