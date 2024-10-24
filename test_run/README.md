## Testing Argonaut

All scripts required to run Argonaut are located within this directory, and the test data can be found in the [v.1.0.2](https://github.com/emilytrybulec/argonaut/releases/tag/v.1.0.2) release as Assets. 

Eastern hoolock gibbon test data was sourced from [GenomeArk](https://www.genomeark.org/genomeark-all/Hoolock_leuconedys.html). The publicly available ONT, PacBio, and Illumina reads have been trimmed down to only contain chromosome 3 hits.

UConn users who would like to run Argonaut on Mantis (much faster than Xanadu!) may copy and paste the [mantis_config](https://github.com/emilytrybulec/argonaut/blob/main/test_run/mantis_config) file into their working directory along with other scripts, indicate that there is a second config file in nextflow.sh, and remove the xanadu profile as shown in the [test nextflow.sh](https://github.com/emilytrybulec/argonaut/blob/v.1.0.2/test_run/nextflow.sh).

Happy assembling!
