include { BIOAWK } from '../../modules/nf-core/bioawk/main'
include { NANOPLOT } from '../../modules/nf-core/nanoplot/main'

workflow LENGTH_FILT {

    take:
  
        decontam_reads  // channel: [ val(meta), [ decontam reads ] ] 
        no_meta_decontam

    main:
    
    ch_versions = Channel.empty()

        // if statement for if min_read_length exists for length filter
        if(params.min_readlength > 0){
            NANOPLOT(decontam_reads)

            BIOAWK(decontam_reads)

            BIOAWK.out.output
                .map { file -> tuple([id:file.baseName, single_end:true], file)  }
                .set { longreads }


            no_meta_longreads = BIOAWK.out.output  // channel: [ val(meta), path(decontam+length filtered fastq) ]
        } else{
            decontam_reads  // channel: [ val(meta), path(decontaminated fastq) ]
                .set{longreads}
            no_meta_decontam
                .view()
                .set{no_meta_longreads}
        }

    emit:
        longreads
        no_meta_longreads

    versions = ch_versions                     // channel: [ versions.yml ]
}
