include { SEQKIT_SEQ } from '../../modules/local/seqkit'
include { NANOPLOT } from '../../modules/nf-core/nanoplot/main'

workflow LENGTH_FILT3 {

    take:
  
        decontam_reads  // channel: [ val(meta), [ decontam reads ] ] 
        decontam_reads_no_meta
    main:
    
    ch_versions = Channel.empty()

        // if statement for if min_read_length exists for length filter
        if(params.min_readlength > 0){
            NANOPLOT(decontam_reads)

            SEQKIT_SEQ(decontam_reads, params.min_readlength)

            SEQKIT_SEQ.out.filter
                .map { file -> tuple([id:file.simpleName, single_end:true], file)  }
                .set { longreads }


            no_meta_longreads = SEQKIT_SEQ.out.filter  // channel: [ val(meta), path(decontam+length filtered fastq) ]
        }
        else{
            longreads = decontam_reads  // channel: [ val(meta), path(decontaminated fastq) ]
            no_meta_longreads = decontam_reads_no_meta
        }

    emit:
        longreads
        no_meta_longreads

    versions = ch_versions                     // channel: [ versions.yml ]
}
