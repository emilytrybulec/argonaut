include { ALIGN } from '../../modules/local/purgehap' 
include { HISTOGRAM } from '../../modules/local/purgehap' 
include { PURGE } from '../../modules/local/purgehap' 

workflow HAPS {

    take:
        assembly
        reads
        
    main:
    ch_versions = Channel.empty()

        if (params.low == null && params.mid == null && params.high == null) {
        println "generating assembly histogram with purge haplotigs!"

        assembly
            .combine(reads)
            .set{align_ch}

        ALIGN(align_ch, params.bam_format, params.cigar_paf_format, params.cigar_bam)

        assembly
            .join(ALIGN.out.bam)
            .set{assembly_alignment}

        HISTOGRAM(assembly_alignment)
        ch_gencov = HISTOGRAM.out.gencov
        assemblies_polished_purged      = Channel.empty()
        purged_assemblies               = Channel.empty()

        } else if (params.low != null && params.mid != null && params.high != null){
        println "purging assemblies with purge haplotigs!"

        assembly
            .combine(reads)
            .set{align_ch}

        ALIGN(align_ch, params.bam_format, params.cigar_paf_format, params.cigar_bam)

        assembly
            .join(ALIGN.out.bam)
            .set{assembly_alignment}

        HISTOGRAM(assembly_alignment)
        ch_gencov = HISTOGRAM.out.gencov

        assembly
            .join(ch_gencov, by: 0)
            .set{ch_purge}

        PURGE(params.low, params.mid, params.high, ch_purge)
        purged_assemblies      = PURGE.out.purged           
        
        purged_assemblies
                .map { file -> tuple(id: file.baseName, file)  }
                .set { assemblies_polished_purged }
        }

    emit:
        assemblies_polished_purged
        purged_assemblies
        
    versions = ch_versions                     // channel: [ versions.yml ]
}
