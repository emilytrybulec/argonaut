include { PILON } from '../../modules/nf-core/pilon' 

workflow PHASING {

    take:
        assemblies // channel: tuple(val(meta), path(fasta))
        alignments // channel: tuple(val(meta_bam), path(bam), path(bai)) 
        reference
        pilon_mode

    main:

    ch_versions = Channel.empty()
        println "running phasing with pilon!"

        // combine assemblies and alignments to form pilon_input assemblies
        .combine(alignments) { tuple(meta, fasta), tuple(meta_bam, bam, bai) -> tuple(meta fasta, meta_bam, bam, bai, pilon_mode)
        }
        .set { pilon_input }
        PILON (pilon_input)
        pilon_improved_assemblies = PILON.out.improved_assembly

        pilon_improved_assemblies
            .map { meta, file -> tuple(meta, file) }
            .set {phased_assemblies } 
        
    emit:
    
        phased_assemblies
        
    versions = ch_versions                     // channel: [ versions.yml ]
}
