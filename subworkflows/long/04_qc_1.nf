include { QUAST } from '../../modules/local/quast'  
include { BUSCO } from '../../modules/nf-core/busco/main' 
include { PYCOQC } from '../../modules/nf-core/pycoqc/main'  
include { MINIMAP2_INDEX } from '../../modules/nf-core/minimap2/index/main' 
include { MINIMAP2_ALIGN } from '../../modules/nf-core/minimap2/align/main'  
include { MERYL_COUNT } from '../../modules/nf-core/meryl/count/main' 
include { MERQURY } from '../../modules/nf-core/merqury/main' 
include { SAMTOOLS_INDEX } from '../../modules/nf-core/samtools/index/main' 
include { BWAMEM2_INDEX } from '../../modules/nf-core/bwamem2/index/main' 
include { BWAMEM2_MEM } from '../../modules/nf-core/bwamem2/mem/main' 

workflow QC_1 {

    take:
        assemblies // channel: [ val(meta), path(assembly.fasta) ]
        fastq_filt // channel: [ val(meta), path(filtered long reads) ]
        summarytxt // channel from params.summarytxt
        shortreads
        genome_size_est
        flattened_lr

    main:

    ch_versions = Channel.empty() 
    ch_index = Channel.empty()
    ch_align_bam = Channel.empty()
    ch_align_paf = Channel.empty()
    ch_combo = Channel.empty()

        if ( params.shortread == true ) {
            BWAMEM2_INDEX(assemblies)
            BWAMEM2_MEM(shortreads, BWAMEM2_INDEX.out.index) }

        if ( params.longread == true ){
            // build index
            MINIMAP2_INDEX(assemblies)
            ch_versions = ch_versions.mix(MINIMAP2_INDEX.out.versions)
            ch_index = MINIMAP2_INDEX.out.index

            // align reads
            MINIMAP2_ALIGN(fastq_filt, assemblies, params.bam_format, params.cigar_paf_format, params.cigar_bam)
            ch_align_bam = MINIMAP2_ALIGN.out.bam
            ch_align_paf = MINIMAP2_ALIGN.out.paf

            fastq_filt
                .map { file -> file }
                .flatten()
                .set { fastq_no_meta }

            assemblies
                .combine(fastq_no_meta)
                .set{ch_combo}
        } else { ch_combo = Channel.empty() }

        // run quast
        QUAST( assemblies // this has to be aggregated because of how QUAST makes the output directory for reporting stats)
        ch_quast = QUAST.out.results
        ch_versions = ch_versions.mix(QUAST.out.versions)

        // run BUSCO
        BUSCO(assemblies, params.busco_lineage, [], [])
        ch_busco = BUSCO.out.short_summaries_txt
        ch_busco_full_table = BUSCO.out.full_table
        ch_versions = ch_versions.mix(BUSCO.out.versions)

        if ( params.longread == true ){
            SAMTOOLS_INDEX (MINIMAP2_ALIGN.out.bam)
            ch_sam = SAMTOOLS_INDEX.out.sam

            ch_combo
                .join(ch_sam)
                .set{racon}

        } else if ( params.shortread == true ){ 
            SAMTOOLS_INDEX (BWAMEM2_MEM.out.bam)
            ch_sam = SAMTOOLS_INDEX.out.sam
            racon = Channel.empty() }

    
        if ( params.summary_txt_file == true ) {
            // create summary txt channel with meta id and run pycoQC
            ch_summarytxt = summarytxt.map { file -> tuple(file.baseName, file) }

            PYCOQC (ch_summarytxt, MINIMAP2_ALIGN.out.bam, SAMTOOLS_INDEX.out.bai)
            ch_versions = ch_versions.mix(PYCOQC.out.versions)

        } else { ch_summarytxt = Channel.empty() }

        if ( params.shortread == true ) {
            MERYL_COUNT ( shortreads, params.kmer_num ) 
        } else {
            MERYL_COUNT ( fastq_filt, params.kmer_num ) }

        MERQURY ( assemblies, MERYL_COUNT.out.meryl_db, genome_size_est, params.tolerable_collision  )
        ch_merqury = MERQURY.out.assembly_qv

    emit:
        ch_index
        ch_align_bam
        ch_align_paf
        ch_quast
        ch_busco
        ch_merqury
        ch_summarytxt
        ch_meryl = MERYL_COUNT.out.meryl_db
        ch_sam
        ch_busco_full_table
        racon

        
    versions = ch_versions                     // channel: [ versions.yml ]
}
