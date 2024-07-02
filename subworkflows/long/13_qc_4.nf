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

workflow QC_4 {

    take:
        assemblies // channel: [ val(meta), path(assemblies) ]
        fastq_filt // channel: [ val(meta), path(filtered reads) ]
        summarytxt // channel from params.summarytxt
        ch_quast
        ch_busco
        ch_merqury
        shortreads
        genome_size_est
        ch_meryl

    main:

    ch_versions = Channel.empty() 

    if ( params.shortread == true ) {
        BWAMEM2_INDEX(polished_assemblies)
        BWAMEM2_MEM(shortreads, BWAMEM2_INDEX.out.index)
    }

    if ( params.longread == true ){
        // build index
        MINIMAP2_INDEX(assemblies)
        ch_versions = ch_versions.mix(MINIMAP2_INDEX.out.versions)
        ch_index = MINIMAP2_INDEX.out.index

        // align reads
        MINIMAP2_ALIGN(fastq_filt, assemblies, params.bam_format, params.cigar_paf_format, params.cigar_bam)
        ch_bam = MINIMAP2_ALIGN.out.bam}
        
        // run quast
        QUAST(
            assemblies
        )
        ch_quast
            .concat(QUAST.out.results)
            .set { ch_quast }
        ch_versions = ch_versions.mix(QUAST.out.versions)

        // run BUSCO
        BUSCO(assemblies, params.busco_lineage, [], [])
        ch_busco
            .concat(BUSCO.out.short_summaries_txt)
            .set { ch_busco }

        ch_busco_full_table = BUSCO.out.full_table

        ch_versions = ch_versions.mix(BUSCO.out.versions)

    if ( params.longread == true ){
        SAMTOOLS_INDEX (MINIMAP2_ALIGN.out.bam)
        ch_sam = SAMTOOLS_INDEX.out.sam
    } else if ( params.shortread == true ){ 
        SAMTOOLS_INDEX (BWAMEM2_MEM.out.bam)
        ch_sam = SAMTOOLS_INDEX.out.sam }

        if ( params.summary_txt_file == true ) {
        ch_summarytxt = summarytxt.map { file -> tuple(file.baseName, file) }

        PYCOQC (
            ch_summarytxt, MINIMAP2_ALIGN.out.bam, SAMTOOLS_INDEX.out.bai
        )
        ch_versions = ch_versions.mix(PYCOQC.out.versions)
        } else {
            ch_summarytxt = Channel.empty()
        }

        MERQURY (
            assemblies, ch_meryl, genome_size_est, params.tolerable_collision
        )
        ch_merqury
            .concat(MERQURY.out.assembly_qv)
            .set { ch_merqury }
        ch_versions = ch_versions.mix(MERQURY.out.versions)

    emit:
        ch_index
        ch_quast
        ch_busco
        ch_merqury
        ch_bam
        ch_busco_full_table

        
    versions = ch_versions                     // channel: [ versions.yml ]
}
