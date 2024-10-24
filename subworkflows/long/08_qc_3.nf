include { QUAST } from '../../modules/local/quast'  
include { BUSCO } from '../../modules/nf-core/busco/main'   
include { MERYL_COUNT } from '../../modules/nf-core/meryl/count/main' 
include { MERQURY } from '../../modules/nf-core/merqury/main' 
include { SAMTOOLS_INDEX } from '../../modules/nf-core/samtools/index/main' 
include { BWAMEM2_INDEX } from '../../modules/nf-core/bwamem2/index/main' 
include { BWAMEM2_MEM } from '../../modules/nf-core/bwamem2/mem/main' 
include { COMPLEASM } from '../../modules/local/compleasm'  
include { WINNOWMAP } from '../../modules/local/winnowmap'  
include { SAMTOOLS_SORT } from '../../modules/nf-core/samtools/sort'
include { MINIMAP2_INDEX } from '../../modules/nf-core/minimap2/index/main' 
include { MINIMAP2_ALIGN } from '../../modules/nf-core/minimap2/align/main' 

workflow QC_3 {

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
        no_meta_fq
        meryl_repk
        merqury_comp

    main:

    ch_versions = Channel.empty() 

    if ( params.shortread == true ) {
        BWAMEM2_INDEX(assemblies)

        shortreads
            .combine(BWAMEM2_INDEX.out.index)
            .set{bwa}

        BWAMEM2_MEM(bwa, params.samtools_sort)
        if(params.longread == false){ch_bam = BWAMEM2_MEM.out.bam}
    }

    if ( params.longread == true ){
        assemblies
            .combine(no_meta_fq)
            .set{align_ch}

        align_ch
            .combine(meryl_repk)
            .set{winnowmap_ch}

        if(params.winnowmap == true){
        WINNOWMAP(winnowmap_ch, params.kmer_num)
        ch_sam = WINNOWMAP.out.sam

        SAMTOOLS_SORT(ch_sam)
        ch_bam = SAMTOOLS_SORT.out.bam

        ch_index = Channel.empty() }

        if(params.minimap2 == true){
        MINIMAP2_INDEX(assemblies)
        ch_versions = ch_versions.mix(MINIMAP2_INDEX.out.versions)
        ch_index = MINIMAP2_INDEX.out.index

        // align reads
        MINIMAP2_ALIGN(align_ch, params.bam_format, params.cigar_paf_format, params.cigar_bam)
        ch_bam = MINIMAP2_ALIGN.out.bam } 
    }
    else {
        ch_index = Channel.empty()  
    }
        
        // run quast
        QUAST(
            assemblies
        )
        ch_quast
            .concat(QUAST.out.results)
            .set { ch_quast }
        ch_versions = ch_versions.mix(QUAST.out.versions)

        // run BUSCO or compleasm
        if (params.busco == true){
            BUSCO(assemblies, params.busco_lineage, [], [])
            ch_busco
                .concat(BUSCO.out.short_summaries_txt)
                .set { ch_busco } 
            ch_busco_full_table = BUSCO.out.full_table
            ch_versions = ch_versions.mix(BUSCO.out.versions)
        }

        if (params.compleasm == true){
            COMPLEASM(assemblies, params.compleasm_lineage)
            ch_busco 
                .concat(COMPLEASM.out.txt)
                .set { ch_busco } 
            if (params.busco == false){
                ch_busco_full_table = Channel.empty() 
            }
        }

    if ( params.longread == true ){
        SAMTOOLS_INDEX (ch_bam)
    } else if ( params.shortread == true ){ 
        SAMTOOLS_INDEX (BWAMEM2_MEM.out.bam) }

        assemblies
            .combine(ch_meryl)
            .set{ch_input_merqury}

        ch_input_merqury
            .combine(genome_size_est)
            .set{merqury_paths}

        MERQURY (
            merqury_paths, params.tolerable_collision
        )
        ch_merqury
            .concat(MERQURY.out.assembly_qv)
            .set { ch_merqury }

        merqury_comp
            .concat(MERQURY.out.stats)
            .set{ch_merqury_comp}
        ch_versions = ch_versions.mix(MERQURY.out.versions)

    emit:
        ch_index = SAMTOOLS_INDEX.out.bai
        ch_quast
        ch_busco
        ch_merqury
        ch_bam
        ch_busco_full_table
        ch_merqury_comp

        
    versions = ch_versions                     // channel: [ versions.yml ]
}
