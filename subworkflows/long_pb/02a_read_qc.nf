include { CUTADAPT } from '../../modules/local/cutadapt'
include { GUNZIP } from '../../modules/nf-core/gunzip/main'
include { GENOMESCOPE2 } from '../../modules/nf-core/genomescope2/main'
include { JELLYFISH_KMER } from '../../modules/local/jellyfish_kmer'
include { JELLYFISH_HIST } from '../../modules/local/jellyfish_hist'
include { KRAKEN2_KRAKEN2_PB } from '../../modules/nf-core/kraken2/kraken2/main'
include { RECENTRIFUGE_KR } from '../../modules/local/recentrifuge/kraken'
include { NANOPLOT } from '../../modules/nf-core/nanoplot/main'
include { TOTAL_BASES_LR } from '../../modules/local/total_bases_lr' 
include { RUN_QC } from '../../modules/local/runqc'

workflow READ_QC3 {

    take:
  
        input_pacbio  // channel: [ val(meta), [ reads ] ]
           
    main:
    
    ch_versions = Channel.empty()

        NANOPLOT(input_pacbio)
        TOTAL_BASES_LR (NANOPLOT.out.txt)

        if (params.pb_xml) {
            RUN_QC(params.pb_xml)
        }

	    CUTADAPT (input_pacbio)
	    ch_versions = ch_versions.mix(CUTADAPT.out.versions)

        CUTADAPT.out.reads
                .map { file -> tuple([id:file.simpleName, single_end:true], file)  }
                .set { adaptor_trimmed }

        if ( params.kraken_db == null ){
            ch_db = Channel.empty()
        } else if (params.kraken_db != null ){
            ch_db = Channel.fromPath(params.kraken_db)
        }
    
        if (params.kraken_pb == true){
            
            //decontamination of trimmed short reads
            KRAKEN2_KRAKEN2_PB(adaptor_trimmed, ch_db, params.save_output_fastqs, params.save_reads_assignment)
            filt_pbhifi = KRAKEN2_KRAKEN2_PB.out.unclassified_reads_fastq   

            filt_pbhifi
                .map { file -> tuple([id:file.simpleName, single_end:true], file)  }
                .set { filtered_fastq }

            if( params.rcf_db ){
                //summarizing and visualizing decontam
                RECENTRIFUGE_KR(KRAKEN2_KRAKEN2_PB.out.classified_reads_assignment, params.rcf_db)
            }

        } else {
            filt_pbhifi = CUTADAPT.out.reads
            adaptor_trimmed
                .set{filtered_fastq}
        }

        
        GUNZIP(filtered_fastq)

        JELLYFISH_KMER(GUNZIP.out.gunzip, params.kmer_num)
        JELLYFISH_HIST(JELLYFISH_KMER.out.shortkmer, params.kmer_num)
        
        GENOMESCOPE2(JELLYFISH_HIST.out.shortkmer_hist, params.kmer_num)

    emit:
        filtered_fastq    // channel: [ val(meta), path(decontaminated fastq) ]
        nanoplot_reads_out   = NANOPLOT.out.html
        nanoplot_report_txt  = NANOPLOT.out.txt
        base_count           = TOTAL_BASES_LR.out.total_bases
        genome_size_est = GENOMESCOPE2.out.summary
        filt_pbhifi
        GUNZIP.out.gunzip
        
    versions = ch_versions                     // channel: [ versions.yml ]
}
