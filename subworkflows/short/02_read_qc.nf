include { GZIP } from '../../modules/local/gzip'
include { FASTP } from '../../modules/nf-core/fastp/main'
include { FASTQC } from '../../modules/nf-core/fastqc/main'
include { KRAKEN2_KRAKEN2_SR } from '../../modules/local/kraken2_sr'
include { FASTQC_2 } from '../../modules/local/fastqc2/main'
include { FASTQC_3 } from '../../modules/local/fastqc3/main'
include { GUNZIP } from '../../modules/nf-core/gunzip/main'
include { GENOMESCOPE2 } from '../../modules/nf-core/genomescope2/main'
include { JELLYFISH_KMER } from '../../modules/local/jellyfish_kmer'
include { JELLYFISH_HIST } from '../../modules/local/jellyfish_hist'
include { RECENTRIFUGE_KR } from '../../modules/local/recentrifuge/kraken'

workflow READ_QC2 {

    take:
        shortreads  // channel: [ val(meta), [ reads ] ] 
           
    main:
    
    ch_versions = Channel.empty()

        //zipping files for downstream processing
        //GZIP(shortreads)

        //qc raw short reads
        FASTQC(shortreads)

        //adapter trimming
        FASTP(shortreads, params.adapter_fasta, params.save_trimmed_fail, params.save_merged, params.dedup, params.quality_trimming, params.length_cutoff)

        //qc adapter trimmed short reads
        FASTQC_2(FASTP.out.reads)

        if ( params.kraken_db == null ){
            ch_db = Channel.empty()
        } else if (params.kraken_db != null ){
            ch_db = Channel.fromPath(params.kraken_db)
        }

        if (params.kraken_ill == true){
        //decontamination of trimmed short reads
        KRAKEN2_KRAKEN2_SR(FASTP.out.reads, ch_db, params.save_output_fastqs, params.save_reads_assignment)
        
            if( params.rcf_db ){
            //summarizing and visualizing decontam
            RECENTRIFUGE_KR(KRAKEN2_KRAKEN2_SR.out.classified_reads_assignment, params.rcf_db)
            }

            filt_shortreads = KRAKEN2_KRAKEN2_SR.out.unclassified_reads_fastq   

            //qc decontaminated short reads
            FASTQC_3(KRAKEN2_KRAKEN2_SR.out.unclassified_reads_fastq)
    
        } else {
            filt_shortreads = FASTP.out.reads
        }

        filt_shortreads
            .map { it[1] }
            .set { filt_sr_no_meta }

        //unzip decontam short reads
        GUNZIP(filt_shortreads)

        GUNZIP.out.gunzip
            .map { it[1] }
            .set { unzip_filt_sr_no_meta }

        JELLYFISH_KMER(GUNZIP.out.gunzip, params.kmer_num)
        JELLYFISH_HIST(JELLYFISH_KMER.out.shortkmer, params.kmer_num)
        
        GENOMESCOPE2(JELLYFISH_HIST.out.shortkmer_hist, params.kmer_num)

    emit:
        filt_shortreads
        GUNZIP.out.gunzip
        fastp_report = FASTP.out.json
        genome_size_est = GENOMESCOPE2.out.summary
        filt_sr_no_meta
        unzip_filt_sr_no_meta

    versions = ch_versions                     // channel: [ versions.yml ]
}
