include { GUNZIP } from '../../modules/nf-core/gunzip/main'
include { CENTRIFUGE_CENTRIFUGE } from '../../modules/nf-core/centrifuge/centrifuge/main'
include { CENTRIFUGE_KREPORT } from '../../modules/nf-core/centrifuge/kreport/main'
include { NANOPLOT } from '../../modules/nf-core/nanoplot/main'
include { KMER_FREQ } from '../../modules/local/kmerfreq'
include { GCE } from '../../modules/local/gce'
include { RECENTRIFUGE_C } from '../../modules/local/recentrifuge/centrifuge'
include { SEQKIT_GREP } from '../../modules/local/seqkit/grep/main'
include { TAR } from '../../modules/local/tar'

workflow READ_QC {

    take:
  
        reads  // channel: [ val(meta), [ reads ] ]
           
    main:
    
    ch_versions = Channel.empty()

        if ( params.tar == true ) { 
            TAR(reads)
            reads = TAR.out.untar }

        if ( params.ONT_lr_herrocorrected == false ) {
            NANOPLOT(reads)
            nanoplot_reads_out   = NANOPLOT.out.html
            nanoplot_report_txt  = NANOPLOT.out.txt
        } else {
            nanoplot_reads_out = Channel.empty()
            nanoplot_report_txt = Channel.empty()
        }
        

        if ( params.centrifuge_db == null ){
            ch_db = Channel.empty() 
            centrifuge_out = Channel.empty() 
            reads
                .set{filtered_fastq}
            filtered_fastq
                .map { it[1] }
                .set { fastq_filt }
            }
        else if (params.centrifuge_db != null ){
            ch_db = Channel.fromPath(params.centrifuge_db)

            // if a centrifuge database is provided, run centrifuge and filter out all classified results
          
            if (params.centrifuge_ont){
                 CENTRIFUGE_CENTRIFUGE        ( reads, ch_db, params.save_unaligned, params.save_aligned, params.sam_format )
                 CENTRIFUGE_KREPORT           ( CENTRIFUGE_CENTRIFUGE.out.results, ch_db )
                 SEQKIT_GREP(CENTRIFUGE_CENTRIFUGE.out.results, reads)
                centrifuge_out       = CENTRIFUGE_KREPORT.out.kreport

                 fastq_filt           = SEQKIT_GREP.out.filter

                 fastq_filt
                    .map { file -> tuple([id:file.simpleName, single_end:true], file)  }
                    .set { filtered_fastq }

                 if( params.rcf_db ){
                    RECENTRIFUGE_C(CENTRIFUGE_CENTRIFUGE.out.results, params.rcf_db)
                 }
            } else {
                 reads
                    .set{filtered_fastq}
                 filtered_fastq
                    .map { it[1] }
                    .set { fastq_filt }
                centrifuge_out = Channel.empty()
            }
        }

        KMER_FREQ(filtered_fastq)

        GCE(KMER_FREQ.out.kmerstat, KMER_FREQ.out.kmernum)
        gce_genome_size      = GCE.out.gce2log

    emit:
        filtered_fastq    // channel: [ val(meta), path(decontaminated fastq) ]
        nanoplot_reads_out   
        centrifuge_out 
        gce_genome_size
        nanoplot_report_txt 
        fastq_filt
        
    versions = ch_versions                     // channel: [ versions.yml ]
}
