//include { GUNZIP } from '../../modules/nf-core/gunzip/main'
include { NANOPLOT } from '../../modules/nf-core/nanoplot/main'
include { TOTAL_BASES_LR } from '../../modules/local/total_bases_lr'
include { COVERAGE_LR } from '../../modules/local/coverage_lr'
include { FLYE } from '../../modules/nf-core/flye/main' 
include { MASURCA } from '../../modules/local/masurca'
include { CANU } from '../../modules/nf-core/canu/main' 
include { HIFIASM } from '../../modules/nf-core/hifiasm/main' 
include { VERKKO } from '../../modules/local/verkko'


workflow ASSEMBLY {

    take:
        longreads   // channel: [ val(meta), [ decontaminated and/or length filtered fastq ] ]
        shortreads
        genome_size_est
        genome_size_est_long
        combined_longreads
        ont_reads
        pacbio_reads
        ont_reads_w_meta

    main:
    ch_versions = Channel.empty() 

        if (params.ONT_lr_herrocorrected == false && params.PacBioHifi_lr == true && params.ONT_lr == true) {
        //makes sure long read input is filtered reads
        NANOPLOT (longreads)
        nanoplot_filtered_out   = NANOPLOT.out.html

        TOTAL_BASES_LR(NANOPLOT.out.txt)

        TOTAL_BASES_LR.out.total_bases
            .combine(genome_size_est_long)
            .set{ch_cov}

        COVERAGE_LR(ch_cov) } 
        else if (params.ONT_lr_herrocorrected == true && params.PacBioHifi_lr == true){
            NANOPLOT (pacbio_reads)
        nanoplot_filtered_out   = NANOPLOT.out.html

        TOTAL_BASES_LR(NANOPLOT.out.txt)

        TOTAL_BASES_LR.out.total_bases
            .combine(genome_size_est_long)
            .set{ch_cov}

        COVERAGE_LR(ch_cov)
        } else if (params.ONT_lr_herrocorrected == false && params.PacBioHifi_lr == false && params.ONT_lr == true){
        NANOPLOT (longreads)
        nanoplot_filtered_out   = NANOPLOT.out.html

        TOTAL_BASES_LR(NANOPLOT.out.txt)

        TOTAL_BASES_LR.out.total_bases
            .combine(genome_size_est_long)
            .set{ch_cov}

        COVERAGE_LR(ch_cov)
        } else{nanoplot_filtered_out = Channel.empty()}

        assemblies = Channel.empty() 
        
        //if statement to run assembly and create channels for each resulting assembly
        if ( params.flye == true ) {
            if (params.flye_mode == "both") {
                println "assembling ont and pb reads with flye!"
                longreads
                    .combine(genome_size_est)
                    .set{ch_flye_input}
            } else if (params.flye_mode == "pb"){
                println "assembling pb reads with flye!"
                pacbio_reads
                    .combine(genome_size_est)
                    .set{ch_flye_input}
            } else if (params.flye_mode == "ont"){
                println "assembling ont reads with flye!"
                ont_reads_w_meta
                    .combine(genome_size_est)
                    .set{ch_flye_input}
            } else {ch_flye_input = Channel.empty()}
            
            FLYE(ch_flye_input)
            flye_assembly      = FLYE.out.fasta   
            flye_assembly
                .map { file -> tuple(id: file.simpleName, file)  }
                .set { f_assembly }      
        } else {
            f_assembly = Channel.empty() 
            flye_assembly = Channel.empty() 
        }

        if ( params.canu == true ) {
            println "assembling long reads with canu!"
            CANU(longreads, genome_size_est)
            canu_assembly      = CANU.out.assembly   

            canu_assembly
                .map { file -> tuple(id: file.simpleName, file)  }
                .set { c_assembly }
        } else {
            canu_assembly = Channel.empty() 
            c_assembly = Channel.empty() 
        }

        if ( params.masurca == true && params.shortread == true) {
            if (params.ONT_lr == true && params.PacBioHifi_lr == true) {
            println "hybrid assembly with maSuRCA!"
            MASURCA(combined_longreads, shortreads)
            masurca_assembly    = MASURCA.out.fasta

            masurca_assembly
                .map { file -> tuple(id: file.simpleName, file)  }
                .set { m_assembly } 
            } else {
            println "hybrid assembly with maSuRCA!"
            MASURCA(longreads, shortreads)
            masurca_assembly    = MASURCA.out.fasta

            masurca_assembly
                .map { file -> tuple(id: file.simpleName, file)  }
                .set { m_assembly } 
        }} else {
            m_assembly = Channel.empty() 
            masurca_assembly = Channel.empty() 
        }

        if (params.hifiasm ==true){
            if (params.ONT_lr_ultralong == true && params.PacBioHifi_lr == true){
                HIFIASM(pacbio_reads, ont_reads)
                println "assembling long reads with hifiasm!"
                hifi_assembly    = HIFIASM.out.assembly_fasta

                hifi_assembly
                    .map { file -> tuple(id: file.simpleName, file)  }
                    .set { h_assembly }
            } else { 
                HIFIASM(pacbio_reads, [])
                println "assembling long reads with hifiasm!"
                hifi_assembly    = HIFIASM.out.assembly_fasta

                hifi_assembly
                    .map { file -> tuple(id: file.simpleName, file)  }
                    .set { h_assembly }           
        }} else {
            h_assembly = Channel.empty() 
            hifi_assembly = Channel.empty() 
        }

        if (params. verkko == true){
            println "assembling long reads with verkko!"
            if(params.ragtag_reference){
                if (params.ONT_lr == true && params.PacBioHifi_lr == true){
                    VERKKO(pacbio_reads, ont_reads, params.ragtag_reference)}
                else {
                    VERKKO(pacbio_reads, [], params.ragtag_reference)}
            } else {
                if (params.ONT_lr == true && params.PacBioHifi_lr == true){
                    VERKKO(pacbio_reads, ont_reads, [])}
                else {
                    VERKKO(pacbio_reads, [], [])}
            }
            
            verkko_assembly = VERKKO.out.fasta

            verkko_assembly
                .map { file -> tuple(id: file.simpleName, file)  }
                .set { v_assembly }           
        } else {
            verkko_assembly = Channel.empty()
            v_assembly = Channel.empty()
        }

        if ( params.ex_assembly == true ) {
            println "inputting existing assembly!"
            existing_assembly = Channel.fromPath(params.existing_assembly)

            existing_assembly
                .map { file -> tuple(id: file.simpleName, file)  }
                .set { ex_assembly }
        } else {
            ex_assembly = Channel.empty() 
            existing_assembly = Channel.empty() 
        }

        no_meta_assemblies = Channel.empty()

        no_meta_assemblies
            .concat(flye_assembly, canu_assembly, masurca_assembly, hifi_assembly, existing_assembly)
            .collect()
            .set { all_assemblies_no_meta }

        no_meta_assemblies
            .concat(flye_assembly, canu_assembly, masurca_assembly, hifi_assembly, existing_assembly)
            .flatten()
            .map { file -> tuple(id: file.simpleName, file) }
            .set { all_assemblies }

        no_meta_assemblies
            .concat(flye_assembly, canu_assembly, masurca_assembly, hifi_assembly, existing_assembly)
            .set{blob_test}

    emit:
        all_assemblies  
        longreads   
        nanoplot_filtered_out  
        f_assembly
        all_assemblies_no_meta
        blob_test
        
    versions = ch_versions                     // channel: [ versions.yml ]
}
