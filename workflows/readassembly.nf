/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowGenomeassembly.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.centrifuge_db ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) }
if (params.shortinput) { ch_shortinput = file(params.shortinput) }
if (params.pb_input) { ch_pb_input = file(params.pb_input)}
if (params.centrifuge_db) { ch_db = file(params.centrifuge_db) }
//if (params.summary_txt) {ch_sequencing_summary = file(params.sequencing_summary) } else { ch_sequencing_summary = []}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
//ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
//ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
//ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// MODULES
include { CAT } from '../modules/local/cat' 
include { OUTPUT } from '../modules/local/output' 
include { TOTAL_BASES_SR } from '../modules/local/total_bases_sr' 
include { TOTAL_BASES_LR } from '../modules/local/total_bases_lr' 
include { COVERAGE_SR } from '../modules/local/coverage_sr'
include { COVERAGE_LR } from '../modules/local/coverage_lr'
include { COVERAGE_LR_PB } from '../modules/local/coverage_lr_pb'
include { MASURCA_SR_ADV } from '../modules/local/masurca_sr_adv'
include { MASURCA_SR } from '../modules/local/masurca_sr'
include { REDUNDANS_A } from '../modules/local/redundans_assembler'
include { FORMAT } from '../modules/local/format_genome_size'
include { EXTRACT_LR } from '../modules/local/extract_genome_size'
include { EXTRACT_SR } from '../modules/local/extract_short_genome_size'
include { EXTRACT_PB } from '../modules/local/extract_pb_genome_size'

// SUBWORKFLOWS
include { INPUT_CHECK } from '../subworkflows/long/01_input_check'
include { READ_QC } from '../subworkflows/long/02a_read_qc'
include { LENGTH_FILT } from '../subworkflows/long/02b_length_filter'
include { ASSEMBLY } from '../subworkflows/long/03_assembly'
include { QC_1 } from '../subworkflows/long/04_qc_1'
include { POLISH } from '../subworkflows/long/05_polish'
include { QC_2 } from '../subworkflows/long/06_qc_2'
include { HAPS } from '../subworkflows/long/07_purge'
include { QC_3 } from '../subworkflows/long/08_qc_3'
include { SCAFFOLD } from '../subworkflows/long/12_scaffold'
include { QC_4 } from '../subworkflows/long/13_qc_4'
include { VISUALIZE } from '../subworkflows/long/14_visualization'


include { INPUT_CHECK2 } from '../subworkflows/short/01_input_check'
include { READ_QC2 } from '../subworkflows/short/02_read_qc'
include { POLISH2 } from '../subworkflows/short/03_polish'
include { PURGE2 } from '../subworkflows/short/04_purge'

include { INPUT_CHECK3 } from '../subworkflows/long_pb/01_input_check'
include { READ_QC3 } from '../subworkflows/long_pb/02a_read_qc'
include { LENGTH_FILT3 } from '../subworkflows/long_pb/02b_length_filter'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GENOMEASSEMBLY {

    ch_versions = Channel.empty()

    if ( params.shortread == true || params.PacBioHifi_lr == true) {
        ch_kraken_db = Channel.fromPath(params.kraken_db)
    }

    if (params.longread == true){
        if (params.ONT_lr == true) {
            ch_data = INPUT_CHECK ( ch_input )
            ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    
            ch_centrifuge_db = Channel.fromPath(params.centrifuge_db)
            //decontamination and quality checking of long reads
            READ_QC (ch_data.reads, ch_centrifuge_db)
            ch_versions = ch_versions.mix(READ_QC.out.versions)
            no_meta_fastq = READ_QC.out[5]

            //optional read filtering by length with bioawk
            LENGTH_FILT (READ_QC.out[0], no_meta_fastq)
            ch_ONTlongreads = LENGTH_FILT.out[0]
            LENGTH_FILT.out[1]
                .set{no_meta_ch_ONT}

            ch_versions = ch_versions.mix(LENGTH_FILT.out.versions)

            if (params.PacBioHifi_lr == false){
                ch_longreads = LENGTH_FILT.out[0]
            }
        } else {
            ch_ONTlongreads = Channel.empty()
            no_meta_fastq = Channel.empty()
            no_meta_ch_ONT = Channel.empty()}
        if (params.PacBioHifi_lr == true) {
            ch_PB_data = INPUT_CHECK3 ( ch_pb_input )
            ch_versions = ch_versions.mix(INPUT_CHECK3.out.versions)    
        
            //decontamination and quality checking of long reads
            READ_QC3 (ch_PB_data.reads, ch_kraken_db)
            no_meta_decontamPB = READ_QC3.out[5]

            ch_versions = ch_versions.mix(READ_QC3.out.versions)

            LENGTH_FILT3 (READ_QC3.out[6], no_meta_decontamPB)
            ch_PacBiolongreads = LENGTH_FILT3.out[0]
            no_meta_ch_PB = LENGTH_FILT3.out[1]

            ch_versions = ch_versions.mix(LENGTH_FILT3.out.versions)

            if (params.ONT_lr == false){
                ch_longreads = LENGTH_FILT3.out[0]
        }} else {
            ch_PacBiolongreads = Channel.empty()
            no_meta_ch_PB = Channel.empty()}
        if (params.PacBioHifi_lr == true || params.ONT_lr == true || params.longread == true) {
            ch_ONTlongreads
                .concat(ch_PacBiolongreads)
                .set{ch_longreads}

            ch_ONTlongreads
                .combine(no_meta_ch_PB)
                .set{ch_combo_longreads}

            no_meta_ch_ONT
                .concat(no_meta_ch_PB)
                .flatten()
                .map{ file -> tuple(id: file.baseName, file) }
                .set{ch_flat_lr}

        }else if (params.PacBioHifi_lr == false && params.ONT_lr == false){ 
            ch_longreads = Channel.empty()
            ch_combo_longreads = Channel.empty()
        }}

    if ( params.shortread == true ) {
        ch_shortdata = INPUT_CHECK2 ( ch_shortinput )
    ch_versions = ch_versions.mix(INPUT_CHECK2.out.versions)

        //adaptor trimming and decontamination of short reads if available
        READ_QC2 (ch_shortdata.reads, ch_kraken_db)
    ch_versions = ch_versions.mix(READ_QC2.out.versions)
        filt_sr_unzip = READ_QC2.out[1]
        filt_sr_nometa = READ_QC2.out[4]
    } else {
        filt_sr_unzip = Channel.empty()
        filt_sr_nometa = Channel.empty()
    }

    // extracting and formatting genome size est

    if(params.shortread == true){
        ill_genome_size = READ_QC2.out[3]
        EXTRACT_SR(ill_genome_size)
        ill_readable_size = EXTRACT_SR.out[0]
        ill_full_size = EXTRACT_SR.out[1]
    }

    if (params.longread == true && params.ONT_lr == true){
        ont_genome_size = READ_QC.out[3]
        EXTRACT_LR(ont_genome_size)
        ont_readable_size = EXTRACT_LR.out[0]
        ont_full_size = EXTRACT_LR.out[1] 
    }

    if (params.longread == true && params.PacBioHifi_lr == true){
        pb_genome_size = READ_QC3.out[4]
        EXTRACT_PB(pb_genome_size)
        pb_readable_size = EXTRACT_PB.out[0]
        pb_full_size = EXTRACT_PB.out[1]
    }

    
    if (params.manual_genome_size){
        genome_size = params.manual_genome_size
        FORMAT(genome_size)
        readable_size = FORMAT.out[0]
        full_size = FORMAT.out[1]
    } else if (params.shortread == true){
        readable_size = EXTRACT_SR.out[0]
        full_size = EXTRACT_SR.out[1]
    } else if (params.longread == true ){
        if (params.PacBioHifi_lr == true) {
            readable_size = EXTRACT_PB.out[0]
            full_size = EXTRACT_PB.out[1]
            }    
        else if (params.ONT_lr == true) {
            readable_size = EXTRACT_LR.out[0]
            full_size = EXTRACT_LR.out[1]
            
        }
    }

    //calculating coverage for long and/or short reads
    if (params.longread == true && params.ONT_lr == true){
        TOTAL_BASES_LR (READ_QC.out[4])
        TOTAL_BASES_LR.out.total_bases
            .combine(full_size)
            .set{ch_ont_cov}
        COVERAGE_LR (ch_ont_cov)
    }   
    if (params.longread == true && params.PacBioHifi_lr == true){
        COVERAGE_LR_PB (full_size, READ_QC3.out[3])
    }
    if (params.shortread == true) {
        TOTAL_BASES_SR (READ_QC2.out[2])
        COVERAGE_SR (full_size, TOTAL_BASES_SR.out.total_bases_before, TOTAL_BASES_SR.out.total_bases_after)
    }

    if (params.ONT_lr == true && params.PacBioHifi_lr == true) {
        CAT(ch_combo_longreads)
        CAT.out.cat_longreads
                .map { file -> tuple(id: file.baseName, file)  }
                .set { combined_lr }
    } else {
        combined_lr = Channel.empty()
    }
    //long read and hybrid assemblies
    if (params.longread == true && params.shortread == true){
        //assembly inputting long & short reads
        ASSEMBLY (ch_longreads, READ_QC2.out[1], readable_size, full_size, combined_lr, no_meta_fastq, ch_PacBiolongreads)
        lr_assemblies   = ASSEMBLY.out[4]
    } else if (params.longread == true && params.shortread == false) {
        ch_shortdata = Channel.empty() 
        //assembly of decontam and length filtered (if specified) long reads
        ASSEMBLY (ch_longreads, [], readable_size, full_size, combined_lr, no_meta_fastq, ch_PacBiolongreads)
        lr_assemblies   = ASSEMBLY.out[4]
    ch_versions = ch_versions.mix(ASSEMBLY.out.versions)   
    } else {
        lr_assemblies = Channel.empty()
    }

    //short read only assembly
    if ( params.shortread == true && params.masurca == true && params.longread == false){
        if (params.masurca_sr_adv == true){
            ch_config = Channel.fromPath(params.masurca_config)
            MASURCA_SR_ADV (ch_config)
            println "assembling short reads with maSuRCA!"
            MASURCA_SR_ADV.out.fasta
                .set{masurca_asm}
            masurca_asm
                .map { file -> tuple(id: file.baseName, file)  }
                .set { masurca_sr_assembly }
        } else {
            MASURCA_SR (READ_QC2.out[1])
            println "assembling short reads with maSuRCA!"
            MASURCA_SR.out.fasta
                .set{masurca_asm}
            masurca_asm
                .map { file -> tuple(id: file.baseName, file)  }
                .set { masurca_sr_assembly }
        }
        
        
    } else {
        masurca_asm = Channel.empty() 
        masurca_sr_assembly = Channel.empty() 
    }

    if ( params.shortread == true && params.redundans == true){
        REDUNDANS_A (ch_shortdata.reads)
        println "assembling short reads with redundans!"
        REDUNDANS_A.out.assembly_fasta
            .set{redundans_asm}
        redundans_asm
            .map { file -> tuple(id: file.baseName, file)  }
            .set { redundans_assembly }
        
    } else {
        redundans_asm = Channel.empty()
        redundans_assembly = Channel.empty() 
    }
    
    if ( params.shortread == true) {
    masurca_asm
        .concat(redundans_asm)
        .collect()
        .set { sr_assemblies }

    } else {
        sr_assemblies = Channel.empty()
    }
    
    lr_assemblies
        .concat(sr_assemblies)
        .flatten()
        .map { file -> tuple(id: file.baseName, file) }
        .set{all_assemblies}

    lr_assemblies
        .concat(sr_assemblies)
        .set{all_assemblies_nm}

    if ( params.summary_txt_file == true) {
        ch_summtxt = Channel.fromPath(params.summary_txt) 
    } else {
        ch_summtxt = Channel.empty() }

    if ( params.shortread == true && params.longread == true ) {
        QC_1 (all_assemblies, ch_longreads, ch_summtxt, READ_QC2.out[0], full_size, ch_flat_lr)
    ch_versions = ch_versions.mix(QC_1.out.versions)
    } else if ( params.longread == true && params.shortread == false ) {
        QC_1 (all_assemblies, ch_longreads, ch_summtxt, [], full_size, ch_flat_lr)
    ch_versions = ch_versions.mix(QC_1.out.versions)
    } else if ( params.shortread == true && params.longread == false ) {
        QC_1 (all_assemblies, READ_QC2.out[0], ch_summtxt, READ_QC2.out[0], full_size, READ_QC2.out[0])
    ch_versions = ch_versions.mix(QC_1.out.versions)}
    
    busco_tsv = QC_1.out[9]
    bam_1 = QC_1.out[1]

    QC_1.out[10].view()

    ASSEMBLY.out[0]
        .join(QC_1.out[8])
        .set{assembly_sam_combo}

    if (params.racon_polish == true && params.ONT_lr == true && params.PacBioHifi_lr == true){
        assembly_sam_combo
            .combine(CAT.out.cat_longreads)
            .set { ch_racon }

    } else if(params.racon_polish == true && params.PacBioHifi_lr == true){
        assembly_sam_combo
            .combine(no_meta_ch_PB)
            .set { ch_racon }

    } else if(params.racon_polish == true && params.ONT_lr == true){
        assembly_sam_combo
            .combine(no_meta_ch_ONT)
            .set { ch_racon }

    } else { ch_racon = Channel.empty() }

    //polish assemblies
     if ( params.longread == true) {
        if ( params.medaka_polish == true || params.racon_polish == true){
            if (params.PacBioHifi_lr == true && params.ONT_lr == true){
                POLISH (ASSEMBLY.out[0], ch_PacBiolongreads, params.model, QC_1.out[8], assembly_sam_combo, no_meta_ch_ONT, no_meta_ch_PB)
            } else if (params.PacBioHifi_lr == false && params.ONT_lr == true){
                POLISH (ASSEMBLY.out[0], ch_ONTlongreads, params.model, QC_1.out[8], assembly_sam_combo, no_meta_ch_ONT, [])
            } else if (params.PacBioHifi_lr == true && params.ONT_lr == false){
                POLISH (ASSEMBLY.out[0], ch_PacBiolongreads, params.model, QC_1.out[8], assembly_sam_combo, [], no_meta_ch_PB)}

        POLISH.out[0] 
            .set{medaka_racon_polish}
        
    ch_versions = ch_versions.mix(POLISH.out.versions) } else {medaka_racon_polish = Channel.empty()}
    } else {
        medaka_racon_polish = Channel.empty()
    }
    //align assemblies to short reads and polish with POLCA if short reads are available
    if ( params.longread == true && params.shortread == true) {
        if (params.medaka_polish == true || params.racon_polish == true){
            POLISH2 (POLISH.out[0], READ_QC2.out[1])
        } else {
            POLISH2 (ASSEMBLY.out[0], READ_QC2.out[1])
        }
        POLISH2.out[0]
            .set{sr_polish}

        sr_polish
                .map { file -> tuple(id: file.baseName, file)  }
                .set { polca_polish }   
        
        //combine polished flye assemblies w other assemblies
        sr_polish
            .concat(medaka_racon_polish)
            .flatten()
            .map { file -> tuple(id: file.baseName, file) }
            .set { polished_assemblies }

    ch_versions = ch_versions.mix(POLISH2.out.versions)
    } else {
        sr_polish   = Channel.empty()
        medaka_racon_polish
            .flatten()
            .map { file -> tuple(id: file.baseName, file) }
            .set { polished_assemblies }


    }

    all_assemblies_nm
        .concat(medaka_racon_polish, sr_polish)
        .flatten()
        .map { file -> tuple(id: file.baseName, file) }
        .set { polished_assemblies_and_no_polish }


    if ( params.medaka_polish == true || params.racon_polish == true || params.shortread == true) {
        if ( params.shortread == true && params.longread == true ) {
            QC_2 (polished_assemblies, ch_longreads, ch_summtxt, QC_1.out[3], QC_1.out[4], QC_1.out[5], READ_QC2.out[0], QC_1.out[2], full_size, QC_1.out[7])
            ch_versions = ch_versions.mix(QC_2.out.versions)
        } else if ( params.longread == true && params.shortread == false ) {
            QC_2 (polished_assemblies, ch_longreads, ch_summtxt, QC_1.out[3], QC_1.out[4], QC_1.out[5], [], QC_1.out[2], full_size, QC_1.out[7])
            ch_versions = ch_versions.mix(QC_2.out.versions)
        } else if ( params.shortread == true && params.longread == false ) {
            QC_2 (polished_assemblies, READ_QC2.out[0], ch_summtxt, QC_1.out[3], QC_1.out[4], QC_1.out[5], READ_QC2.out[0], QC_1.out[2], full_size, QC_1.out[7])
    } 
    busco_tsv
        .concat(QC_2.out[6]) 
        .set{busco_tsv}

    bam_2 = QC_2.out[1]
    qc_quast = QC_2.out[3]
    qc_busco = QC_2.out[4]
    qc_merqury = QC_2.out[5]
    } else {
    bam_2 = Channel.empty()
    qc_quast = QC_1.out[3]
    qc_busco = QC_1.out[4]
    qc_merqury = QC_1.out[5]
    }
        

    purged_assemblies_common = Channel.empty()

    if (params.longread == true && params.purge == true) {
        HAPS (polished_assemblies_and_no_polish, ch_longreads)
        lr_purge = HAPS.out[0]
        lr_purge
            .concat(purged_assemblies_common)
            .set { purged_assemblies_common }
        HAPS.out[1]
            .set{no_meta_lr_purge}
        ch_versions = ch_versions.mix(HAPS.out.versions)
    } else {
        lr_purge = Channel.empty()
        no_meta_lr_purge = Channel.empty()
    }

    if (params.shortread == true && params.purge == true) {
        sr_assemblies
            .flatten()
            .map { file -> tuple(file.baseName, file) }
            .set{assemblies_sr_meta}
        PURGE2 (assemblies_sr_meta, READ_QC2.out[1])
        sr_purge = PURGE2.out[0]
        sr_purge
            .concat(purged_assemblies_common)
            .set{purged_assemblies_common}
        PURGE2.out[1]
            .set{no_meta_sr_purge}
    } else {
        sr_purge = Channel.empty()
        no_meta_sr_purge = Channel.empty()
    }


    if ( params.purge == true ) {
    if ( params.shortread == true && params.longread == true) {
        QC_3 (purged_assemblies_common, ch_longreads, ch_summtxt, qc_quast, qc_busco, qc_merqury, READ_QC2.out[0], full_size, QC_1.out[7])
    ch_versions = ch_versions.mix(QC_3.out.versions)
    } else if ( params.longread == true && params.shortread == false) {
        QC_3 (purged_assemblies_common, ch_longreads, ch_summtxt, qc_quast, qc_busco, qc_merqury, [], full_size, QC_1.out[7])  
    } else if ( params.shortread == true && params.longread == false) {
        QC_3 (purged_assemblies_common, READ_QC2.out[0], ch_summtxt, qc_quast, qc_busco, qc_merqury, READ_QC2.out[0], full_size, QC_1.out[7])
    } 
    busco_tsv
        .concat(QC_3.out[5]) 
        .set{busco_tsv}
    bam_3 = QC_3.out[4] 
    } else {
    bam_3 = Channel.empty()
    }
            
    if (params.ragtag_scaffold == true) {
        if (params.purge == true){
        qc_quast = QC_3.out[1]
        qc_busco = QC_3.out[2]
        qc_merqury = QC_3.out[3]} 

        ch_reference = Channel.fromPath(params.ragtag_reference)
        SCAFFOLD (purged_assemblies_common, ch_reference)
    ch_versions = ch_versions.mix(SCAFFOLD.out.versions)

        final_assemblies = SCAFFOLD.out[0]
        if ( params.shortread == true && params.longread == true ) {
            QC_4 (SCAFFOLD.out[0], ch_longreads, ch_summtxt, qc_quast, qc_busco, qc_merqury, READ_QC2.out[0], full_size, QC_1.out[7]) 
        } else if ( params.longread == true && params.shortread == false ) {
            QC_4 (SCAFFOLD.out[0], ch_longreads, ch_summtxt, qc_quast, qc_busco, qc_merqury, [], full_size, QC_1.out[7])  
        } else if ( params.shortread == true && params.longread == false ) {
            QC_4 (SCAFFOLD.out[0], READ_QC2.out[0], ch_summtxt, qc_quast, qc_busco, qc_merqury, READ_QC2.out[0], full_size, QC_1.out[7]) }

        busco_tsv
            .concat(QC_4.out[5]) 
            .set{busco_tsv}

        bam_4 = QC_4.out[4]

        SCAFFOLD.out[0]
            .concat(purged_assemblies_common, polished_assemblies, all_assemblies)
            .collect()
            .set{final_assemblies}

        SCAFFOLD.out[1]
            .concat(no_meta_lr_purge, no_meta_sr_purge, medaka_racon_polish, sr_polish, ASSEMBLY.out[4], masurca_asm, redundans_asm)
            .flatten()
            .map { file -> tuple(id: file.baseName, file)  }
            .set{ch_all_assemblies}

    } else {
        bam_4 = Channel.empty()

        ASSEMBLY.out[5]
            .concat(no_meta_lr_purge, no_meta_sr_purge, medaka_racon_polish, sr_polish, masurca_asm, redundans_asm)
            .flatten()
            .map { file -> tuple(id: file.baseName, file)  }
            .set{ch_all_assemblies}

        all_assemblies
            .concat(polished_assemblies, purged_assemblies_common)
            .collect() 
            .set {final_assemblies}
    }

    if ( params.ragtag_scaffold == true ) {
        ch_quast = QC_4.out[1]
        ch_busco = QC_4.out[2]
        ch_merqury = QC_4.out[3]
    } else if (params.purge == true){
        ch_quast = QC_3.out[1]
        ch_busco = QC_3.out[2]
        ch_merqury = QC_3.out[3]
    } else if (params.medaka_polish || params.racon_polish == true || params.shortread == true) {
        ch_quast = QC_2.out[3]
        ch_busco = QC_2.out[4]
        ch_merqury = QC_2.out[5]
    } else {
        ch_quast = QC_1.out[3]
        ch_busco = QC_1.out[4]
        ch_merqury = QC_1.out[5]
    }

    bam_1
        .concat(bam_2, bam_3, bam_4)
        .set{qc_bam}

    if (params.blobtools_visualization == true){
        ch_all_assemblies.view { "Final Assemblies: $it" }

        if(params.shortread == true && params.ONT_lr == true && params.PacBioHifi_lr == true){
            VISUALIZE(ch_all_assemblies, no_meta_ch_ONT, no_meta_ch_PB, filt_sr_nometa, qc_bam, busco_tsv)
        } else if (params.shortread == false && params.ONT_lr == true && params.PacBioHifi_lr == true){
            VISUALIZE(ch_all_assemblies, no_meta_ch_ONT, no_meta_ch_PB, [], qc_bam, busco_tsv)
        } else if (params.shortread == true && params.ONT_lr == false && params.PacBioHifi_lr == true){
            VISUALIZE(ch_all_assemblies, [], no_meta_ch_PB, filt_sr_nometa, qc_bam, busco_tsv)
        } else if (params.shortread == true && params.ONT_lr == true && params.PacBioHifi_lr == false){
            VISUALIZE(ch_all_assemblies, no_meta_ch_ONT, [], filt_sr_nometa, qc_bam, busco_tsv)
        } else if (params.shortread == true && params.ONT_lr == false && params.PacBioHifi_lr == false){
            VISUALIZE(ch_all_assemblies, [], [], filt_sr_nometa, qc_bam, busco_tsv)
        } else if (params.shortread == false && params.ONT_lr == true && params.PacBioHifi_lr == false){
            VISUALIZE(ch_all_assemblies, no_meta_ch_ONT, [], [], qc_bam, busco_tsv)
        } else if (params.shortread == false && params.ONT_lr == false && params.PacBioHifi_lr == true){
            VISUALIZE(ch_all_assemblies, [], no_meta_ch_PB, [], qc_bam, busco_tsv)
        }
    }

    ch_quast
        .join(ch_busco, ch_merqury)
        .set{ch_output}
    OUTPUT (ch_output)

    assembly_stats  =   OUTPUT.out.assemblyStats

    //
    // MODULE: MultiQC
    //
   //workflow_summary    = WorkflowGenomeassembly.paramsSummaryMultiqc(workflow, summary_params)
   //ch_workflow_summary = Channel.value(workflow_summary)

    //methods_description    = WorkflowGenomeassembly.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    //ch_methods_description = Channel.value(methods_description)

    //ch_multiqc_files = Channel.empty()
    //ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    //ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    //ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    //ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    //MULTIQC (
    //    ch_multiqc_files.collect(),
    //    ch_multiqc_config.toList(),
    //    ch_multiqc_custom_config.toList(),
    //    ch_multiqc_logo.toList()
   //)
    //multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
       NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
   }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
