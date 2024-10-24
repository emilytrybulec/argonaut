include { MEDAKA } from '../../modules/nf-core/medaka/main'  
include { RACON } from '../../modules/nf-core/racon/main'  

workflow POLISH {

    take:
        assembly
        fastq_filt
        model
        paf
        ch_racon
        no_meta_ont
        no_meta_pb
    main:

    ch_versions = Channel.empty() 

        if (params.racon_polish == true) {
            println "polishing assemblies with racon!"

            if (params.ONT_lr  == true && params.PacBioHifi_lr == true){
                ch_racon
                    .combine(no_meta_pb)
                    .set{ch_racon_full}
            } else if (params.ONT_lr  == true && params.PacBioHifi_lr == false){
                ch_racon
                    .combine(no_meta_ont)
                    .set{ch_racon_full}
            } else if (params.ONT_lr  == false && params.PacBioHifi_lr == true){
                ch_racon
                    .combine(no_meta_pb)
                    .set{ch_racon_full}
            }
            RACON (ch_racon_full)
                
            no_meta_polished_assembly = RACON.out.improved_assembly

            no_meta_polished_assembly
                .map { file -> tuple(id: file.baseName.replaceAll(/\\.fasta$/, ''), file) }
                .set { racon_polished_assembly }

        ch_versions = ch_versions.mix(RACON.out.versions)
        }
        
        if (params.medaka_polish == true && params.racon_polish == true) {
            println "polishing racon-polished assemblies with medaka!"
            MEDAKA (fastq_filt, racon_polished_assembly, model)
        ch_versions = ch_versions.mix(MEDAKA.out.versions)
 
            assembly_polished = MEDAKA.out.assembly
        
            assembly_polished
                .map { file -> tuple(id: file.baseName, file)  }
                .set { medaka_assembly_polished }

            assembly_polished
                .concat(no_meta_polished_assembly)
                .set{no_meta_polished_assembly}
            
        } else if (params.medaka_polish == true && params.racon_polish == false){

            println "polishing assemblies with medaka!"
            MEDAKA (fastq_filt, assembly, model)
        ch_versions = ch_versions.mix(MEDAKA.out.versions)
 
        no_meta_polished_assembly = MEDAKA.out.assembly

        no_meta_polished_assembly
                .map { file -> tuple(id: file.baseName, file)  }
                .set { assembly_polished }
        } 
        
    emit:
    
        no_meta_polished_assembly      
        
    versions = ch_versions                     // channel: [ versions.yml ]
}
