include { MEDAKA } from '../../modules/nf-core/medaka/main'  
include { RACON } from '../../modules/nf-core/racon/main'
include { PILON } from '../../modules/nf-core/pilon/main'

workflow POLISH {

    take:
        assembly
        fastq_filt
        model
        paf
        ch_racon
        no_meta_ont
        no_meta_pb
        alignments
        pilon_mode
    main:

    ch_versions = Channel.empty() 

        if (params.racon_polish == true) {
            println "polishing assemblies with racon!"

            if (params.ONT_lr  == true && params.PacBioHifi_lr == true){
                RACON (ch_racon, no_meta_ont, no_meta_pb)
            } else if (params.ONT_lr  == true && params.PacBioHifi_lr == false){
                RACON (ch_racon, no_meta_ont, [])
            } else if (params.ONT_lr  == false && params.PacBioHifi_lr == true){
                RACON (ch_racon, [], no_meta_pb)
            }
                
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

        if (params.pilon_polish == true && params.racon_polish == true && params.medaka_polish == false) { 
           println "polishing racon-polished assemblies with pilon!"
           racon_polished_assembly.join(alignments)
             .combine(pilon_mode)
             .set{pilon_input}

           PILON(pilon_input)
 
        ch_versions = ch_versions.mix(PILON.out.assembly)
        
           pilon_improved_assembly = PILON.out.improved_assembly

           pilon_improved_assembly
               .map { file -> tuple(id: file.baseName, file)  }
               .set { pilon_polished_assembly } 
       
           pilon_improved_assembly 
               .concat(no_meta_polished_assembly)
               .set{no_meta_polished_assembly}
      
       } else if (params.pilon_polish == true && params.racon_polish == false && params.medaka_polish == true) {
           println "polishing medaka-polished assemblies with pilon!"
           assembly_polished.join(alignments)
             .combine(pilon_mode)
             .set{pilon_input}
          
           PILON(pilon_input)
        
        ch_versions = ch_versions.mix(PILON.out.assembly)
        
           pilon_improved_assembly = PILON.out.improved_assembly
        
           pilon_improved_assembly
               .map { file -> tuple(id: file.baseName, file)  }
               .set { pilon_polished_assembly }
  
           pilon_improved_assembly
               .concat(no_meta_polished_assembly)
               .set{no_meta_polished_assembly}

       } else if (params.pilon == true && params.racon_polish == true && params.medaka_polish == true) {
           println "polishing racon-and-medaka-polished assemblies with pilon!"
           medaka_assembly_polished.join(alignments) 
             .combine(pilon_mode)
             .set{pilon_input}
             
           PILON(pilon_input)
           
        ch_versions = ch_versions.mix(PILON.out.assembly)
        
           pilon_improved_assembly = PILON.out.improved_assembly

           pilon_improved_assembly
               .map { file -> tuple(id: file.baseName, file)  }
               .set { pilon_polished_assembly }

           pilon_improved_assembly
               .concat(no_meta_polished_assembly)
               .set{no_meta_polished_assembly}

        } else if (params.pilon == true && params.racon_polish == false && params.medaka_polish == false) { 
            println "polishing assemblies with pilon!"
            assembly.join(alignments) 
             .combine(pilon_mode)
             .set{pilon_input}
             
           PILON(pilon_input)

        ch_versions = ch_versions.mix(PILON.out.assembly)

           pilon_improved_assembly = PILON.out.improved_assembly

           pilon_improved_assembly 
               .map { file -> tuple(id: file.baseName, file)  }
               .set { pilon_polished_assembly }
        }
        
    emit:
    
        no_meta_polished_assembly      
        
    versions = ch_versions                     // channel: [ versions.yml ]
}
