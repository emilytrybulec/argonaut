include { POLCA } from '../../modules/local/polca'  
include { NextPolish2 } from '../../modules/local/nextpolish2'

workflow POLISH2 {

    take:
        flye_assembly
        shortreads
        
    main:

    ch_versions = Channel.empty() 

        println "polishing assemblies with short reads using POLCA!"

        POLCA (flye_assembly, shortreads)


    emit:
    
        flye_assembly_polished      = POLCA.out.sr_polished_assembly           
        
    versions = ch_versions                     // channel: [ versions.yml ]
}
