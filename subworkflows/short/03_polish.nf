include { POLCA } from '../../modules/local/polca'  

workflow POLISH2 {

    take:
        flye_assembly
        shortreads
        
    main:

    ch_versions = Channel.empty() 

    if(params.polca == true){
        println "polishing assemblies with short reads using POLCA!"

        flye_assembly
            .combine(shortreads)
            .set{ch_polca}

        POLCA (ch_polca) } else { flye_assembly_polished = Channel.empty()}


    emit:
    
        flye_assembly_polished      = POLCA.out.sr_polished_assembly           
        
    versions = ch_versions                     // channel: [ versions.yml ]
}
