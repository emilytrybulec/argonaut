include { POLCA } from '../../modules/local/polca'  
include { NEXTPOLISH2 } from '../../modules/local/nextpolish2'

workflow POLISH2 {

    take:
        flye_assembly
        shortreads
        hifi_reads
        
    main:

    ch_versions = Channel.empty() 

    if (hifi_reads != null && !hifi_reads.isEmpty()) {
        println "polishing assemblies with short reads using Nextpolish2!"

        flye_assembly
            .combine(hifi_reads)
            .combine(shortreads)
            .set { ch_nextpolish2 }

        NextPolish2 (ch_nextpolish2)

        flye_assembly_polished = NextPolish2.out.np2_polished_assembly
        ch_versions = NextPolish2.out.versions
    } else {
        println "polishing assemblies with short reads using POLCA!"

        flye_assembly
            .combine(shortreads)
            .set { ch_polca) }

        POLCA (ch_polca)

        flye_assembly_polished = POLCA.out.sr_polished_assembly
        ch_versions = POLCA.out.versions
    }

    emit:
    
        flye_assembly_polished
        versions = ch_versions
}
