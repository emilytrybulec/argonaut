include { RAGTAG } from '../../modules/local/ragtag' 

workflow SCAFFOLD {

    take:
        assemblies // channel: val(meta), path(polished assembly alignment to reads)
        reference

    main:

    ch_versions = Channel.empty()
        println "scaffolding assemblies with ragtag!"

        //optional scaffolding with the same species or most closely related species available
        RAGTAG (assemblies, reference)
        scaffolded_assemblies = RAGTAG.out.scaffolded_assembly           

        scaffolded_assemblies
                .map { file -> tuple([id: file.baseName], file)  }
                .set { assembly_polished_purged_scaffolded }

        
    emit:
    
        assembly_polished_purged_scaffolded      
        scaffolded_assemblies
        
    versions = ch_versions                     // channel: [ versions.yml ]
}
