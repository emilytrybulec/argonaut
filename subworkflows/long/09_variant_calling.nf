include { LONGSHOT } from '../../modules/local/longshot'

workflow VARIANT_CALLING {

    take:
        alignments // channel: val(meta), path(polished assembly bam)
        reference 

    main: 

    ch_versions = Channel.empty()
        println "variant calling with longshot!"

    alignments = alignments
        .map { tuple(meta, bam) -> tuple(meta, bam, bam + ".bai") }
        .each { meta, bam, bai -> "samtools index $bam".execute().waitFor() }
        .map { meta, bam -> tuple(meta, bam, bam + ".bai") }
        .map { tuple(meta, bam, bai) -> tuple(meta, bam, reference) }
        .set {longshot_input}

    LONGSHOT(longshot_input)
    vcf = LONGSHOT.out.vcf
    hap_bam = LONGSHOT.out.hap_bam

    // process outputs 
    vcf 
        .map { meta, file -> tuple(meta, file) }
        .set { phased_variants }

    hap_bam
        .map { meta, bam, bai -> tuple(meta, bam, bai) }
        .set { haplotype_bam }

    emit: 
        phased_variants
        haplotype_bam

    versions = ch_versions                    // channel: [ versions.yml ]
}
