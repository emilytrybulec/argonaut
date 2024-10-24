process PILON {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container 'staphb/pilon'

    input:
    tuple val(meta), path(fasta), path(bam), path(bai)

    output:
    path("*.fasta") , emit: improved_assembly
    tuple val(meta), path("*.vcf")   , emit: vcf               , optional : true
    tuple val(meta), path("*.change"), emit: change_record     , optional : true
    tuple val(meta), path("*.bed")   , emit: tracks_bed        , optional : true
    tuple val(meta), path("*.wig")   , emit: tracks_wig        , optional : true
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "pilon_${meta.id}"
    def reads = "${fasta}"
    def auto_bam_mode = reads.contains('ont') ? '--nanopore' : reads.contains('pb') ? '--pacbio' : '--bam'
    """
    java -Xmx500G -jar /pilon/pilon.jar \\
        --genome $fasta \\
        --output ${prefix} \\
        $args \\
        $auto_bam_mode $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pilon: \$(echo \$(pilon --version) | sed 's/^.*version //; s/ .*\$//' )
    """
}
