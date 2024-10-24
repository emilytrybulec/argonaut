process CANU {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::canu=2.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/canu:2.2--ha47f30e_0':
        'biocontainers/canu:2.2--ha47f30e_0' }"

    input:
    tuple val(meta), path(reads)
    path(genome_size_est)

    output:
    tuple val(meta), path("*.report")                   , emit: report
    tuple val(meta), path("*.contigs.fasta.gz")         , emit: assembly                , optional: true
    tuple val(meta), path("*.unassembled.fasta.gz")     , emit: contigs
    tuple val(meta), path("*.correctedReads.fasta.gz")	, emit: corrected_reads         , optional: true
    tuple val(meta), path("*.trimmedReads.fasta.gz")	, emit: corrected_trimmed_reads , optional: true
    tuple val(meta), path("*.contigs.layout")           , emit: metadata                , optional: true
    tuple val(meta), path("*.contigs.layout.readToTig") , emit: contig_position         , optional: true
    tuple val(meta), path("*.contigs.layout.tigInfo")   , emit: contig_info             , optional: true
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "canu_${meta.id}"
    def genomesize 
    def read_name = "${reads}"
    def auto_ont_mode = read_name.contains('ont') ? '-nanopore' : ''
    def auto_pb_mode = read_name.contains('pb') ? '-pacbio-hifi' : ''

    """
    genomesize=\$(echo "\$(<${genome_size_est})")

    canu \\
        -p ${prefix} \\
        genomeSize=\$genomesize \\
        $args \\
        maxThreads=$task.cpus \\
        ${auto_ont_mode} \\
        ${auto_pb_mode} \\ 
        $reads
        

    gzip *.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        canu: \$(echo \$(canu --version 2>&1) | sed 's/^.*canu //; s/Using.*\$//' )
    END_VERSIONS
    """
}
