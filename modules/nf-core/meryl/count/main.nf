process MERYL_COUNT {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::meryl=1.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/meryl:1.3--h87f3376_1':
        'biocontainers/meryl:1.3--h87f3376_1' }"

    input:
    tuple val(meta), path(reads)
    val kmer

    output:
    path("*filtered.meryl")           , emit: meryl_db
    path("*.txt")           , emit: repetitive_k
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def kmernum = "k=${kmer}"
    """
    for READ in $reads; do
        meryl count \\
            threads=$task.cpus \\
            $kmernum \\
            $args \\
            $reads \\
            output kmer_db.meryl

        meryl greater-than 1 \\
            threads=$task.cpus \\
            $kmernum \\
            output kmer_db.filtered.meryl kmer_db.meryl

        meryl print greater-than distinct=0.9998 kmer_db.filtered.meryl > repetitive_k${kmer}.txt

    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        meryl: \$( meryl --version |& sed 's/meryl //' )
    END_VERSIONS
    """
}
