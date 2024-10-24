process SEQKIT_GREP {
    tag "$meta.id"
    label 'process_low'


    conda "bioconda::seqkit=2.4.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.4.0--h9ee0642_0':
        'biocontainers/seqkit:2.4.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(centrifuge_out)
    tuple val(meta), path(longreads)

    output:
    path("*.fastq")  , emit: filter
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    grep -vw "unclassified" $centrifuge_out > contaminated_reads.txt
    awk NF=1 contaminated_reads.txt > contaminated_read_ids.txt
    sort -u contaminated_read_ids.txt > no_dup_centrifuge_contaminated_read_ids.txt

    seqkit grep -v -f no_dup_centrifuge_contaminated_read_ids.txt $longreads > dc_${prefix}.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$( seqkit version | sed 's/seqkit v//' )
    END_VERSIONS
    """
}
