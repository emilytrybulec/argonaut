process WINNOWMAP {
    tag "$meta.id"
    label 'process_medium'

    // Note: the versions here need to match the versions used in the mulled container below and minimap2/index

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/winnowmap:2.03--h5b5514e_1' :
        'quay.io/biocontainers/winnowmap:2.03--h5b5514e_1' }"

    input:
    tuple val(meta), path(reference), path(reads), path(repetitive_txt)
    val(kmernum)

    output:
    tuple val(meta), path("*.sam"), emit: sam

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def read_name = "${reads}"
    def auto_ont_mode = read_name.contains('ont') ? 'map-ont' : ''
    def auto_pb_mode = read_name.contains('pb') ? 'map-pb' : ''
    """
    winnowmap \\
      -W $repetitive_txt \\
      -k $kmernum \\
      -ax ${auto_pb_mode}${auto_ont_mode} \\
      $reference \\
      $reads > ${prefix}.sam
    """
}
