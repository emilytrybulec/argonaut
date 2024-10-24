process BLOBTOOLS_CONFIG {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path('*config.yaml'), emit: config

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def configContent = """
    assembly:
      accession: ${meta.id}
      file: $assembly
      level: scaffold
      prefix: ${meta.id}
    busco:
      download_dir: ${params.busco_lineage}
      lineages:
        - ${params.busco_lineage}
      basal_lineages:
        - ${params.busco_lineage}
        
    """

    // Write the combined configuration to the output file
    """
    echo '$configContent' > ${prefix}_config.yaml
    """
}
