process RUN_QC {
    tag "$meta.id"
    label 'process_low'

    container 'kathrin1414/runqc'

    input:
    path dataset_xml

    output:
    tuple val(meta), path("*.gz")        , emit: untar

    script: 
    def args = task.ext.args ?: ''
    """
    runqc-reports $dataset_xml 
    """
}
