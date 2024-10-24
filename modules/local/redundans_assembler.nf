process REDUNDANS_A {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/redundans:2.01--py310pl5321h43eeafb_0' :
        'biocontainers/redundans:2.01--py310pl5321h43eeafb_0' }"


    input:
    tuple val(meta), path(shortreads)

    output:
    tuple val(meta), path("*filled.fa")               , emit: assembly_fasta

    script:
    """
    redundans.py -v -i $shortreads -t $task.cpus --noreduction --nogapclosing --noscaffolding -o redundans_assembly

    """
}
