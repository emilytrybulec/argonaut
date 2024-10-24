process REDUNDANS_P {
    tag "$meta"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/redundans:2.01--py310pl5321h43eeafb_0' :
        'biocontainers/redundans:2.01--py310pl5321h43eeafb_0' }"


    input:
    tuple val(meta), path(fasta), path(shortread1), path(shortread2)  

    output:
    path("*redundans_purge/*redundans_purge_contig.fasta")               , emit: assembly_fasta

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    redundans.py -v -f $fasta -i $shortread1 $shortread2 -t $task.cpus --nogapclosing --noscaffolding -o ${prefix}_redundans_purge
    
    cd ${prefix}_redundans_purge
    mv scaffolds.reduced.fa ${prefix}_redundans_purge_scaf.fasta
    mv contigs.reduced.fa ${prefix}_redundans_purge_contig.fasta
    """
}
