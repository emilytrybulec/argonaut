process VERKKO {
    tag "$meta.id"
    label 'process_high_memory', 'error_ignore'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/verkko:2.2.1--h45dadce_0' :
        'quay.io/biocontainers/verkko:2.2.1--h45dadce_0' }"
        
    input:
    tuple val(meta), path(pb)
    path(ont)
    path(ref)

    output:
    path("verkko*/*${meta.id}.fasta")        , emit: fasta
    tuple val(meta), path("verkko*/*.gfa")                    , emit: gfa

    script:
    def VERSION = '4.1.0'
    def prefix = task.ext.prefix ?: "verkko_${meta.id}"
if(ont){
    if(ref){
    """
    mkdir cache
    export XDG_CACHE_HOME=./cache
    verkko -d ${prefix} \\
    --hifi $pb \\
    --nano $ont \\
    --ref $ref \\
    --threads $task.cpus

    cd ${prefix}
    mv assembly.fasta ${prefix}.fasta
    mv assembly.homopolymer-compressed.gfa ${prefix}_homopolymer-compressed.gfa
    
    """
    } else {
        """
    mkdir cache
    export XDG_CACHE_HOME=./cache
    verkko -d ${prefix} \\
    --hifi $pb \\
    --nano $ont \\
    --threads $task.cpus

    cd ${prefix}
    mv assembly.fasta ${prefix}.fasta
    mv assembly.homopolymer-compressed.gfa ${prefix}_homopolymer-compressed.gfa
    
    """
    }
} else {
    if(ref){
    """
    mkdir cache
    export XDG_CACHE_HOME=./cache
    verkko -d ${prefix} \\
    --hifi $pb \\
    --ref $ref \\
    --threads $task.cpus

    cd ${prefix}
    mv assembly.fasta ${prefix}.fasta
    mv assembly.homopolymer-compressed.gfa ${prefix}_homopolymer-compressed.gfa
    
    """ } else {
    """
    mkdir cache
    export XDG_CACHE_HOME=./cache
    verkko -d ${prefix} \\
    --hifi $pb \\
    --threads $task.cpus

    cd ${prefix}
    mv assembly.fasta ${prefix}.fasta
    mv assembly.homopolymer-compressed.gfa ${prefix}_homopolymer-compressed.gfa
    """
    }
}

    
}
