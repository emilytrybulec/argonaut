process MASURCA {
    tag "$meta.id"
    label 'process_high_memory', 'error_ignore'

    container 'staphb/masurca'
    
    input:
    tuple val(meta), path(longreads) //path_to/longreads.gz
    tuple val(meta), path(shortreads) //path_to/pe_R1.fa /path_to/pe_R2.fa

    output:
    path("hybrid_masurca*")                , emit: fasta
    path ("versions.yml")                , emit: versions

    script:
    def VERSION = '4.1.0'
    def sr
    def prefix = task.ext.prefix ?: "masurca_${meta.id}"
    """
    sr=\$(echo '${shortreads}' | sed -e "s/ /,/g")
    masurca -t $task.cpus -i \$sr -r $longreads

    mv CA*/primary.genome.scf.fasta hybrid_${prefix}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MaSuRCA: $VERSION
    END_VERSIONS
    """
}
