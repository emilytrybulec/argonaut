process MASURCA_LR_ADV {
    label 'process_high_memory', 'error_ignore'

    container 'staphb/masurca'
    
    input:
    path config

    output:
    path("masurca*")                , emit: fasta
    path ("*/versions.yml")                , emit: versions

    script:
    def VERSION = '4.1.0'
    def prefix = task.ext.prefix ?: "masurca"
    """
    masurca $config
    ./assemble.sh > masurca_advanced.log

    mv CA*/primary.genome.scf.fasta ${prefix}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MaSuRCA: $VERSION
    END_VERSIONS
    """
}
