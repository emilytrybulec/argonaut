process POLCA {
tag "$meta.id"
    label 'process_medium'

    container 'staphb/masurca:4.1.0'

    input:
    tuple val(meta), path(assembly) //path_to/assembly.fasta
    tuple val(meta), path(shortreads) //path_to/pe_R1.fa,/path_to/pe_R2.fa

    output:
    path("*.fasta")                    , emit: sr_polished_assembly
    path "versions.yml"                , emit: versions

    script:
    def VERSION = '4.1.0'
    def prefix = task.ext.prefix ?: "polca_${meta.id}"
    """
    polca.sh -t 6 -a $assembly -r '$shortreads'

    mv *PolcaCorrected.fa ${prefix}.fasta
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        POLCA/MaSuRCA: $VERSION
    END_VERSIONS
    """
}