process MASURCA_SR_ADV {
    label 'process_high_memory', 'error_ignore'

    container 'staphb/masurca'
    
    input:
    path config

    output:
    path("SOAP_a*/masurca*")                , emit: fasta
    path ("SOAP_a*/versions.yml")                , emit: versions

    script:
    def VERSION = '4.1.0'
    def prefix = task.ext.prefix ?: "masurca_sr_only"
    """
    masurca $config
    ./assemble.sh > masurca_advanced.log

    cd SOAP_assembly
    mv asm2.scafSeq2 ${prefix}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MaSuRCA: $VERSION
    END_VERSIONS
    """
}
