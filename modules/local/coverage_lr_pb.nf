process COVERAGE_LR_PB {
    label 'process_low'
    tag "$meta.id"

    input:
    val genome_size
    tuple val(meta), path(total_bases)

    output:
    tuple val(meta), path("*estimatedCoverageLR.txt")        , emit: coverage_est

    script: 
    def est_size
    def est_bases
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    est_size=\$(grep -o '[0-9]\\+' $genome_size)
    est_bases=\$(grep -o '[0-9]\\+' $total_bases | head -1)

    printf "%.2f\n" \$((10**2 * \$est_bases/\$est_size))e-2 > ${prefix}_estimatedCoverageLR.txt
    """
}
