process COVERAGE_SR {
    label 'process_low'

    input:
    val genome_size
    tuple val(meta), path(total_bases_before)
    tuple val(meta), path(total_bases_after)

    output:
    tuple val(meta), path("estimatedCoverageSR_before.txt")        , emit: coverage_est_before
    tuple val(meta), path("estimatedCoverageSR_after.txt")        , emit: coverage_est_filtered


    script: 
    def est_size
    def bases_before
    def bases_after
    """
    est_size=\$(grep -o '[0-9]\\+' $genome_size)
    bases_before=\$(grep -o '[0-9]\\+' $total_bases_before)
    bases_after=\$(grep -o '[0-9]\\+' $total_bases_after)

    printf "%.2f\n" \$((10**2 * \$bases_before/\$est_size))e-2 > estimatedCoverageSR_before.txt
    printf "%.2f\n" \$((10**2 * \$bases_after/\$est_size))e-2 > estimatedCoverageSR_after.txt

    """
}
