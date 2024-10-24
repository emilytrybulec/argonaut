process TOTAL_BASES_LR {
    label 'process_low'
    tag "$meta.id"

    input:
    tuple val(meta), path(nanoplot_report)

    output:
    tuple val(meta), path("*totalBasesLR.txt")        , emit: total_bases
    tuple val(meta), path("*totalBasesLR_pretty.txt")        , emit: pretty_total_bases

    script: 
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    sed -n '9p' < $nanoplot_report | awk '{print "\\""\$3"\\""}' | sed -e 's/^"//' -e 's/"\$//' | sed 's/,//g' > ${prefix}_totalBasesLR.txt
    head -1 ${prefix}_totalBasesLR.txt | numfmt --to=si > ${prefix}_totalBasesLR_pretty.txt
    """
}
