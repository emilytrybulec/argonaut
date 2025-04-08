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
    def herro = "${nanoplot_report}"
    def auto_herro = herro.contains('herro') ? '7p' : herro.contains('hifiasm') ? '7p' : '9p'
    """
    sed -n '${auto_herro}' < $nanoplot_report | awk '{print "\\""\$3"\\""}' | sed -e 's/^"//' -e 's/"\$//' | sed 's/,//g' > ${prefix}_totalBasesLR.txt
    head -1 ${prefix}_totalBasesLR.txt | numfmt --to=si > ${prefix}_totalBasesLR_pretty.txt
    """
}
