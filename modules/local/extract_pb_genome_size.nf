process EXTRACT_PB {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(genomescope_report)

    output:
    path("shortenedSizeFinal.txt")        , emit: genome_size_est
    path("standardSizeFinal.txt")     , emit: standard_fmt_est    

    script: 
    def number
    """
    sed -n '11p' < $genomescope_report | awk '{print \$6}' | sed 's/,//g' > standardSize.txt
    
    head -1 standardSize.txt > standardSizeFinal.txt
    
    head -1 standardSizeFinal.txt | numfmt --to=si --format %.2f --round=nearest > shortenedSizeFinal.txt
    """
}
