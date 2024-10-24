process FORMAT {
    label 'process_low'

    input:
    val genome_size

    output:
    path("shortenedSizeFinal.txt")        , emit: genome_size_est
    path("standardSize.txt")     , emit: standard_fmt_est    

    script: 
    """
    echo $genome_size > standardSize.txt 
    
    head -1 standardSize.txt | numfmt --to=si > shortenedSizeFinal.txt
    """
}