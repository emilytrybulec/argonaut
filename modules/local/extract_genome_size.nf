process EXTRACT_LR {
    tag "$meta.id"
    label 'process_low'

    container 'emilytrybulec/genassembly:kmer'

    input:
    tuple val(meta), path(gce2log)

    output:
    path("shortenedSizeFinal.txt")        , emit: genome_size_est
    path("standardSizeFinal.txt")     , emit: standard_fmt_est    

    script: 
    def number
    """
    tail -3 $gce2log | awk '{print \$5}' > scientificSize.txt
    less scientificSize.txt | awk -F"E" 'BEGIN{OFMT="%10.10f"} {print \$1 * (10 ^ \$2)}' > standardSize.txt 
    
    head -1 standardSize.txt | awk '{printf("%d\\n",\$0+=\$0<0?-0.5:0.5)}' > standardSizeFinal.txt
    
    head -1 standardSizeFinal.txt | numfmt --to=si --format %.2f --round=nearest > shortenedSizeFinal.txt
    """
}
