process NEXTPOLISH2 {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nextpolish2:0.2.0--hdcf5f25_0' :
        'quay.io/biocontainers/nextpolish2:0.2.0--hdcf5f25_0' }"
  
input:
    tuple val(meta), path(assembly)
    tuple val(meta), path(hifi)
    tuple val(meta), path(shortreads)
    

    output:
    path("*.np2.fa")              , emit: np2_polished_assembly
    path "versions.yml"           , emit: versions

    script:
    def VERSION = '0.2.0'

   """
   // prepare HiFi mapping file using winnowmap
    meryl count k=15 output merylDB ${assembly}
    meryl print greater-than distinct=0.9998 merylDB > repetitive_k15.txt
    winnowmap -t 5 -W repetitive_k15.txt -ax map-pb ${assembly} ${hifi} | samtools sort -o hifi.map.sort.bam -

    // index sorted BAM file
    samtools index hifi.map.sort.bam

    // prepare k-mer dataset files (21-mer and 31-mer)
    yak count -o k21.yak -k 21 -b 37 <(zcat ${shortreads[0]}) <(zcat ${shortreads[1]})
    yak count -o k31.yak -k 31 -b 37 <(zcat ${shortreads[0]}) <(zcat ${shortreads[1]})

    // run nextpolish2
    nextPolish2 -t 5 hifi.map.sort.bam ${assembly} k21.yak k31.yak > np2_${meta.id}.np2.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        NextPolish2: $VERSION
    END_VERSIONS
    """
}
