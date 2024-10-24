process GCE {
    tag "$meta.id"
    label 'process_medium', 'error_ignore'

    container 'emilytrybulec/genassembly:kmer'

    input:
    tuple val(meta), path(kmerstat)
    tuple val(meta), path(kmernum)

    output:
    tuple val(meta), path("*.log")                 , emit: gce2log
    path  "versions.yml"                           , emit: versions
   
    script: 
    def VERSION = '1.0.2'
    def number
    """
    grep -Po '\\d+' $kmernum > extractedNum.txt
    number=\$(grep -o '[0-9]\\+' extractedNum.txt)
    /GCE-master/gce-1.0.2/gce -g \$number -f $kmerstat >gce2.table 2>gce2.log
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gce: $VERSION
    END_VERSIONS
    """
}
