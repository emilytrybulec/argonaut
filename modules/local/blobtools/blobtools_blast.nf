process BLOBTOOLS_BLAST {
    tag "$meta.id"
    label 'process_high'

    container 'ncbi/blast'

    input:
    tuple val(meta), path(assembly)

    output:
    path('*.out'), emit: hits

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    blastn -db ${params.blast_db} \\
       -query $assembly \\
       -outfmt "6 qseqid staxids bitscore std" \\
       -max_target_seqs 10 \\
       -max_hsps 1 \\
       -evalue 1e-25 \\
       -num_threads $task.cpus \\
       -out ${prefix}_blast_hits.out

    """
}    
