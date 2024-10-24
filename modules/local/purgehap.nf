process ALIGN {
    tag "$meta.id"
    label 'process_medium'
    time '36h'
    // Note: the versions here need to match the versions used in the mulled container below and minimap2/index
    conda "bioconda::minimap2=2.24 bioconda::samtools=1.14"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:1679e915ddb9d6b4abda91880c4b48857d471bd8-0' :
        'biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:1679e915ddb9d6b4abda91880c4b48857d471bd8-0' }"

    input:
    tuple val(meta), path(reference), path(reads)
    val bam_format
    val cigar_paf_format
    val cigar_bam

    output:
    tuple val(meta), path("*.paf"), optional: true, emit: paf
    tuple val(meta), path("*.bam"), optional: true, emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bam_output = bam_format ? "-a | samtools sort | samtools view -@ ${task.cpus} -b -h -o ${prefix}.bam" : "-o ${prefix}.paf"
    def cigar_paf = cigar_paf_format && !bam_format ? "-c" : ''
    def set_cigar_bam = cigar_bam && bam_format ? "-L" : ''
    """
    minimap2 \\
        $args \\
        -t $task.cpus \\
        "${reference ?: reads}" \\
        "$reads" \\
        $cigar_paf \\
        $set_cigar_bam \\
        $bam_output


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """
}

process HISTOGRAM {
    tag "$meta"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/purge_haplotigs:1.1.2--hdfd78af_0' :
        'biocontainers/purge_haplotigs:1.1.2--hdfd78af_0' }"

    input:
    tuple val(meta), path(contigs), path(aligned_bam)

    output:
    tuple val(meta), path ("*.bam.gencov"), emit: gencov
    tuple val(meta), path ("*.bam.histogram.png"), emit: hist

    """
    purge_haplotigs hist -t $task.cpus -b $aligned_bam -g $contigs
    """
}

process PURGE {
    tag "$meta"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/purge_haplotigs:1.1.2--hdfd78af_0' :
        'biocontainers/purge_haplotigs:1.1.2--hdfd78af_0' }"
        
    input:
    val low
    val mid
    val high
    tuple val(meta), path(assembly), path(gencov)

    output:
    path("purge_*"), emit: purged

    """
    purge_haplotigs cov -in $gencov \
        -low $low -mid $mid -high $high -o ${assembly.baseName}.coverage_stats.csv
    purge_haplotigs purge -t $task.cpus -g $assembly -c ${assembly.baseName}.coverage_stats.csv

    mv curated.fasta purge_${assembly.baseName}.fasta
    """
}
