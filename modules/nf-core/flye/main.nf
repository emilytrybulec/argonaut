process FLYE {
    tag "$meta.id"
    label 'process_high_memory'

    conda "bioconda::flye=2.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/flye:2.9--py39h6935b12_1' :
        'biocontainers/flye:2.9--py39h6935b12_1' }"

    input:
    tuple val(meta), path(reads), path(genome_size_est)

    output:
    path("*.fasta"), emit: fasta
    tuple val(meta), path("*.gfa")  , emit: gfa
    tuple val(meta), path("*.gv")   , emit: gv
    tuple val(meta), path("*.txt")     , emit: txt
    tuple val(meta), path("*.log")     , emit: log
    tuple val(meta), path("*.json")    , emit: json
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "flye_${meta.id}"
    def size
    def read_name = "${reads}"
    def auto_ont_mode = read_name.contains('ont') ? '--nano-raw' : ''
    def auto_pb_mode = read_name.contains('pb') ? '--pacbio-hifi' : ''
   
    """
    size=\$(echo "\$(<${genome_size_est})")
    flye \\
        ${auto_ont_mode} \\
        ${auto_pb_mode} \\
        $reads \\
        --out-dir . \\
        --threads \\
        $task.cpus \\
        --genome-size \$size \\
        $args

    mv assembly.fasta ${prefix}.fasta
    mv assembly_info.txt ${prefix}.assembly_info.txt
    mv flye.log ${prefix}.flye.log
    mv params.json ${prefix}.params.json
    mv assembly_graph.gfa ${prefix}_assembly_graph.gfa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flye: \$( flye --version )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo stub > assembly.fasta | gzip -c assembly.fasta > ${prefix}.assembly.fasta.gz
    echo stub > assembly_graph.gfa | gzip -c assembly_graph.gfa > ${prefix}.assembly_graph.gfa.gz
    echo stub > assembly_graph.gv | gzip -c assembly_graph.gv > ${prefix}.assembly_graph.gv.gz
    echo contig_1 > ${prefix}.assembly_info.txt
    echo stub > ${prefix}.flye.log
    echo stub > ${prefix}.params.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flye: \$( flye --version )
    END_VERSIONS
    """
}
