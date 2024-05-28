process VERKKO {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/verkko:2.0--h45dadce_0' :
        'quay.io/biocontainers/verkko:2.0--h45dadce_0' }"

    input:
    tuple val(meta), path(hifi_reads)
    path ont_reads, optional: true

    output:
    path "asm/assembly.fasta", emit: final_assembly
    path "asm/assembly.homopolymer-compressed.gfa", emit: final_graph
    path "asm/assembly*csv", emit: coverage_files
    path "asm/assembly.scfmap", emit: scaffold_map
    path "versions.yml", emit: versions

    script:
    def VERSION = '2.0'

   """
    verkko -d asm --hifi ${hifi_reads} ${ont_reads ? '--nano ' + ont_reads : ''}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Verkko: $VERSION
    END_VERSIONS
    """
}
