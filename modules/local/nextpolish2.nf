process NEXTPOLISH2 {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nextpolish2:0.2.0--hdcf5f25_0' :
        'quay.io/biocontainers/nextpolish2:0.2.0--hdcf5f25_0' }"

   \\  STILL NEED TO EDIT THIS  
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
