process HIFIASM {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::hifiasm=0.18.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hifiasm:0.19.8--h43eeafb_0' :
        'biocontainers/hifiasm:0.19.8--h43eeafb_0' }"

    input:
    tuple val(meta), path(hifi_reads)
    path ont


    output:
    tuple val(meta), path("*.r_utg.gfa")       , emit: raw_unitigs
    tuple val(meta), path("*.ec.bin")          , emit: corrected_reads
    tuple val(meta), path("*.ovlp.source.bin") , emit: source_overlaps
    tuple val(meta), path("*.ovlp.reverse.bin"), emit: reverse_overlaps
    tuple val(meta), path("*.p_utg.gfa")       , emit: processed_unitigs, optional: true
    tuple val(meta), path("*.asm.p_ctg.gfa")   , emit: primary_contigs  , optional: true
    tuple val(meta), path("*.asm.a_ctg.gfa")   , emit: alternate_contigs, optional: true
    path("*.fasta")           , emit: assembly_fasta
    path  "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if(ont){
       """
        hifiasm \\
            $args \\
            -o ${prefix}.asm \\
            --primary \\
            -t $task.cpus \\
            --ul $ont \\
            $hifi_reads

        awk '/^S/{print ">"\$2;print \$3}' *.p_ctg.gfa > hifiasm_${prefix}.asm.p_ctg.fasta

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            hifiasm: \$(hifiasm --version 2>&1)
        END_VERSIONS
        """ 
    } else {
        """
        hifiasm \\
            $args \\
            -o ${prefix}.asm \\
            --primary \\
            -t $task.cpus \\
            $hifi_reads

        awk '/^S/{print ">"\$2;print \$3}' *.p_ctg.gfa > hifiasm_${prefix}.asm.p_ctg.fasta

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            hifiasm: \$(hifiasm --version 2>&1)
        END_VERSIONS
    
    """}
}
