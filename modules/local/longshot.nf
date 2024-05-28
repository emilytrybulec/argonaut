process LONGSHOT {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/longshot:1.0.0--hd4f2111_1' :
        'quay.io/biocontainers/longshot:1.0.0--hd4f2111_1' }"

    input:
    tuple val(meta), path(bam_file)
    path reference_fasta

    output:
    path "output.vcf", emit: vcf
    path "haplotype_seperated.bam", emit: hap_bam, optional: true
    path "versions.yml", emit: versions

    script:
    def VERSION = '1.0.0'

   """
    longshot \\
        --bam ${bam_file} \\
        --ref ${reference_fasta} \\
        --out output.vcf 

  \\ to get haplotype-seperated BAM
  longshot --bam ${bam)file} --ref ${reference_fasta} --out output.vcf --out_bam haplotype_seperated.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Longshot: $VERSION
    END_VERSIONS
    """
}
