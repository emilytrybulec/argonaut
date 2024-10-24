process MERQURY {
    tag "$meta"
    label 'process_low'

    conda "bioconda::merqury=1.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/merqury:1.3--hdfd78af_1':
        'biocontainers/merqury:1.3--hdfd78af_1' }"

    input:
    tuple val(meta), path(assembly), path(meryl_db), path(genome_size_est)
    val tolerable_collision

    output:
    tuple val(meta), path("*_only.bed")          , emit: assembly_only_kmers_bed
    tuple val(meta), path("*_only.wig")          , emit: assembly_only_kmers_wig
    tuple val(meta), path("*.completeness.stats"), emit: stats
    tuple val(meta), path("*.dist_only.hist")    , emit: dist_hist
    tuple val(meta), path("*.spectra-cn.fl.png") , emit: spectra_cn_fl_png
    tuple val(meta), path("*.spectra-cn.hist")   , emit: spectra_cn_hist
    tuple val(meta), path("*.spectra-cn.ln.png") , emit: spectra_cn_ln_png
    tuple val(meta), path("*.spectra-cn.st.png") , emit: spectra_cn_st_png
    tuple val(meta), path("*.spectra-asm.fl.png"), emit: spectra_asm_fl_png
    tuple val(meta), path("*.spectra-asm.hist")  , emit: spectra_asm_hist
    tuple val(meta), path("*.spectra-asm.ln.png"), emit: spectra_asm_ln_png
    tuple val(meta), path("*.spectra-asm.st.png"), emit: spectra_asm_st_png
    tuple val(meta), path("${prefix}.qv")        , emit: assembly_qv
    tuple val(meta), path("${prefix}.*.qv")      , emit: scaffold_qv
    tuple val(meta), path("*.hist.ploidy")       , emit: read_ploidy
    tuple val(meta), path("best_kmer_num.txt")   , emit: kmer
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = 1.3
    def genome_size
    """
    # Nextflow changes the container --entrypoint to /bin/bash (container default entrypoint: /usr/local/env-execute)
    # Check for container variable initialisation script and source it.
    if [ -f "/usr/local/env-activate.sh" ]; then
        set +u  # Otherwise, errors out because of various unbound variables
        . "/usr/local/env-activate.sh"
        set -u
    fi
    # limit meryl to use the assigned number of cores.
    export OMP_NUM_THREADS=$task.cpus

    genome_size=\$(echo "\$(<${genome_size_est})")

    /usr/local/share/merqury/best_k.sh \$genome_size $tolerable_collision > best_kmer_num.txt

    merqury.sh \\
        $meryl_db \\
        $assembly \\
        $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        merqury: $VERSION
    END_VERSIONS
    """
}
