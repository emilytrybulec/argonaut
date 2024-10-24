process KMER_FREQ {
    tag "$meta.id"
    label 'process_high'

    container 'emilytrybulec/genassembly:kmer'

    input:
    tuple val(meta), path(ontfile)

    output:
    tuple val(meta), path("*.2colum")              , emit: kmerstat
    tuple val(meta), path("*.txt")                 , emit: kmernum
    path  "versions.yml"                           , emit: versions

    script:
    def VERSION = '4.0'

    """
    ls "$ontfile" > read_files.lib 
    /kmerfreq/kmerfreq -t 10 read_files.lib

    #calculate kmer number
    less read_files.lib.kmer.freq.stat | grep "#Kmer indivdual number" > kmernum.txt
    less read_files.lib.kmer.freq.stat | perl -ne 'next if(/^#/ || /^\s/); print; ' | awk '{print \$1"\t"\$2}' > read_files.lib.kmer.freq.stat.2colum 
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kmerfreq: $VERSION
    END_VERSIONS
    """
}