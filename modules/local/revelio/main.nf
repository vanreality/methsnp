process REVELIO {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(bai)
    tuple val(meta3), path(fasta)
    tuple val(meta4), path(fai)

    output:
    tuple val(meta), path("*.preprocessed.bam")  , emit: bam
    path  "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python3 ${projectDir}/scripts/revelio/revelio.py \\
        -f ${fasta} \\
        -T ${task.cpus-1} \\
        -t . \\
        $args \\
        $bam \\
        ${prefix}.preprocessed.bam
    """
}
