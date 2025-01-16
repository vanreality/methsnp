process REVELIO {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path("*.preprocessed.bam")  , emit: bam

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
