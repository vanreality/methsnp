process REVELIO {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(input)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("${prefix}.preprocessed.bam")  , emit: bam
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
        $input \\
        ${prefix}.preprocessed.bam
    """
}
