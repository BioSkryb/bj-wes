nextflow.enable.dsl=2
params.timestamp = ""

process PICARD_COLLECTHSMETRICS {
    tag "${sample_name}"
    publishDir "${publish_dir}_${params.timestamp}/", enabled:"$enable_publish"
    
    
    input:
    tuple val(sample_name), path(bam), path(bai)
    path reference
    path intervals
    val(publish_dir)
    val(enable_publish)
    
    output:
    tuple val(sample_name), file("*sentieonmetrics*"), emit: metrics_tuple
    path("picard_version.yml"), emit: version
    
    script:
    def avail_mem = 3
    if(task.memory){
        avail_mem = task.memory.giga
    }
    def near_distance_param =     params.NEAR_DISTANCE ? "--NEAR_DISTANCE ${params.NEAR_DISTANCE}" : ''
    
    """
    picard \\
            -Xmx${avail_mem}g \\
            CollectHsMetrics \\
            --INPUT ${bam} \\
            --OUTPUT ${sample_name}.hsmetricalgo.sentieonmetrics.txt \\
            --BAIT_INTERVALS ${intervals} \\
            --TARGET_INTERVALS ${intervals} \\
            --REFERENCE_SEQUENCE ${reference}/genome.fa \\
            ${near_distance_param}
        
    export PICARD_VER=\$(echo \$(picard CollectHsMetrics --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
    echo Picard: \$PICARD_VER > picard_version.yml
    """
}

workflow PICARD_COLLECTHSMETRICS_WF {
    take:
        ch_bam
        ch_reference
        ch_intervals
        ch_publish_dir
        ch_enable_publish
    main:
        PICARD_COLLECTHSMETRICS (
                                    ch_bam,
                                    ch_reference,
                                    ch_intervals,
                                    ch_publish_dir,
                                    ch_enable_publish
                                )
    emit:
        report = PICARD_COLLECTHSMETRICS.out.metrics
        version = PICARD_COLLECTHSMETRICS.out.version
}

workflow{
    

    log.info """\
      reference                       : ${ params.reference }
      calling_intervals_filename      : ${ params.calling_intervals_filename }
      \n
     """

    
    ch_bam_raw = Channel.fromFilePairs(params.bam, size: -1)
    ch_bam_raw
              .map{ it -> it.flatten().collect() }
              .set{ ch_bam }
              
    // ch_bam.view()
    
    
    
    PICARD_COLLECTHSMETRICS_WF (
                                    ch_bam,
                                    params.reference,
                                    params.calling_intervals_filename,
                                    params.publish_dir,
                                    params.enable_publish
                                )
}

