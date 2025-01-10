nextflow.enable.dsl=2

include { printHeader; helpMessage } from './help' params ( params )
include { WES_WF } from './workflow/wes.nf' params ( params )

if ( params.help ) {
    helpMessage()
    exit 0
}


workflow {
    
    printHeader()
    
    if ( params.mode == "exome") {
        if ( params.exome_panel == "" ) {
            exit 1, "ERROR: exome_panel is not defined.\nPlease add --exome_panel. Available options are 'xGen Exome Hyb Panel v2' aka IDT, 'TruSight One', and 'TWIST'."
        }
    }
    
    // setting ch_reads from reads or input_csv
    if ( params.reads ) {
                
        ch_reads = Channel.fromFilePairs( params.reads , size: -1 , checkExists: true )
        
        ch_reads.ifEmpty{ exit 1, "ERROR: cannot find any fastq files matching the pattern: ${params.reads}\nMake sure that the input file exists!" }

    } else if ( params.input_csv ) {
        
        if ( params.platform.equalsIgnoreCase("Ultima") ) {
            ch_reads = Channel.fromPath( params.input_csv ).splitCsv( header:true )
                            .map { row -> [ row.biosampleName, [ row.cram, row.crai ] ] }

        } else {
                
            ch_reads = Channel.fromPath( params.input_csv ).splitCsv( header:true )
                        .map { row -> [ row.biosampleName, [ row.read1, row.read2 ] ] }
                    
        }
        
        ch_reads.ifEmpty{ exit 1, "ERROR: Input csv file is empty." }
    }
    
    
    ch_reads.view()
    ch_dummy_file = Channel.fromPath("$projectDir/assets/dummy_file.txt", checkIfExists: true).collect()
    // ch_dummy_file2 = Channel.fromPath("$projectDir/assets/dummy_file2.txt", checkIfExists: true).collect()

    
    WES_WF( 
                params.input_csv,
                ch_reads,
                ch_dummy_file,
                params.min_reads
             )
    
    
}

// OnComplete
workflow.onComplete{
    println( "\nPipeline completed successfully.\n\n" )
}

// OnError
workflow.onError{
    println( "\nPipeline failed.\n\n" )
}
