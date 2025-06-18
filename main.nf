include { printHeader; helpMessage } from './help'
include { WES_WF } from './workflow/wes.nf'

workflow {

    if ( params.help ) {
        helpMessage()
        System.exit(0)
    }
    
    printHeader()
    
    if ( params.mode == "exome") {
        if ( params.exome_panel == "" ) {
            error "ERROR: exome_panel is not defined.\nPlease add --exome_panel. Available options are 'xGen Exome Hyb Panel v2' aka IDT, 'TruSight One', and 'TWIST'."
        }
    }
    
    // setting ch_reads from reads or input_csv
    if ( params.reads ) {
                
        ch_reads = Channel.fromFilePairs( params.reads , size: -1 , checkExists: true )
        
        ch_reads.ifEmpty{ error "ERROR: cannot find any fastq files matching the pattern: ${params.reads}\nMake sure that the input file exists!" }

    } else if ( params.input_csv ) {
        
        if ( params.platform.equalsIgnoreCase("Ultima") ) {
            ch_reads = Channel.fromPath( params.input_csv ).splitCsv( header:true )
                            .map { row -> [ row.biosampleName, [ row.cram, row.crai ] ] }

        } else {
                
            ch_reads = Channel.fromPath( params.input_csv ).splitCsv( header:true )
                        .map { row -> [ row.biosampleName, [ row.read1, row.read2 ] ] }
                    
        }
        
        ch_reads.ifEmpty{ error "ERROR: Input csv file is empty." }
    }
    
    
    ch_reads.view()
    ch_dummy_file = Channel.fromPath("$projectDir/assets/dummy_file.txt", checkIfExists: true).collect()

    
    WES_WF( 
                params.input_csv,
                ch_reads,
                ch_dummy_file,
                params.min_reads,
                params.genomes,
                params.genome,
                params.platform,
                params.dnascope_model_selection,
                params.mode,
                params.giab_reference_name,
                params.exome_panel,
                params.variant_caller,
                params.skip_subsampling,
                params.subsample_array,
                params.read_length,
                params.seqtk_sample_seed,
                params.reference,
                params.dbsnp,
                params.dbsnp_index,
                params.mills,
                params.mills_index,
                params.onekg_omni,
                params.onekg_omni_index,
                params.ploidy,
                params.pcrfree,
                params.tmp_dir,
                params.timestamp,
                params.skip_ado,
                params.ado_infer_germ_het,
                params.ado_slice_het_sites,
                params.ado_jg_file,
                params.ado_min_samples,
                params.ado_variant_type,
                params.ado_min_af,
                params.ado_max_af,
                params.ado_sample_prop,
                params.ado_cov_cutoff,
                params.base_metrics_intervals,
                params.skip_gene_coverage,
                params.roi_intervals,
                params.refseq,
                params.skip_variant_annotation,
                params.skip_vcfeval,
                params.vcfeval_output_mode,
                params.vcfeval_score_metric,
                params.vcfeval_other_options,
                params.snpeff_db,
                params.snpeff_genome_name,
                params.hgvs_old,
                params.multiqc_config,
                params.NEAR_DISTANCE,
                params.sigprofilermatrixgenerator_reference,
                params.report_s3_dir,
                params.skip_sigprofile,
                params.publish_dir,
                params.enable_publish,
                params.disable_publish
             )
    
    
}
