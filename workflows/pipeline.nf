
// nextflow run -profile lsf -resume main.nf --input samplesheet.csv --outdir /hps/nobackup/rdf/metagenomics/service-team/users/chrisata/asv_nf_testing/samplesheet

include { READS_QC } from '../subworkflows/ebi-metagenomics/reads_qc/main.nf'
include { READS_QC as READS_QC_MERGE } from '../subworkflows/ebi-metagenomics/reads_qc/main.nf'
include { CMSEARCH_SUBWF } from '../subworkflows/local/cmsearch_swf.nf'
include { ITS_SWF } from '../subworkflows/local/its_swf.nf'
include { MAPSEQ_OTU_KRONA as MAPSEQ_OTU_KRONA_SSU} from '../subworkflows/local/mapseq_otu_krona_swf.nf'
include { MAPSEQ_OTU_KRONA as MAPSEQ_OTU_KRONA_LSU} from '../subworkflows/local/mapseq_otu_krona_swf.nf'
include { MAPSEQ_OTU_KRONA as MAPSEQ_OTU_KRONA_PR2} from '../subworkflows/local/mapseq_otu_krona_swf.nf'
include { MAPSEQ_OTU_KRONA as MAPSEQ_OTU_KRONA_UNITE} from '../subworkflows/local/mapseq_otu_krona_swf.nf'
include { MAPSEQ_OTU_KRONA as MAPSEQ_OTU_KRONA_ITSONEDB} from '../subworkflows/local/mapseq_otu_krona_swf.nf'
include { AMP_REGION_INFERENCE } from '../subworkflows/local/amp_region_inference_swf.nf'
include { PRIMER_IDENTIFICATION } from '../subworkflows/local/primer_identification_swf.nf'
include { AUTOMATIC_PRIMER_PREDICTION } from '../subworkflows/local/automatic_primer_prediction.nf'
include { CONCAT_PRIMER_CUTADAPT } from '../subworkflows/local/concat_primer_cutadapt.nf'
include { PRIMER_VALIDATION } from '../subworkflows/local/primer_validation_swf.nf'
include { DADA2_KRONA as DADA2_KRONA_SILVA} from '../subworkflows/local/dada2_krona_swf.nf'
include { DADA2_KRONA as DADA2_KRONA_PR2} from '../subworkflows/local/dada2_krona_swf.nf'

// Initialise different database inputs for MapSeq+Krona
ssu_mapseq_krona_tuple = tuple(file(params.ssu_db_fasta), file(params.ssu_db_tax), file(params.ssu_db_otu), file(params.ssu_db_mscluster), params.ssu_label)
lsu_mapseq_krona_tuple = tuple(file(params.lsu_db_fasta), file(params.lsu_db_tax), file(params.lsu_db_otu), file(params.lsu_db_mscluster), params.lsu_label)
itsonedb_mapseq_krona_tuple = tuple(file(params.itsone_db_fasta), file(params.itsone_db_tax), file(params.itsone_db_otu), file(params.itsone_db_mscluster), params.itsone_label)
unite_mapseq_krona_tuple = tuple(file(params.unite_db_fasta), file(params.unite_db_tax), file(params.unite_db_otu), file(params.unite_db_mscluster), params.unite_label)
pr2_mapseq_krona_tuple = tuple(file(params.pr2_db_fasta), file(params.pr2_db_tax), file(params.pr2_db_otu), file(params.pr2_db_mscluster), params.pr2_label)

// Initialise database inputs for DADA2+Krona
silva_dada2_db = file(params.silva_dada2_db)
dada2_krona_silva_tuple = tuple(file(params.ssu_db_fasta), file(params.ssu_db_tax), file(params.ssu_db_otu), file(params.ssu_db_mscluster), params.dada2_silva_label)

pr2_dada2_db = file(params.pr2_dada2_db)
dada2_krona_pr2_tuple = tuple(file(params.ssu_db_fasta), file(params.ssu_db_tax), file(params.ssu_db_otu), file(params.ssu_db_mscluster), params.dada2_pr2_label)

// Read input samplesheet
samplesheet = Channel.fromSamplesheet( "input" )

workflow AMPLICON_PIPELINE_V6 {

    // Organise input tuple channel
    groupReads = { meta, fq1, fq2 ->
        if (fq2 == []) {
            return tuple(meta, [fq1])
        }
        else {
            return tuple(meta, [fq1, fq2])
        }
    }

    ch_input = samplesheet
                  .map(groupReads)

    // Quality control
    READS_QC_MERGE(
        ch_input,
        true
    )
    // Run it again without merging to keep PE files unmerged for primer trimming+DADA2
    READS_QC(
        ch_input,
        false
    )

    // Cmsearch subworkflow to find rRNA reads for SSU+LSU
    CMSEARCH_SUBWF(
        READS_QC_MERGE.out.reads_fasta,
        file(params.rfam),
        file(params.claninfo)
    )

    // Masking subworkflow to find rRNA reads for ITS
    ITS_SWF(
        READS_QC_MERGE.out.reads_fasta,
        CMSEARCH_SUBWF.out.concat_ssu_lsu_coords
    )

    // Next five subworkflow calls are MapSeq annotation + Krona generation for SSU+LSU+ITS
    MAPSEQ_OTU_KRONA_SSU(
        CMSEARCH_SUBWF.out.ssu_fasta,
        ssu_mapseq_krona_tuple
    )

    MAPSEQ_OTU_KRONA_PR2(
        CMSEARCH_SUBWF.out.ssu_fasta,
        pr2_mapseq_krona_tuple
    )  

    MAPSEQ_OTU_KRONA_LSU(
        CMSEARCH_SUBWF.out.lsu_fasta,
        lsu_mapseq_krona_tuple
    )     

    MAPSEQ_OTU_KRONA_ITSONEDB(
        ITS_SWF.out.its_masked_out,
        itsonedb_mapseq_krona_tuple
    )    

    MAPSEQ_OTU_KRONA_UNITE(
        ITS_SWF.out.its_masked_out,
        unite_mapseq_krona_tuple
    )

    // Infer amplified variable regions for SSU, extract reads for each amplified region if there are more than one
    AMP_REGION_INFERENCE(
        CMSEARCH_SUBWF.out.cmsearch_deoverlap_out,
        READS_QC_MERGE.out.reads_se_and_merged
    )

    // Identify whether primers exist or not in reads, separated by different amplified regions if more than one exists in a run
    PRIMER_IDENTIFICATION(
        AMP_REGION_INFERENCE.out.extracted_var_out
    )

    // Join primer identification flags with reads belonging to each run+amp_region
    auto_trimming_input = PRIMER_IDENTIFICATION.out.conductor_out
                          .join(AMP_REGION_INFERENCE.out.extracted_var_out, by: [0])

    // Run subworkflow for automatic primer prediction
    // Outputs empty fasta file if no primers, or fasta file containing predicted primers
    AUTOMATIC_PRIMER_PREDICTION(
        auto_trimming_input
    )

    // Concatenate the different combinations of stranded std/auto primers for each run+amp_region
    concat_input = PRIMER_IDENTIFICATION.out.std_primer_out
                   .join(AUTOMATIC_PRIMER_PREDICTION.out.auto_primer_trimming_out, by: [0])

    // Concatenate all primers for for a run, send them to cutadapt with original QCd reads for primer trimming
    CONCAT_PRIMER_CUTADAPT(
        concat_input,
        READS_QC.out.reads
    )

    primer_validation_input = CONCAT_PRIMER_CUTADAPT.out.final_concat_primers_out
                              .map{ tuple(it[0], it[1]) }

    // Verify that any identified primers (both std+auto) actually match to regions of the SSU gene (for Bacteria/Archaea/Eukaryotes)
    // Output of this (a .tsv file) will go to CDCH
    // TODO THIS SUBWORKFLOW NEEDS REFACTORING
    PRIMER_VALIDATION(
        primer_validation_input
    )

    cutadapt_channel = CONCAT_PRIMER_CUTADAPT.out.cutadapt_out
                       .map{ meta, reads -> 
                         [ meta.subMap('id', 'single_end'), meta['var_region'], meta['var_regions_size'], reads ]
                        }

    dada2_input = concat_input
                  .map{ meta, std_primer, auto_primer -> 
                    key = groupKey(meta.subMap('id', 'single_end'), meta['var_regions_size'])
                    [ key, meta['var_region'], meta['var_regions_size'] ]
                   }
                  .groupTuple(by: 0)
                  .join(READS_QC.out.reads, by: 0)
                  .mix(cutadapt_channel)
                  .map { meta, var_region, var_regions_size, reads -> 
                    [ groupKey(meta.subMap('id', 'single_end'), 2), var_region, var_regions_size, reads ]
                  }
                  .groupTuple(by:0) 
                  .map{ meta, var_region, var_regions_size, reads -> 
                    final_var_region = var_region.unique()[0]
                    final_var_region_size = var_regions_size[0][0]
                    final_meta = meta + ['var_region': final_var_region, 'var_regions_size': final_var_region_size]

                    fastp_reads = reads[0]
                    cutadapt_reads = reads[1]
                    final_reads = ''
                    if (meta['single_end'] ){
                        if (cutadapt_reads.size() > 0){
                            final_reads = cutadapt_reads
                        }
                        else{
                            final_reads = fastp_reads
                        }
                    }
                    else{
                        if (cutadapt_reads[0].size() > 0){
                            final_reads = cutadapt_reads
                        }
                        else{
                            final_reads = fastp_reads
                        }
                    }
                    [ final_meta, final_reads ]
                  }


    // Run DADA2 ASV generation + generate Krona plots for each run+amp_region 
    DADA2_KRONA_SILVA(
        dada2_input,
        AMP_REGION_INFERENCE.out.concat_var_regions,
        AMP_REGION_INFERENCE.out.extracted_var_path,
        READS_QC.out.reads,
        silva_dada2_db,
        dada2_krona_silva_tuple,
    )

    DADA2_KRONA_PR2(
        dada2_input,
        AMP_REGION_INFERENCE.out.concat_var_regions,
        AMP_REGION_INFERENCE.out.extracted_var_path,
        READS_QC.out.reads,
        pr2_dada2_db,
        dada2_krona_pr2_tuple,
    )

}