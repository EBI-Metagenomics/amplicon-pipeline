
include { READS_QC                                      } from '../subworkflows/ebi-metagenomics/reads_qc/main.nf'
include { READS_QC as READS_QC_MERGE                    } from '../subworkflows/ebi-metagenomics/reads_qc/main.nf'
include { RRNA_EXTRACTION                               } from '../subworkflows/ebi-metagenomics/rrna_extraction/main'
include { MASK_FASTA_SWF                                } from '../subworkflows/local/mask_fasta_swf.nf'
include { MAPSEQ_OTU_KRONA as MAPSEQ_OTU_KRONA_SSU      } from '../subworkflows/ebi-metagenomics/mapseq_otu_krona/main'
include { MAPSEQ_OTU_KRONA as MAPSEQ_OTU_KRONA_LSU      } from '../subworkflows/ebi-metagenomics/mapseq_otu_krona/main'
include { MAPSEQ_OTU_KRONA as MAPSEQ_OTU_KRONA_PR2      } from '../subworkflows/ebi-metagenomics/mapseq_otu_krona/main'
include { MAPSEQ_OTU_KRONA as MAPSEQ_OTU_KRONA_UNITE    } from '../subworkflows/ebi-metagenomics/mapseq_otu_krona/main'
include { MAPSEQ_OTU_KRONA as MAPSEQ_OTU_KRONA_ITSONEDB } from '../subworkflows/ebi-metagenomics/mapseq_otu_krona/main'

include { AMP_REGION_INFERENCE                          } from '../subworkflows/local/amp_region_inference_swf.nf'
include { PRIMER_IDENTIFICATION                         } from '../subworkflows/local/primer_identification_swf.nf'
include { AUTOMATIC_PRIMER_PREDICTION                   } from '../subworkflows/local/automatic_primer_prediction.nf'
include { CONCAT_PRIMER_CUTADAPT                        } from '../subworkflows/local/concat_primer_cutadapt.nf'
include { PRIMER_VALIDATION                             } from '../subworkflows/local/primer_validation_swf.nf'
include { DADA2_SWF                                     } from '../subworkflows/local/dada2_swf.nf'
include { MAPSEQ_ASV_KRONA as MAPSEQ_ASV_KRONA_SILVA    } from '../subworkflows/local/mapseq_asv_krona_swf.nf'
include { MAPSEQ_ASV_KRONA as MAPSEQ_ASV_KRONA_PR2      } from '../subworkflows/local/mapseq_asv_krona_swf.nf'

include { dada2_input_preparation_function              } from '../lib/nf/dada2_input_preparation_function.nf'

include { samplesheetToList } from 'plugin/nf-schema'


// Initialise different database inputs for taxonomic assignments with regular taxonomy resolution method
ssu_mapseq_krona_tuple = Channel.value([
    file(params.ssu_db_fasta, checkIfExists: true),
    file(params.ssu_db_tax, checkIfExists: true),
    file(params.ssu_db_otu, checkIfExists: true),
    file(params.ssu_db_mscluster, checkIfExists: true),
    params.ssu_label
])
lsu_mapseq_krona_tuple = Channel.value([
    file(params.lsu_db_fasta, checkIfExists: true),
    file(params.lsu_db_tax, checkIfExists: true),
    file(params.lsu_db_otu, checkIfExists: true),
    file(params.lsu_db_mscluster, checkIfExists: true),
    params.lsu_label
])
itsonedb_mapseq_krona_tuple = Channel.value([
    file(params.itsone_db_fasta, checkIfExists: true),
    file(params.itsone_db_tax, checkIfExists: true),
    file(params.itsone_db_otu, checkIfExists: true),
    file(params.itsone_db_mscluster, checkIfExists: true),
    params.itsone_label
])
unite_mapseq_krona_tuple = Channel.value([
    file(params.unite_db_fasta, checkIfExists: true),
    file(params.unite_db_tax, checkIfExists: true),
    file(params.unite_db_otu, checkIfExists: true),
    file(params.unite_db_mscluster, checkIfExists: true),
    params.unite_label
])
pr2_mapseq_krona_tuple = Channel.value([
    file(params.pr2_db_fasta, checkIfExists: true),
    file(params.pr2_db_tax, checkIfExists: true),
    file(params.pr2_db_otu, checkIfExists: true),
    file(params.pr2_db_mscluster, checkIfExists: true),
    params.pr2_label
])

// Initialise different database inputs for taxonomic assignments with ASV resolution method
dada2_krona_silva_tuple = tuple(
    file(params.ssu_db_fasta, checkIfExists: true),
    file(params.ssu_db_tax, checkIfExists: true),
    file(params.ssu_db_otu, checkIfExists: true),
    file(params.ssu_db_mscluster, checkIfExists: true),
    params.dada2_silva_label
)
dada2_krona_pr2_tuple = tuple(
    file(params.pr2_db_fasta, checkIfExists: true),
    file(params.pr2_db_tax, checkIfExists: true),
    file(params.pr2_db_otu, checkIfExists: true),
    file(params.pr2_db_mscluster, checkIfExists: true),
    params.dada2_pr2_label
)

// Standard primer library
std_primer_library = file(params.std_primer_library, type: 'dir', checkIfExists: true)

// Read input samplesheet
// samplesheet = Channel.fromSamplesheet( "input" )

samplesheet = Channel.fromList(samplesheetToList(params.input, "./assets/schema_input.json"))

workflow AMPLICON_PIPELINE {

    // Organise input tuple channel
    groupReads = { meta, fq1, fq2 ->
        if (fq2 == []) {
            return tuple(meta, [fq1])
        }
        else {
            return tuple(meta, [fq1, fq2])
        }
    }

    ch_input = samplesheet.map(groupReads)

    // Quality control
    READS_QC_MERGE(
        ch_input,
        true // merge
    )
    // Run it again without merging to keep PE files unmerged for primer trimming+DADA2
    READS_QC(
        ch_input,
        false // merge
    )

    // rRNA extraction subworkflow to find rRNA reads for SSU+LSU
    RRNA_EXTRACTION(
        READS_QC_MERGE.out.reads_fasta,
        file( params.rfam, checkIfExists: true ),
        file( params.claninfo, checkIfExists: true )
    )

    // Masking subworkflow to find rRNA reads for ITS
    MASK_FASTA_SWF(
        READS_QC_MERGE.out.reads_fasta,
        RRNA_EXTRACTION.out.concat_ssu_lsu_coords
    )

    // Next five subworkflow calls are MapSeq annotation + Krona generation for SSU+LSU+ITS
    MAPSEQ_OTU_KRONA_SSU(
        RRNA_EXTRACTION.out.ssu_fasta,
        ssu_mapseq_krona_tuple
    )

    MAPSEQ_OTU_KRONA_PR2(
        RRNA_EXTRACTION.out.ssu_fasta,
        pr2_mapseq_krona_tuple
    )  

    MAPSEQ_OTU_KRONA_LSU(
        RRNA_EXTRACTION.out.lsu_fasta,
        lsu_mapseq_krona_tuple
    )     

    MAPSEQ_OTU_KRONA_ITSONEDB(
        MASK_FASTA_SWF.out.masked_out,
        itsonedb_mapseq_krona_tuple
    )    

    MAPSEQ_OTU_KRONA_UNITE(
        MASK_FASTA_SWF.out.masked_out,
        unite_mapseq_krona_tuple
    )

    // Infer amplified variable regions for SSU, extract reads for each amplified region if there are more than one
    AMP_REGION_INFERENCE(
        RRNA_EXTRACTION.out.cmsearch_deoverlap_out,
        READS_QC_MERGE.out.reads_se_and_merged
    )

    // Identify whether primers exist or not in reads, separated by different amplified regions if more than one exists in a run
    PRIMER_IDENTIFICATION(
        AMP_REGION_INFERENCE.out.extracted_var_out,
        std_primer_library
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
                              .map{ meta, primers ->
                                if (primers.size() > 0){
                                    [ meta, primers ]
                                }    
                              }

    // Verify that any identified primers (both std+auto) actually match to regions of the SSU gene (for Bacteria/Archaea/Eukaryotes)
    // Output of this (a .tsv file) will go to CDCH
    // TODO THIS SUBWORKFLOW NEEDS REFACTORING
    PRIMER_VALIDATION(
        primer_validation_input
    )

    cutadapt_channel = CONCAT_PRIMER_CUTADAPT.out.cutadapt_out
                       .map { meta, reads -> 
                         [ meta.subMap('id', 'single_end'), meta['var_region'], meta['var_regions_size'], reads ]
                       }

    dada2_input = dada2_input_preparation_function(concat_input, READS_QC.out.reads, cutadapt_channel)

    // Run DADA2 ASV generation
    DADA2_SWF(
        dada2_input
    )

    // ASV taxonomic assignments + generate Krona plots for each run+amp_region
    MAPSEQ_ASV_KRONA_SILVA(
        DADA2_SWF.out.dada2_out,
        AMP_REGION_INFERENCE.out.concat_var_regions,
        AMP_REGION_INFERENCE.out.extracted_var_path,
        dada2_krona_silva_tuple,
    )

    MAPSEQ_ASV_KRONA_PR2(
        DADA2_SWF.out.dada2_out,
        AMP_REGION_INFERENCE.out.concat_var_regions,
        AMP_REGION_INFERENCE.out.extracted_var_path,
        dada2_krona_pr2_tuple,
    )

}