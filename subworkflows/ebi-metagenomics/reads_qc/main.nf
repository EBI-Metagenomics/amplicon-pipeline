
include { SEQFU_CHECK            } from '../../../modules/ebi-metagenomics/seqfu/check/main'
include { ASSESSMCPPROPORTIONS   } from '../../../modules/ebi-metagenomics/assessmcpproportions/main'
include { FASTQSUFFIXHEADERCHECK } from '../../../modules/ebi-metagenomics/fastqsuffixheadercheck/main'
include { FASTP                  } from '../../../modules/ebi-metagenomics/fastp/main'
include { SEQTK_SEQ              } from '../../../modules/ebi-metagenomics/seqtk/seq/main'

workflow  READS_QC {

    take:
    filter_amplicon // channel: val(boolean)
    ch_reads    // channel: [ val(meta), [ fastq ] ]
    save_merged // channel:  val(boolean)

    main:
    ch_versions = Channel.empty()

    SEQFU_CHECK(ch_reads)
    ch_versions = ch_versions.mix(SEQFU_CHECK.out.versions.first())

    passed_seqfu_reads = SEQFU_CHECK.out.tsv
        .splitCsv(sep: "\t", elem: 1)
        .filter { meta, seqfu_res ->
            seqfu_res[0] == "OK"
        }
        .map { map, seqfu_res -> map }
        .join(ch_reads)

    FASTQSUFFIXHEADERCHECK(passed_seqfu_reads)
    ch_versions = ch_versions.mix(FASTQSUFFIXHEADERCHECK.out.versions.first())

    passed_suffixheader_reads = FASTQSUFFIXHEADERCHECK.out.json
        .filter { meta, sufhd_res ->
            sufhd_res.countLines() == 0
        }
        .map { meta, _ -> [ meta ] }
        .join(ch_reads)

    if ( filter_amplicon ) {
        assess_mcp_proportions_input = passed_suffixheader_reads
                        .map { meta, fastq ->
                            if ( meta.single_end ) {
                                [ meta, "auto", "none", fastq ]
                            } else {
                                [ meta, "auto", "none", fastq[0] ]
                            }
                        }

        ASSESSMCPPROPORTIONS(assess_mcp_proportions_input, true)
        ch_versions = ch_versions.mix(ASSESSMCPPROPORTIONS.out.versions.first())

        fastp_input = ASSESSMCPPROPORTIONS.out.library_check_out
                        .filter { meta, strategy ->
                            strategy == "AMPLICON"
                        }
                        .map { meta, _ -> [ meta ] }
                        .join(ch_reads)

        amplicon_check = ASSESSMCPPROPORTIONS.out.library_check_out
    } else {
        fastp_input = passed_suffixheader_reads
        amplicon_check = Channel.empty()
    }

    FASTP ( fastp_input, params.save_trimmed_fail, save_merged )
    ch_versions = ch_versions.mix(FASTP.out.versions.first())

    ch_se_fastp_reads = FASTP
        .out.reads
        .filter { it[0].single_end }

    ch_reads_se_and_merged = ch_se_fastp_reads
        .mix(FASTP.out.reads_merged)

    SEQTK_SEQ(ch_reads_se_and_merged)
    ch_versions = ch_versions.mix(SEQTK_SEQ.out.versions.first())

    emit:
    seqfu_check           = SEQFU_CHECK.out.tsv                        // channel: [ val(meta), tsv  ]
    suffix_header_check   = FASTQSUFFIXHEADERCHECK.out.json            // channel: [ val(meta), json  ]
    amplicon_check        = amplicon_check                             // channel: [ val(meta), env  ]
    reads                 = FASTP.out.reads                            // channel: [ val(meta), [ fastq ] ]
    reads_se_and_merged   = ch_reads_se_and_merged                     // channel: [ val(meta), [ fastq ] ]
    fastp_summary_json    = FASTP.out.json                             // channel: [ val(meta), [ json ] ]
    reads_fasta           = SEQTK_SEQ.out.fastx                        // channel: [ val(meta), [ fasta ] ]
    versions              = ch_versions                                // channel: [ versions.yml ]
}

