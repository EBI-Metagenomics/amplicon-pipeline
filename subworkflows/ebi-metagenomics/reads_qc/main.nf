
include { SEQFU_CHECK            } from '../../../modules/ebi-metagenomics/seqfu/check/main'
include { PIMENTO_GENERATEBCV    } from '../../../modules/ebi-metagenomics/pimento/generatebcv/main'
include { LIBRARYSTRATEGYCHECK    } from '../../../modules/ebi-metagenomics/librarystrategycheck/main'
include { FASTQSUFFIXHEADERCHECK } from '../../../modules/ebi-metagenomics/fastqsuffixheadercheck/main'
include { BBMAP_REFORMAT_STANDARDISE } from '../../../modules/local/bbmap/reformat_standardise/main'
include { FASTP                  } from '../../../modules/ebi-metagenomics/fastp/main'
include { SEQTK_SEQ              } from '../../../modules/ebi-metagenomics/seqtk/seq/main'

workflow  READS_QC {

    take:
    filter_amplicon // channel: val(boolean)
    ch_reads    // channel: [ val(meta), [ fastq ] ]
    save_merged // channel:  val(boolean)

    main:
    ch_versions = Channel.empty()

    BBMAP_REFORMAT_STANDARDISE(ch_reads, 'fastq.gz')

    SEQFU_CHECK(BBMAP_REFORMAT_STANDARDISE.out.reformated)
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
        generatebcv_input = passed_suffixheader_reads
                        .map { meta, fastq ->
                            if ( meta.single_end ) {
                                [ meta, "auto", "auto", fastq ]
                            } else {
                                [ meta, "auto", "auto", fastq[0] ]
                            }
                        }

        PIMENTO_GENERATEBCV(generatebcv_input)
        ch_versions = ch_versions.mix(PIMENTO_GENERATEBCV.out.versions.first())

        LIBRARYSTRATEGYCHECK(PIMENTO_GENERATEBCV.out.tsv)
        ch_versions = ch_versions.mix(LIBRARYSTRATEGYCHECK.out.versions.first())

        fastp_input = LIBRARYSTRATEGYCHECK.out.library_check_out
                        .filter { meta, strategy ->
                            strategy == "AMPLICON"
                        }
                        .map { meta, _ -> [ meta ] }
                        .join(ch_reads)

        amplicon_check = LIBRARYSTRATEGYCHECK.out.library_check_out
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

