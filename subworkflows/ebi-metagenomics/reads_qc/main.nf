
include { FASTP      } from '../../../modules/ebi-metagenomics/fastp/main'
include { SEQTK_SEQ     } from '../../../modules/ebi-metagenomics/seqtk/seq/main'

workflow  READS_QC {

    take:
    ch_reads    // channel: [ val(meta), [ fastq ] ]
    save_merged // channel:  val(boolean)

    main:

    ch_versions = Channel.empty()

    FASTP ( ch_reads, params.save_trimmed_fail, save_merged )
    ch_versions = ch_versions.mix(FASTP.out.versions.first())

    ch_se_fastp_reads = FASTP
                        .out.reads
                        .filter { it[0].single_end }

    ch_reads_se_and_merged = ch_se_fastp_reads
                            .mix(FASTP.out.reads_merged)

    SEQTK_SEQ(ch_reads_se_and_merged)
    ch_versions = ch_versions.mix(SEQTK_SEQ.out.versions.first())

    emit:
    reads               = FASTP.out.reads           // channel: [ val(meta), [ fastq ] ]
    reads_se_and_merged = ch_reads_se_and_merged    // channel: [ val(meta), [ fastq ] ]
    fastp_summary_json  = FASTP.out.json            // channel: [ val(meta), [ json ] ]
    reads_fasta         = SEQTK_SEQ.out.fastx       // channel: [ val(meta), [ fasta ] ]
    versions            = ch_versions               // channel: [ versions.yml ]
}

