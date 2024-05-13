
include { INFERNAL_CMSEARCH           } from '../../../modules/ebi-metagenomics/infernal/cmsearch/main'
include { CMSEARCHTBLOUTDEOVERLAP     } from '../../../modules/ebi-metagenomics/cmsearchtbloutdeoverlap/main'
include { EASEL_ESLSFETCH             } from '../../../modules/ebi-metagenomics/easel/eslsfetch/main'
include { EXTRACTCOORDS               } from '../../../modules/ebi-metagenomics/extractcoords/main'

workflow RRNA_EXTRACTION {

    take:
    ch_fasta     // channel: [ val(meta), [ fasta ] ]
    rfam         // file: rfam for cmsearch
    claninfo     // file: claninfo for cmsearchtbloutdeoverlap

    main:

    ch_versions = Channel.empty()

    INFERNAL_CMSEARCH(
        ch_fasta,
        rfam
    )
    ch_versions = ch_versions.mix(INFERNAL_CMSEARCH.out.versions.first())

    CMSEARCHTBLOUTDEOVERLAP(
        INFERNAL_CMSEARCH.out.cmsearch_tbl,
        claninfo
    )
    ch_versions = ch_versions.mix(CMSEARCHTBLOUTDEOVERLAP.out.versions.first())

    ch_easel = ch_fasta
                .join(CMSEARCHTBLOUTDEOVERLAP.out.cmsearch_tblout_deoverlapped)

    EASEL_ESLSFETCH(
        ch_easel
    )
    ch_versions = ch_versions.mix(EASEL_ESLSFETCH.out.versions.first())

    EXTRACTCOORDS(
        EASEL_ESLSFETCH.out.easel_coords,
        EASEL_ESLSFETCH.out.matched_seqs_with_coords
    )
    ch_versions = ch_versions.mix(EXTRACTCOORDS.out.versions.first())

    emit:
    cmsearch_deoverlap_out = CMSEARCHTBLOUTDEOVERLAP.out.cmsearch_tblout_deoverlapped   // channel: [ val(meta), [ deoverlapped ] ]
    easel_out = EASEL_ESLSFETCH.out.easel_coords                                        // channel: [ val(meta), [ fasta ] ]
    ssu_fasta = EXTRACTCOORDS.out.ssu_fasta                                             // channel: [ val(meta), [ fasta ] ]
    lsu_fasta = EXTRACTCOORDS.out.lsu_fasta                                             // channel: [ val(meta), [ fasta ] ]
    concat_ssu_lsu_coords = EXTRACTCOORDS.out.concat_ssu_lsu_coords                     // channel: [ val(meta), [ txt ] ]
    versions = ch_versions                                                              // channel: [ versions.yml ]
}

