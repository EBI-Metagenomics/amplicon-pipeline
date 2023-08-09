
include { std_primer_flag } from '../modules/std_primer_flag.nf'
include { general_primer_flag } from '../modules/general_primer_flag.nf'
include { trimming_conductor } from '../modules/trimming_conductor.nf'
include { parse_conductor } from '../modules/parse_conductor.nf'

workflow PRIMER_IDENTIFICATION {
    
    // Subworkflow that attempts to identify the presence of primers in sequencing files (fastq)
    // Checks for the presence of standard primers from a known library in ""./data/standard_primers"
    // Also checks for the presence of primers in general in case a presence is present but is not part of the current library
    // Outputs strand-based flags outlining whether primers were found, and if so, how they should be removed 
    // This can be either removal through matches to known standard primers, or the need to predict the used primer if it's not part of the current library

    take:
        merged_reads
        outdir

    main:
        // Check for presence of primers
        std_primer_flag(merged_reads, outdir) // Standard Library primers
        general_primer_flag(merged_reads, outdir) // Primers in general

        // Combining std and general primer outputs and parsing them to guide 
        // samples through automatic primer identification and trimming by cutadapt
        comb_flags = general_primer_flag.out.general_primer_out
        .join(std_primer_flag.out.std_primer_out)
        trimming_conductor(comb_flags, outdir) // Generate a flags file with vals of 'none/std/auto' for both fwd and rev
        parse_conductor(trimming_conductor.out.trimming_conductor_out, outdir) // Parse flags file into env variables

    
    emit:
        conductor_out = parse_conductor.out.conductor_out
        std_primer_out = std_primer_flag.out.std_primer_out
    
}