
process CONCAT_PRIMERS {
    // Concatenate result of both the std primer identification and auto primer identification
    // This is important in case for a run the fwd strand has a known std primer but the reverse
    // strand has an auto primer

    label 'light'
    
    input:
    tuple val(project), val(sampleId), val(var_region), path(std_primers), path(auto_primers)
    val outdir

    output:
    tuple val(project), val(sampleId), val(var_region), path("*concat_primers.fasta") , optional: true, emit: concat_primers_out

    shell:
    '''
    if [[ -s !{auto_primers} ]]; then
        cat !{std_primers} > ./!{sampleId}_!{var_region}_concat_primers.fasta
        echo '\n' >> ./!{sampleId}_!{var_region}_concat_primers.fasta
        cat !{auto_primers} >>./!{sampleId}_!{var_region}_concat_primers.fasta
    else
        cat !{std_primers} > ./!{sampleId}_!{var_region}_concat_primers.fasta
    fi
    '''
}
