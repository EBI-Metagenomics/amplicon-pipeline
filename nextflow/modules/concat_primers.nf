
process CONCAT_PRIMERS {
    // Concatenate result of both the std primer identification and auto primer identification
    // This is important in case for a run the fwd strand has a known std primer but the reverse
    // strand has an auto primer

    label 'light'
    publishDir "${outdir}/${project}", mode : "copy"
    
    input:
    tuple val(project), val(sampleId), val(var_region), path(std_primers), path(auto_primers)
    val outdir

    output:
    tuple val(project), val(sampleId), val(var_region), path("*concat_primers.fasta") , optional: true, emit: concat_primers_out
    

    // script:
    // def auto_filter = auto_primers.name != "NO_FILE" ? "$auto_primers" : ''
    
    // TODO fix the cat
    shell:
    '''
    if [[ -s !{auto_primers} ]]; then
        cat !{std_primers} > ./!{var_region}_concat_primers.fasta
        echo '\n' >> ./!{var_region}_concat_primers.fasta
        cat !{auto_primers} >>./!{var_region}_concat_primers.fasta
    else
        cat !{std_primers} > ./!{var_region}_concat_primers.fasta
    fi
    '''
}
