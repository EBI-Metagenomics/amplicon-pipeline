
process concat_primers {

    label 'light'
    publishDir "${outdir}/merged/${project}", mode : "copy"
    
    input:
    tuple val(project), path(std_primers), path(auto_primers)
    val outdir

    output:
    tuple val(project), path("concat_primers.fasta") , optional: true, emit: concat_primers_out
    

    // script:
    // def auto_filter = auto_primers.name != "NO_FILE" ? "$auto_primers" : ''
    
    // TODO fix the cat
    shell:
    '''
    if [[ -s !{auto_primers} ]]; then
        cat !{std_primers} > ./concat_primers.fasta
        echo '\n' >> ./concat_primers.fasta
        cat !{auto_primers} >>./concat_primers.fasta
    else
        cat !{std_primers} > ./concat_primers.fasta
    fi
    '''
}
