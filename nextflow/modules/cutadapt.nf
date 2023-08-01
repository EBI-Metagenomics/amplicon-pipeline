
process cutadapt {

    label 'light'
    publishDir "${outdir}/merged/${project}", mode : "copy"
    
    input:
    tuple val(project), path(std_primers), path(auto_primers),  val(sampleId), path(fastq_1), path(fastq_2)
    val outdir

    output:
    tuple val(project), path("*_1.cutadapt.fastq"), path("*_2.cutadapt.fastq"), optional: true, emit: cutadapt_out
    path "concat_primers.fasta", optional: true, emit: cutadapt_primers
    

    // script:
    // def auto_filter = auto_primers.name != "NO_FILE" ? "$auto_primers" : ''
    
    shell:
    '''
    if [[ -s !{auto_primers} ]]; then
        cat !{std_primers} !{auto_primers} > ./concat_primers.fasta
    else
        cat !{std_primers} > ./concat_primers.fasta
    fi

    fwd_primer=''
    rev_primer=''

    if [[ -s ./concat_primers.fasta ]]; then

        fwd_flag="$(grep -A 1 '^>.*F' ./concat_primers.fasta || true)"
        rev_flag="$(grep -A 1 '^>.*R' ./concat_primers.fasta || true)"

        if [[ ! -z $fwd_flag ]]; then
            grep -A 1 '^>.*F' ./concat_primers.fasta > ./fwd_primer.fasta
        else
            touch ./fwd_primer.fasta
        fi

        if [[ ! -z $rev_flag ]]; then
            grep -A 1 '^>.*R' ./concat_primers.fasta > ./rev_primer.fasta
        else
            touch ./rev_primer.fasta
        fi

        if [[ -s ./fwd_primer.fasta ]]; then
            fwd_primer="$(sed '2q;d' ./fwd_primer.fasta)"
        fi
        
        if [[ -s ./rev_primer.fasta ]]; then
            rev_primer="$(sed '2q;d' ./rev_primer.fasta)"
        fi

        if [[ ! -z $fwd_primer ]] & [[ ! -z $rev_primer ]]; then
            cutadapt -g $fwd_primer -G $rev_primer -o "!{sampleId}"_1.cutadapt.fastq -p "!{sampleId}"_2.cutadapt.fastq "!{fastq_1}" "!{fastq_2}"

        elif [[ ! -z $fwd_primer ]]; then
            cutadapt -g $fwd_primer -o "!{sampleId}"_1.cutadapt.fastq -p "!{sampleId}"_2.cutadapt.fastq "!{fastq_1}" "!{fastq_2}"
        
        elif [[ ! -z $rev_primer ]]; then
            cutadapt -G $rev_primer -o "!{sampleId}"_1.cutadapt.fastq -p "!{sampleId}"_2.cutadapt.fastq "!{fastq_1}" "!{fastq_2}"
        fi

    fi
    '''
}
