
process CUTADAPT {
    // Trim the given primers using cutadapt for two paired-end read files

    label 'light'
    publishDir "${outdir}/${project}", mode : "copy"
    
    input:
    tuple val(project), val(sampleId), path(concat_primers), path(fastq_1), path(fastq_2)
    val outdir

    output:
    tuple val(project), val(sampleId), path("*_1.cutadapt.fastq.gz"), path("*_2.cutadapt.fastq.gz"), optional: true, emit: cutadapt_out
        
    shell:
    '''
    fwd_primer=''
    rev_primer=''

    if [[ -s !{concat_primers} ]]; then

        fwd_flag="$(grep -A 1 '^>.*F' !{concat_primers} || true)"
        rev_flag="$(grep -A 1 '^>.*R' !{concat_primers} || true)"

        if [[ ! -z $fwd_flag ]]; then
            grep -A 1 '^>.*F' !{concat_primers} > ./fwd_primer.fasta
        else
            touch ./fwd_primer.fasta
        fi

        if [[ ! -z $rev_flag ]]; then
            grep -A 1 '^>.*R' !{concat_primers} > ./rev_primer.fasta
        else
            touch ./rev_primer.fasta
        fi

        if [[ -s ./fwd_primer.fasta ]]; then
            fwd_primer="$(sed '2q;d' ./fwd_primer.fasta)"
        fi
        
        if [[ -s ./rev_primer.fasta ]]; then
            rev_primer="$(sed '2q;d' ./rev_primer.fasta)"
        fi

        if [[ ! -z $fwd_primer ]] && [[ ! -z $rev_primer ]]; then
            cutadapt -g $fwd_primer -G $rev_primer -o "!{sampleId}"_1.cutadapt.fastq.gz -p "!{sampleId}"_2.cutadapt.fastq.gz "!{fastq_1}" "!{fastq_2}"

        elif [[ ! -z $fwd_primer ]]; then
            cutadapt -g $fwd_primer -o "!{sampleId}"_1.cutadapt.fastq.gz -p "!{sampleId}"_2.cutadapt.fastq.gz "!{fastq_1}" "!{fastq_2}"
        
        elif [[ ! -z $rev_primer ]]; then
            cutadapt -G $rev_primer -o "!{sampleId}"_1.cutadapt.fastq.gz -p "!{sampleId}"_2.cutadapt.fastq.gz "!{fastq_1}" "!{fastq_2}"
        fi

    fi
    '''
}
