process DOWNLOAD_FROM_FIRE {

    secret 'FIRE_ACCESS_KEY'
    secret 'FIRE_SECRET_KEY'

    tag "${meta.id}"

    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/boto3:1.35.37--a82b4d378d332259' :
        'community.wave.seqera.io/library/pip_boto3:501beb4bd409b3e1' }"

    input:
    tuple val(meta), val(input_paths)

    output:
    tuple val(meta), path("downloaded_files/*fa*.gz"), emit: downloaded_files
    path "versions.yml"                              , emit: versions

    script:
    """
    s3fire_downloader.py \\
        --access-key \${FIRE_ACCESS_KEY} \\
        --secret-key \${FIRE_SECRET_KEY} \\
        --ftp-paths ${input_paths.join(" ")} \\
        --outdir downloaded_files

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        boto3: \$(python -c "import boto3; print(boto3.__version__)")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p downloaded_files
    touch downloaded_files/${meta.id}_1.fastq
    touch downloaded_files/${meta.id}_2.fastq
    gzip downloaded_files/*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        boto3: \$(python -c "import boto3; print(boto3.__version__)")
    END_VERSIONS
    """
}