
process EXTRACT_PROJECT_NAME {

    input:
    val proj_path

    output:
    val project, emit: projects 

    script:
    project = "${proj_path}".split('/').last()

    """
    echo $project
    """
}