process CLONE_SUBWORKFLOW_REPO {
    tag "$meta.id"
    
    conda "${moduleDir}/environment.yml"

    output:
    path 'subworkflows'
    
    script:
    """
    git clone --branch ${params.subworkflowBranch} ${params.subworkflowRepo} subworkflows
    """
}