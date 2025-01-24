process SED {
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*.fa"),     emit: fa

    shell:
    '''
	sed "/^>/ s/$/!{meta.id}/" "!{input}" > "!{meta.id}_new_draft_seqs.fa"
    '''
}