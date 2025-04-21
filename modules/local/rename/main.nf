process RENAME {
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*.fa"),     emit: fa

    shell:
    '''
    awk '/^>/ {
        match($0, /best10_refseq\\.([0-9]+)/, arr);
        if (arr[1] != "") {
            print ">viral_consensus_" arr[1];
        } else {
            print $0;
        }
        next
    } { 
        print 
    }' "!{input}" > "!{meta.id}"_new_draft_seqs.fa
    '''
}