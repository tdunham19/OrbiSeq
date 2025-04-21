process RENAME {
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*.fa"),     emit: fa

    shell:
     '''
     awk 'NR==1 {
        match($0, /_refseq\\.([0-9]+)/, arr);
        if (arr[1] != "") {
            print ">viral_consensus_best10_refseq." arr[1];
        } else {
            print $0;
        }
        next;
    } { 
        print;
    }' "!{input}" > temp.fa && mv temp.fa "!{input}"
    '''
}