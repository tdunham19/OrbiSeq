// The only difference between each of these processes is the name of the output file

process MINIMAP2_ALIGN_TO_EXISTING {
	tag "$meta.id"

conda "${moduleDir}/environment.yml"
container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:3161f532a5ea6f1dec9be5667c9efc2afdac6104-0' :
        'biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:3161f532a5ea6f1dec9be5667c9efc2afdac6104-0' }"

input:
tuple val(meta), path(reads)
tuple val(meta2), path(reference)

output:
tuple val(meta), path("*.sam")                       , emit: sam

script:
"""
minimap2 -ax map-ont $reference $reads > ${meta.id}.sam 

cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
"""
}


process MINIMAP2_ALIGN_TO_NEW_DRAFT {
	tag "$meta.id"
	// label "no_publish"

conda "${moduleDir}/environment.yml"
container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:3161f532a5ea6f1dec9be5667c9efc2afdac6104-0' :
        'biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:3161f532a5ea6f1dec9be5667c9efc2afdac6104-0' }"

input:
tuple val(meta),  path(reads), path(input)

output:
tuple val(meta), path("*.sam")                       , emit: sam

script:
"""
minimap2 -ax map-ont $input $reads > ${meta.id}_alignment_best10.sam 

cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
"""
}

process MINIMAP2_ALIGN_TO_FINAL {
	tag "$meta.id"

conda "${moduleDir}/environment.yml"
container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:3161f532a5ea6f1dec9be5667c9efc2afdac6104-0' :
        'biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:3161f532a5ea6f1dec9be5667c9efc2afdac6104-0' }"

input:
tuple val(meta),  path(reads), path(input)

output:
tuple val(meta), path("*.sam")                       , emit: sam

script:
"""
minimap2 -ax map-ont $input $reads > ${meta.id}_alignment_new_draft_seq.sam 

cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
"""
}