#!/usr/bin/env nextflow

// If the user uses the --help flag, print the help text below
params.help = false

// Function which prints help message text
def helpMessage() {
  log.info"""
    QC for MAGs

    Required Arguments:
      --manifest				Location where 1 or more fasta files are stored.
      --project 	            Folder to place analysis outputs (default: ${params.project})

    Options
      --outdir              Base directory for output files (default: ${params.outdir})
    """.stripIndent()
}

log.info"""Starting""".stripIndent()

// Show help message if the user specifies the --help flag at runtime
if (params.help) {
  // Invoke the function above which prints the help message
  helpMessage()
  // Exit out and do not run anything else
  exit 0
}

// // Show help message if the user specifies a fasta file but not makedb or db
if (params.manifest == "") {
  // Invoke the function above which prints the help message
  helpMessage()
  // Exit out and do not run anything else
  exit 0
}

def outputBase = "${params.outdir}/${params.project}"

/*
 * Defines the pipeline inputs parameters (giving a default value for each for them)
 * Each of the following parameters can be specified as command line options
 */
Channel
	.fromPath(params.manifest)
	.ifEmpty { exit 1, "Cannot find any seed file matching: ${params.manifest}." }
	.splitCsv(header:true)
	.map{ row -> tuple(row.sample_name, file(row.bam_file), file(row.bai_file), file(row.binmap_file), row.library_size, row.target_genome)}
	.set { manifest_ch }

process ninjamap2 {
	tag "${prefix}"

	container params.docker_container_ninjamap

	publishDir "$outputBase/${prefix}/", pattern: "*.{bedgraph,json,log}"
	// "$outputBase/${prefix}/Logs", mode: 'copy', pattern: "*.{log}"
	// "$outputBase/${prefix}/Bedgraph", mode: 'copy', pattern: "*.{bedgraph}"
	// "$outputBase/${prefix}/Stats", mode: 'copy', pattern: "*.{tsv,json}"

	input:
	tuple val(prefix), path(bam_file), path(bai_file), path(binmap_file), val(library_size), val(target_genome) from manifest_ch

	output:
	path "${prefix}*json"
	path "${prefix}*bedgraph" optional true
	path ".command.log" optional true

	script:
	"""
	ls -lhtra
	ninjamap2.py $bam_file $binmap_file $library_size $prefix $target_genome
	ls -lhtra
	"""

}