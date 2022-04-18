#!/usr/bin/env nextflow

/*
 Gene-centered analysis of paired-end and single-end RNA-seq data
 Authors: 
	- Yuvia A. PEREZ RICO <yuvia.perez-rico@embl.de>
	- Lena CLERQUIN lena.clerquin@embl.de>
*/

log.info "            standard RNA-seq - version 0.0.1           "
log.info "#######################################################"
log.info "Sample description file	= ${params.sampleInfo}"
log.info "Reads per split file		= ${params.chunkSize}"
log.info "Output directory		= ${params.outDir}"
log.info "Path to genome index		= ${params.genomeDirPath}"
log.info "Gene annotations		= ${params.annotations}"
log.info "Script to calculate TPM	= ${params.rsTPM}"
log.info "Paired-end data		= ${params.pairedEnd}"
log.info "Number of threads		= ${params.numCPUs}"
log.info "\n"

/*
 Validate input parameters
*/

if( !(params.chunkSize instanceof Number) ){
	exit 1, "Invalid chunk size = ${params.chunkSize}"
}

if( !(params.numCPUs instanceof Number) ){
	exit 1, "Invalid number of CPUs = ${params.numCPUs}"
}

/*
 Validate input files
*/

sdFile = file(params.sampleInfo)
if( !sdFile.exists() ){
	exit 1, "The specified sample description file does not exist = ${params.sampleInfo}"
}
log.info "Checking sample description file = $sdFile"

featuresFile = file(params.annotations)
if( !featuresFile.exists() ){
	exit 1, "The specified feature annotation file does not exist = ${params.annotations}"
}
log.info "Checking feature annotation file = $featuresFile"

calTPM = file(params.rsTPM)
if( !calTPM.exists() ){
	exit 1, "The specified R script file does not exist = ${params.rsTPM}"
}
log.info "Checking R script to calculate TPM values = $calTPM"

gnmDir = file(params.genomeDirPath)
if( !gnmDir.exists() ){
	exit 1, "The specified genome annotation directory does not exist = ${params.genomeDirPath}"
}
log.info "Checking genome annotation directory = $gnmDir"

resDir = file(params.outDir)
if( !resDir.exists() && !resDir.mkdirs() ){
	exit 1, "The specified directory to save results cannot be created = ${params.outDir}\n Check file system access permission"
}
log.info "Checking results directory = $resDir"

if( !(params.pairedEnd in ['yes', 'no']) ){
	exit 1, "Invalid value to indicate if data is paired-end or single-end = ${params.pairedEnd}"
}

/*
 Program versions
*/

process get_program_versions{
	publishDir "${resDir}/software", mode: 'move'

	output:
	file('programs_version.txt') into programs_version

	"""
	echo nextflow ${nextflow.version} > tmp_version.txt
	fastqc -v >> tmp_version.txt
	echo trim_galore \$(trim_galore -v | awk '/version/{print\$2}') >> tmp_version.txt
	echo STAR \$(STAR --version) >> tmp_version.txt
	samtools --version | grep samtools >> tmp_version.txt
	picard SortSam --version 2>&1 | awk '{sub(/-SNAPSHOT/,"");print"picard "\$1}' >> tmp_version.txt
	featureCounts -v 2>&1 | grep feature >> tmp_version.txt
	R --version | head -1 | cut -d "(" -f 1 >> tmp_version.txt
	bamCoverage --version >> tmp_version.txt
	sort tmp_version.txt > programs_version.txt
	"""

}

/*
 Create channels with the fastq files for processing and quality control
*/

if(params.pairedEnd == "yes"){

	// Reads1
	Channel
		.fromPath(sdFile)
		.splitCsv(header:true)
		.map{ row -> tuple(row.sample, file(row.reads1)) }
		.set { samples_r1 }

	// Reads2

	Channel
		.fromPath(sdFile)
		.splitCsv(header:true)
		.map{ row -> tuple(row.sample, file(row.reads2)) }
		.set { samples_r2 }

	// Both reads

	Channel
		.fromPath(sdFile)
		.splitCsv(header:true)
		.map{ row -> def key = row.sample; def reads = row[1, 2] 
		return tuple(key.toString(), reads) }
		.set { samples_quality }

} else if(params.pairedEnd == "no"){

	// Reads1
	Channel
		.fromPath(sdFile)
		.splitCsv(header:true)
		.map{ row -> tuple(row.sample, file(row.reads1)) }
		.into { samples_r1; samples_r2; samples_quality }

}

/*
 Step 0. Quality check
*/

process fastq_quality{
	publishDir "${resDir}/qc/fastqc", mode: 'copy'

	input:
	set val(name), file(reads) from samples_quality

	output:
	file("${name}_fastqc") into fastqc_results

	"""
	mkdir ${name}_fastqc
	fastqc -o ${name}_fastqc -q -t ${params.numCPUs} ${reads}
	"""

}

/*
 Step 1. Split fastq files
*/

process split_Reads1{
	label 'fastq_splitting'

	input:
	set val(name), file(reads) from samples_r1

	output:
	set val(name), file('*.gz') into samples_r1_split

	"""
	zcat ${reads} | split --numeric-suffixes=1 -a 4 -l \$((${params.chunkSize} * 4)) --filter='gzip > \$FILE.gz' - "${name}_reads1_"
	"""

}

process split_Reads2{
	label 'fastq_splitting'

	input:
	set val(name), file(reads) from samples_r2

	when:
	params.pairedEnd == "yes"

	output:
	set val(name), file('*.gz') into samples_r2_split

	"""
	zcat ${reads} | split --numeric-suffixes=1 -a 4 -l \$((${params.chunkSize} * 4)) --filter='gzip > \$FILE.gz' - "${name}_reads2_"
	"""

}

if(params.pairedEnd == "yes"){

	samples_r1_split
		.combine(samples_r2_split, by: 0)
		.transpose()
		.map{ row -> def key = row[0]; def reads = row[1, 2] 
		return tuple(key.toString(), reads) }		
		.set { samples }

} else if(params.pairedEnd == "no"){

	samples_r1_split
		.transpose()
		.set { samples }

}

/*
 Step 2. Trimming
*/

process trim_reads{
	publishDir "${resDir}/qc/trimming", mode: 'copy', pattern: "${name}/*report.txt"

	input:
	set val(name), file(reads) from samples

	output:
	set val(name), file("${name}/*.fq.gz") into trimmed_samples
	file("${name}/*report.txt") into trim_stats

	script:
	if(params.pairedEnd == "yes"){
		"""
		trim_galore --paired ${reads} -o ${name}
		"""
	} else if(params.pairedEnd == "no"){
		"""
		trim_galore ${reads} -o ${name}
		"""
	}

}

/*
 Step 3. Mapping
*/

process read_mapping{
	publishDir "${resDir}/qc/mapping", mode: 'copy', pattern: '*Log.final.out'

	input:
	set val(name), file(trimmed_reads) from trimmed_samples

	output:
	set val(name), file('*.bam') into split_mapping
	file('*Log.final.out') into mapping_stats

	script:
	if(params.pairedEnd == "yes"){
		"""
		file=\$(echo ${trimmed_reads[0]} | sed s/_1_sequence.fastq.gz//)
		STAR --readFilesIn ${trimmed_reads} --outFileNamePrefix \$file. --outTmpDir ${params.tmpOutDir}/\$file \
		--readFilesCommand zcat --runThreadN ${params.numCPUs} --genomeDir ${gnmDir} --sjdbGTFfile ${featuresFile} \
		--sjdbOverhang 150 --outFilterMultimapNmax 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.06 \
		--alignIntronMax 500000 --alignMatesGapMax 500000 --alignEndsType EndToEnd \
		--outSAMattributes NH HI NM MD --outSAMtype BAM Unsorted
		"""
	} else if(params.pairedEnd == "no"){
		"""
		file=\$(echo ${trimmed_reads[0]} | sed s/_1_sequence.fastq.gz//)
		STAR --readFilesIn ${trimmed_reads} --outFileNamePrefix \$file. --outTmpDir ${params.tmpOutDir}/\$file \
		--readFilesCommand zcat --runThreadN ${params.numCPUs} --genomeDir ${gnmDir} --sjdbGTFfile ${featuresFile} \
		--sjdbOverhang 51 --outFilterMultimapNmax 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.06 \
		--alignIntronMax 500000 --alignEndsType EndToEnd \
		--outSAMattributes NH HI NM MD --outSAMtype BAM Unsorted
		"""

	}

}

/*
 Step 4. Filter bam files to remove singletons and/or mitochondrial reads
*/

process mapped_filter{
	publishDir "${resDir}/qc/bam_filter", mode: 'copy', pattern: '*stats'

	input:
	set val(name), file(mapped_reads) from split_mapping

	output:
	set val(name), file('*rmChrM.bam') into split_mapping_filtered
	file('*stats') into filtering_stats

	script:
	if(params.pairedEnd == "yes"){
		"""
		samtools flagstat ${mapped_reads} > ${mapped_reads.baseName}_rmSE.stats
		samtools view -b -f 0x2 ${mapped_reads} > ${mapped_reads.baseName}_rmSE.bam
		samtools flagstat ${mapped_reads.baseName}_rmSE.bam > ${mapped_reads.baseName}_rmSE_post.stats
		samtools view -h ${mapped_reads.baseName}_rmSE.bam | grep -v chrM | samtools view -bS - > ${mapped_reads.baseName}_rmChrM.bam
		samtools flagstat ${mapped_reads.baseName}_rmChrM.bam > ${mapped_reads.baseName}_rmChrM.stats
		"""
	} else if(params.pairedEnd == "no"){
		"""
		samtools flagstat ${mapped_reads} > ${mapped_reads.baseName}_rmSE.stats
		samtools view -h ${mapped_reads} | grep -v chrM | samtools view -bS - > ${mapped_reads.baseName}_rmChrM.bam
		samtools flagstat ${mapped_reads.baseName}_rmChrM.bam > ${mapped_reads.baseName}_rmChrM.stats
		"""
	}

}

/*
 Step 5. Merge and sort BAM files per sample
*/

// Organise files

split_mapping_filtered.map { row -> def key = row[0]; def bamFiles = row[1]
	return tuple(key.toString(), bamFiles) }
	.groupTuple()
	.set { filtered_bams }

// Process files

process sort_bams{
	publishDir "${resDir}/mapping", mode: 'copy', pattern: '*bam'
	publishDir "${resDir}/qc/bam_filter", mode: 'copy', pattern: '*stats'

	input:
	set val(name), file(bam_files) from filtered_bams

	output:
	set val(name), file('*sorted.bam') into bam_merge, bam_track_fwd, bam_track_rv
	file('*stats') into merging_stats

	"""
	# Merge
	samtools merge -f ${name}.bam ${bam_files}
	# Sort
	picard -Xmx5g -Djava.io.tmpdir=${params.tmpOutDir} SortSam I=${name}.bam O=${name}_sorted.bam SORT_ORDER=coordinate
	# Final bam stats
	samtools flagstat ${name}_sorted.bam > ${name}_postMerging.stats
	"""

}

/*
 Step 6. Read counting
*/

process read_counts{
	publishDir "${resDir}/quant", mode: 'copy', pattern: '*.counts'
	publishDir "${resDir}/qc/read_counting", mode: 'copy', pattern: '*.summary'
	label 'counting'

	input:
	set val(name), file(bam_fC) from bam_merge

	output:
	file("*.counts") into count_files
	file("*.summary") into counting_summary

	script:
	if(params.pairedEnd == "yes"){
		"""
		featureCounts -C -p -s 2 -t exon -g gene_id -a ${featuresFile} -o ${bam_fC.baseName}.counts ${bam_fC}
		"""
	} else if(params.pairedEnd == "no"){
		"""
		featureCounts -s 2 -t exon -g gene_id -a ${featuresFile} -o ${bam_fC.baseName}.counts ${bam_fC}
		"""
	}

}

// Normalize read counts

process normalized_counts{
	publishDir "${resDir}/quant", mode: 'copy'
	label 'counting'

	input:
	file(raw_counts) from count_files.collect()

	output:
	file("*counts*.txt") into count_tables

	"""
	awk 'NR > 1 {print \$1"\t"\$6}' ${raw_counts[0]} > table_raw_counts.txt
	for i in ${raw_counts}; do cp table_raw_counts.txt tmp_table.txt; \
		awk 'NR > 1 {print\$7}' \$i | paste tmp_table.txt - > table_raw_counts.txt; \
		rm tmp_table.txt; done
	Rscript ${calTPM} table_raw_counts.txt
	"""

}

/*
 Step 7. Generate normalized signal tracks
*/

// Generate tracks

process signal_tracks_fwd{
	publishDir "${resDir}/signal", mode: 'copy'
	label 'signal_tracks'

	input:
	set val(name), file(bam_cov) from bam_track_fwd

	output:
	file("*.bw") into bigWig_files_fwd

	"""
	samtools index ${bam_cov}
	bamCoverage --bam ${bam_cov} --outFileName ${bam_cov.baseName}_fwd.bw --binSize 1 --normalizeUsing CPM \
		--filterRNAstrand forward --numberOfProcessors ${params.numCPUs_Dtools} --outFileFormat bigwig
	"""

}

process signal_tracks_rv{
	publishDir "${resDir}/signal", mode: 'copy'
	label 'signal_tracks'

	input:
	set val(name), file(bam_cov) from bam_track_rv

	output:
	file("*.bw") into bigWig_files_rv

	"""
	samtools index ${bam_cov}
	bamCoverage --bam ${bam_cov} --outFileName ${bam_cov.baseName}_rv.bw --binSize 1 --normalizeUsing CPM \
		--filterRNAstrand reverse --numberOfProcessors ${params.numCPUs_Dtools} --outFileFormat bigwig
	"""

}

workflow.onComplete {
	log.info ( workflow.success ? "\nDone! Main results are saved in ${resDir}\n" : "There was an error during the execution, check log files." )
}


