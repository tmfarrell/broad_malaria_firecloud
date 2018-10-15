#
# pasec.wdl
#  
#   An amplicon sequencing data processing workflow for use on FireCloud  
# 
#   For each bam_path in bam_paths_file, in parallel: 
# 		a) sorts bam and converts bam to fastq 
#    	b) merges paired-end reads using FLASH 
#    	c) realigns with bwa and generates sam 
#    	d) converts sam to bam, sorts and indexes 
#    	e) computes haplotype coverage for each amplicon in bed file 
# 
#
# Tim Farrell
# tfarrell@broadinstitute.org
# Broad Institute, IDMP, Malaria 
# 20181009
#

#import "/seq/plasmodium/tfarrell/rtss/code/pasec/src/wdl/generate_haplotypes.wdl" as generate_haplotypes
import "https://api.firecloud.org/ga4gh/v1/tools/broad-malaria-firecloud-methods:pasec-generate_haplotypes/versions/3/plain-WDL/descriptor" as generate_haplotypes

## WORKFLOW DEFINITION
workflow pasec {
	# data input
	File ref
    File ref_amb
    File ref_ann
    File ref_bwt
    File ref_pac
    File ref_sa
    File ref_index
	File bam_paths_file
	File amplicon_bed_file
	File run_seq_metadata
	File run_sample_metadata
	Array[String] amplicons
	Map[String, String] amplicon_bed_index_map

	# optional 
	File? mask_file
	File? known_haplotypes

	# output params
	Boolean do_compute_read_metrics

	# filtering/ clustering configs
	# each a map like { amplicon: value } 
	Map[String, Int] min_haplotype_cov
	Map[String, Float] min_haplotype_freq
	Map[String, Int] haplotype_cluster_dist
	Map[String, Float] haplotype_cluster_cov_ratio

	# dockers/ tool paths
	String gitc_docker         # genomes-in-the-cloud (gitc) docker
    String gitc_path_to_bwa
    String flash_docker
	String flash_path_to_flash
	String flash_path_to_samtools
	String pasec_docker 
	String generate_haplotypes_script_path 
	String compute_read_metrics_script_path
	String process_haplotypes_script_path 

	# memory/ disk sizes 
	Int reg_mem_size_gb
	Int large_mem_size_gb
	Int reg_disk_size
	Int large_disk_size

	# workflow control
	Boolean realign

	# parse bam_paths_file 
	Array[File] bam_paths = read_lines(bam_paths_file) 

	if (realign == true) { 
		# run alignment pipeline on each bam_path in parallel 
		scatter (bam_path in bam_paths) { 
			String id = basename(bam_path, '.bam')
			# merge paired-end	       
		 	call merge_paired_end_reads { 
		 	    input: 
		 	    id = id,
		 	    bam = bam_path, 
		 	    flash = flash_path_to_flash,
		 	    samtools = flash_path_to_samtools, 
		 	    flash_docker = flash_docker,
		 	    mem_size_gb = reg_mem_size_gb, 
				disk_size = reg_disk_size
		 	}  
		 	# alignment 
			call bwa_mem { 
			    input: 
			    id = id,
			    ref = ref,
                amb = ref_amb,
                ann = ref_ann, 
                bwt = ref_bwt, 
                pac = ref_pac, 
                sa = ref_sa, 
                ref_index = ref_index, 
			    gitc_docker = gitc_docker, 
			    gitc_path_to_bwa = gitc_path_to_bwa,
				fastq = merge_paired_end_reads.merged_fastq,
				mem_size_gb = large_mem_size_gb, 
				disk_size = large_disk_size
			}
		}
	}

	Array[File] aligned_bams = select_first([bwa_mem.aligned_bam, bam_paths])

	# loop over amplicons 
	scatter (amplicon in amplicons) { 
		# generate haplotypes
		call generate_haplotypes.generate_amplicon_haplotypes { 
			input: 
			amplicon = amplicon,
			bwa_aligner = "mem",  
			mask_file = mask_file,
			bed = amplicon_bed_file,
			bam_paths = aligned_bams,	
			pasec_docker = pasec_docker, 
			bed_index = amplicon_bed_index_map[amplicon],
			min_haplotype_cov = min_haplotype_cov[amplicon],
			min_haplotype_freq = min_haplotype_freq[amplicon],
			haplotype_cluster_dist = haplotype_cluster_dist[amplicon],
			haplotype_cluster_cov_ratio = haplotype_cluster_cov_ratio[amplicon],
			do_compute_read_metrics = do_compute_read_metrics,
			generate_haplotypes_script_path = generate_haplotypes_script_path, 
			compute_plot_read_metrics_path = compute_read_metrics_script_path,
			reg_mem_size_gb = reg_mem_size_gb,
			large_mem_size_gb = large_mem_size_gb,
			reg_disk_size = reg_disk_size,
			large_disk_size = large_disk_size
		}
		# analyze run 
		call analyze_run { 
			input: 
			amplicon = amplicon,
			pasec_docker = pasec_docker,
			seq_metadata = run_seq_metadata,
			sample_metadata = run_sample_metadata, 
			known_haplotypes = known_haplotypes, 
			analysis_script_path = process_haplotypes_script_path, 
			sample_haplotypes_files = generate_amplicon_haplotypes.haplotypes_files,
			mem_size_gb = large_mem_size_gb,
			disk_size = large_disk_size
		}
	}
	# pool different amplicon analyses together
	call pool_amplicon_analyses { 
		input: 
		seq_index_files = analyze_run.seq_index, 
		filter_summary_files = analyze_run.filter_summary,
		haplotype_index_files = analyze_run.haplotype_index,
		read_metrics_files = flatten(generate_amplicon_haplotypes.read_metrics_files),
        gitc_docker = gitc_docker, 
		mem_size_gb = large_mem_size_gb, 
		disk_size = large_disk_size
	}

	output { 
       File seq_index = pool_amplicon_analyses.pooled_seq_index
       File filter_summary = pool_amplicon_analyses.pooled_filter_summary
       File haplotype_index_file = pool_amplicon_analyses.pooled_haplotype_index
       File read_metrics_file = pool_amplicon_analyses.read_metrics_file
	}	 
}


## TASK DEFINITIONS
# flash task 
# merge paired end reads using FLASH (https://ccb.jhu.edu/software/FLASH/)
task merge_paired_end_reads { 
	File bam
	String id
	String flash
	String samtools 
	String flash_docker

	Int mem_size_gb
    Int disk_size
	
	command { 
		# sort bam and bam2fq
		${samtools} sort -o ${id}.sorted.bam ${bam}
		${samtools} bam2fq ${id}.sorted.bam > ${id}.raw.fastq
		# merge
		${flash} -I ${id}.raw.fastq -o ${id} -d . -M 200
	} 

	output { 
		File merged_fastq = "${id}.extendedFrags.fastq"
	} 

	runtime { 
        preemptible: 4
        docker: flash_docker
        memory: mem_size_gb + " GB" 
        disks: "local-disk " + disk_size + " HDD"
    }
} 

# bwa alignment 
task bwa_mem { 
	File ref
	File sa 
	File amb
	File ann 
	File pac 
	File bwt
    File ref_index 
	String id
	File fastq
	String gitc_docker
	String gitc_path_to_bwa
	String gitc_path_to_samtools

	Int mem_size_gb
    Int disk_size
	
	command {
		# bwa mem 
		${gitc_path_to_bwa} mem -M ${ref} ${fastq} > ${id}.sam
		# sam2bam, sort and index
		${gitc_path_to_samtools} view ${id}.sam -bS > ${id}.unsorted.bam
		${gitc_path_to_samtools} sort -o ${id}.bam ${id}.unsorted.bam
		${gitc_path_to_samtools} index ${id}.bam
	} 

	output { 
		File aligned_bam = "${id}.bam"
		File aligned_bam_index = "${id}.bam.bai"
	} 

	runtime { 
        preemptible: 4
        docker: gitc_docker
        memory: mem_size_gb + " GB" 
        disks: "local-disk " + disk_size + " HDD"
    }
} 

# join sample haplotypes
task analyze_run { 
	String run_id
	String amplicon
	File seq_metadata
	File sample_metadata
	String analysis_script_path
	Array[File] sample_haplotypes_files

	# optional 
	File? known_haplotypes
	Float? min_population_freq

	String haplotypes_file = "run.haplotypes." + amplicon + ".raw.tsv"

	String pasec_docker 
	Int mem_size_gb
    Int disk_size

	command { 
		# pool raw seq haplotypes for amplicon 
		head -n 1 ${sample_haplotypes_files[0]} > ${haplotypes_file}
        for f in ${sep=" " sample_haplotypes_files} ; do 
			tail -q -n +2 $f >> ${haplotypes_file}
       	done
		# run 
		python ${analysis_script_path} \
			--run_id ${run_id} \
			--amplicon ${amplicon} \
			--run_metadata ${seq_metadata} \
			--sample_metadata ${sample_metadata} \
			--known_haplotypes_file ${default='""' known_haplotypes} \
			--min_population_freq ${default="0.00001" min_population_freq} \
			--run_haplotype_coverage_file ${haplotypes_file} \
			--output_dir . 
	}

	output { 
		File seq_index = "run.${amplicon}.index.tsv"
		File filter_summary = "${amplicon}.filter.cluster.summary.tsv"
		File haplotype_index = "haplotypes.${amplicon}.index.tsv"
	}

	runtime { 
        preemptible: 4
        docker: pasec_docker
        memory: mem_size_gb + " GB" 
        disks: "local-disk " + disk_size + " HDD"
    }
}

task pool_amplicon_analyses { 
	Array[File] seq_index_files
	Array[File] filter_summary_files
	Array[File] haplotype_index_files

	# lists of optional files 
	Array[File?] read_metrics_files

	String gitc_docker
	Int mem_size_gb
	Int disk_size 

	command { 
		# pool solexa indexes
		head -n 1 ${seq_index_files[0]} > run.index.tsv 
		tail -q -n +2 ${sep=" " seq_index_files} >> run.index.tsv
		# pool filter summary 
		head -n 1 ${filter_summary_files[0]} > filter.cluster.summary.tsv 
		tail -q -n +2 ${sep=" " filter_summary_files} >> filter.cluster.summary.tsv
		# pool haplotype indexes
		head -n 1 ${haplotype_index_files[0]} > haplotypes.index.tsv 
		tail -q -n +2 ${sep=" " haplotype_index_files} >> haplotypes.index.tsv
		# pool read metrics 
		head -n 1 ${read_metrics_files[0]} > read_metrics.csv
		tail -q -n +2 ${sep=" " read_metrics_files} >> read_metrics.csv 
	}

	output { 
		File pooled_seq_index = "run.index.tsv"
		File read_metrics_file = "read_metrics.csv"
		File pooled_filter_summary = "filter.cluster.summary.tsv"
		File pooled_haplotype_index = "haplotypes.index.tsv"
	}

	runtime { 
        preemptible: 4
        docker: gitc_docker
        memory: mem_size_gb + " GB" 
        disks: "local-disk " + disk_size + " HDD"
    }
}
