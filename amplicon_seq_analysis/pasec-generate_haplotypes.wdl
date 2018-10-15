##
##	generate_haplotypes.wdl
## 
## 	  A subworkflow called by pasec.wdl that generates 
##    haplotypes for each input sample bam, in parallel  
##
## 	tfarrell@broadinstitute.org
##  20180312 
## 

workflow generate_amplicon_haplotypes { 
	# inputs 
	File bed
	String amplicon
	String bed_index
	String bwa_aligner
	Array[File] bam_paths

	# optional 
	File? mask_file 
	
	# compute read metrics?  
	Boolean do_compute_read_metrics 

	# clustering/ filtering params 
	Int min_haplotype_cov
	Float min_haplotype_freq
	Int haplotype_cluster_dist
	Float haplotype_cluster_cov_ratio

	# docker/ tool paths  
	String pasec_docker
	String generate_haplotypes_script_path
	String compute_plot_read_metrics_path

	# memory/ disk sizes 
	Int reg_mem_size_gb
	Int large_mem_size_gb
	Int reg_disk_size
	Int large_disk_size  

	scatter (bam_path in bam_paths) { 

		if (!do_compute_read_metrics) { 
		   call generate_haplotypes { 
			   	input:
				bed = bed,
				bam = bam_path, 
				amplicon = amplicon, 
				bed_index = bed_index,
				mask_file = mask_file,
				bwa_aligner = bwa_aligner, 
				pasec_docker = pasec_docker,
				bam_index = bam_path + ".bai",
				min_haplotype_cov = min_haplotype_cov,
				min_haplotype_freq = min_haplotype_freq,
				haplotype_cluster_dist = haplotype_cluster_dist,
				haplotype_cluster_cov_ratio = haplotype_cluster_cov_ratio,
				script_path = generate_haplotypes_script_path, 
				mem_size_gb = reg_mem_size_gb,
				disk_size = reg_disk_size
		    }
		}
		if (do_compute_read_metrics) { 
		   call generate_haplotypes_reads { 
			   	input:
				bed = bed,
				bam = bam_path, 
				amplicon = amplicon, 
				bed_index = bed_index,
				mask_file = mask_file, 
				bwa_aligner = bwa_aligner,
				pasec_docker = pasec_docker, 
				bam_index = bam_path + ".bai",
				min_haplotype_cov = min_haplotype_cov,
				min_haplotype_freq = min_haplotype_freq,
				haplotype_cluster_dist = haplotype_cluster_dist,
				haplotype_cluster_cov_ratio = haplotype_cluster_cov_ratio,
				script_path = generate_haplotypes_script_path,
				mem_size_gb = reg_mem_size_gb,
				disk_size = reg_disk_size
		    }
		    call compute_read_metrics { 
				input:
				amplicon = amplicon, 
                pasec_docker = pasec_docker, 
				reads_file = generate_haplotypes_reads.reads_file, 
				haplotypes_file = generate_haplotypes_reads.haplotypes_file,
				script_path = compute_plot_read_metrics_path,
				mem_size_gb = large_mem_size_gb,
				disk_size = large_disk_size
		    }
		}  
	}

	output { 
		Array[File?] reads_files = generate_haplotypes_reads.reads_file
		Array[File?] read_metrics_files = compute_read_metrics.read_metrics_file
		Array[File] haplotypes_files = select_all(flatten([generate_haplotypes.haplotypes_file, generate_haplotypes_reads.haplotypes_file]))
	} 
}

task generate_haplotypes { 
	File bam
	File bed
	File bam_index
	String amplicon
	String bed_index 
	String script_path
	String bwa_aligner

	File? mask_file 
	
	String min_haplotype_cov
	String min_haplotype_freq
	Int haplotype_cluster_dist
	Float haplotype_cluster_cov_ratio 

	String id = basename(bam, ".bam")

	String pasec_docker
	Int mem_size_gb
	Int disk_size 

	command { 
		# if mask file provided
		if [ -f "${mask_file}" ] ; then 
		   ruby ${script_path} haplotypes \
			--bam ${bam} \
			--mask_bed_file ${mask_file} \
			--bwa_aligner_used ${bwa_aligner} \
			--bed ${bed} --bed_index ${bed_index} \
			--min_haplotype_coverage ${min_haplotype_cov} \
			--min_haplotype_freq ${min_haplotype_freq} \
			--haplotype_clustering-edit_dist ${haplotype_cluster_dist} \
			--haplotype_clustering-coverage_ratio ${haplotype_cluster_cov_ratio} \
			> ${id}.${amplicon}.haplotypes.tsv			   
		else
		   ruby ${script_path} haplotypes \
			--bam ${bam} \
			--bwa_aligner_used ${bwa_aligner} \
			--bed ${bed} --bed_index ${bed_index} \
			--min_haplotype_coverage ${min_haplotype_cov} \
			--min_haplotype_freq ${min_haplotype_freq} \
			--haplotype_clustering-edit_dist ${haplotype_cluster_dist} \
			--haplotype_clustering-coverage_ratio ${haplotype_cluster_cov_ratio} \
			> ${id}.${amplicon}.haplotypes.tsv	
		fi
		# if no haplotypes generated, touch empty file 
		if [ ! -f ${id}.${amplicon}.haplotypes.tsv ] ; then 
		   touch ${id}.${amplicon}.haplotypes.tsv
		fi 
	} 

	output { 
		File haplotypes_file = "${id}.${amplicon}.haplotypes.tsv"
	}

	runtime { 
        preemptible: 4
        docker: pasec_docker 
        memory: mem_size_gb + " GB" 
        disks: "local-disk " + disk_size + " HDD"
    }
}

task generate_haplotypes_reads { 
	File bam
	File bed
	File bam_index
	String bed_index 
	String amplicon
	String script_path

	File? mask_file 

	String bwa_aligner
	String min_haplotype_cov
	String min_haplotype_freq
	Int haplotype_cluster_dist
	Float haplotype_cluster_cov_ratio 

	String id = basename(bam, ".bam")

	String pasec_docker
	Int mem_size_gb
	Int disk_size 

	command { 
		if [ -f "${mask_file}" ] ; then 
		   ruby ${script_path} haplotypes \
			--bam ${bam} \
			--print_reads --print_to . \
			--mask_bed_file ${mask_file} \
			--bwa_aligner_used ${bwa_aligner} \
			--bed ${bed} --bed_index ${bed_index} \
			--min_haplotype_coverage ${min_haplotype_cov} \
			--min_haplotype_freq ${min_haplotype_freq} \
			--haplotype_clustering-edit_dist ${haplotype_cluster_dist} \
			--haplotype_clustering-coverage_ratio ${haplotype_cluster_cov_ratio}
		else
		   ruby ${script_path} haplotypes \
			--bam ${bam} \
			--print_reads --print_to . \
			--bwa_aligner_used ${bwa_aligner} \
			--bed ${bed} --bed_index ${bed_index} \
			--min_haplotype_coverage ${min_haplotype_cov} \
			--min_haplotype_freq ${min_haplotype_freq} \
			--haplotype_clustering-edit_dist ${haplotype_cluster_dist} \
			--haplotype_clustering-coverage_ratio ${haplotype_cluster_cov_ratio}
		fi
		if [ ! -f ${id}.${amplicon}.haplotypes.tsv ] ; then 
		   touch ${id}.${amplicon}.haplotypes.tsv
		fi	 
		if [ ! -f ${id}.${amplicon}.reads.tsv ] ; then
		   touch ${id}.${amplicon}.reads.tsv
		fi 
	} 

	output { 
		File haplotypes_file = "${id}.${amplicon}.haplotypes.tsv"
		File reads_file = "${id}.${amplicon}.reads.tsv"
	}

	runtime { 
        preemptible: 4
        docker: pasec_docker 
        memory: mem_size_gb + " GB" 
        disks: "local-disk " + disk_size + " HDD"
    }
} 
 
task compute_read_metrics { 
	String amplicon 
	String script_path     

	File reads_file
	File haplotypes_file

	String id = basename(reads_file, "." + amplicon + ".reads.tsv")
	String outfile = id + "." + amplicon + ".read_metrics.csv"

	String pasec_docker
	Int mem_size_gb
	Int disk_size 
     
	command { 
		python ${script_path} \
			--reads_file ${reads_file} \
			--haplotypes_file ${haplotypes_file} \
			--output_dir . \
			--do_not_plot
		if [ ! -f ${outfile} ] ; then 
			touch ${outfile}
		fi 	     
	} 

    output { 
     	File read_metrics_file = "${outfile}"
    }  

    runtime { 
        preemptible: 4
        docker: pasec_docker 
        memory: mem_size_gb + " GB" 
        disks: "local-disk " + disk_size + " HDD"
    }  
} 
