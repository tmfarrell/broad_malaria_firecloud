## 
## gatk4-germline_cnv.wdl 
## 

## WORKFLOW
workflow GATK4_GermlineCNV { 
    File ref
    File ref_dict
    File ref_index
    String run_id
    Int interval_size 
    String output_dir
    File sample_bam_paths_file
    File contig_ploidy_priors_file

    Boolean align 

    String gitc_docker
    String gatk4_docker

    # get sample_sam_paths 
    Array[File] sample_bams = read_lines(sample_bam_paths_file)

    # get intervals file
    call CreateIntervalsFile { 
        input: 
        ref = ref, 
        interval_size = interval_size
    }

    # align, if applicable
    if (align) { 
        scatter (sample_bam in sample_bams) { 
            call AlignBam { 
                input: 
                ref = ref, 
                run_id = run_id, 
                ref_dict = ref_dict,
                ref_index = ref_index,
                bam = sample_bam,
                gitc_docker = gitc_docker
            }
        }
    }

    # assumes if not aligning, then sample_bams are already aligned 
    Array[File] aligned_bams = select_first([AlignBam.aligned_bam, sample_bams])

    # for each sample_sam
    scatter(sample_bam in aligned_bams) {  
    	call CleanBam { 
        	input: 
        	bam = sample_bam, 
        	docker_image = gatk4_docker
    	} 

    	call CollectFragmentCounts { 
    		input:
        	bam = CleanBam.cleaned_bam, 
        	bam_index = CleanBam.cleaned_bam_index,
        	intervals_file = CreateIntervalsFile.intervals_file,
        	docker_image = gatk4_docker,
        	output_dir = output_dir
    	}  
    } 

    call DeterminePloidyAndCallCNVs { 
    	input: 
        run_id = run_id, 
        count_files = CollectFragmentCounts.count_file,
        contig_ploidy_priors_table = contig_ploidy_priors_file, 
        intervals_file = CreateIntervalsFile.intervals_file,  
        docker_image = gatk4_docker,  
        output_dir = output_dir
    } 
}

## TASKS
task CreateIntervalsFile { 
    File ref
    Int interval_size
    File create_intervals_script

    command { 
    	pip install biopython
        python ${create_intervals_script} -r ${ref} \
            --interval_size ${interval_size} > intervals_file.intervals
    }

    runtime { 
        docker: "continuumio/anaconda"
    }

    output { 
        File intervals_file = "intervals_file.intervals"
    }
}

task AlignBam { 
    File bam 
    File ref
    File ref_dict
    File ref_index
    String run_id
    File sample_name = basename(bam, ".bam")

    Int disk_size
    String mem_size
    String gitc_docker

    String read_group_str = "@RG\tID:${run_id}\tSM:${sample_name}\tPL:ILLUMINA"

    command { 
        # bam to fastq
        samtools bam2fq ${bam} > ${sample_name}.fastq

        # sort fastq
        cat ${sample_name}.fastq | paste - - - - \
            | sort -k1,1 -t " " | tr "\t" "\n" > ${sample_name}.sorted.fastq
        rm ${sample_name}.fastq

        # align bam 
        bwa mem -p -R ${read_group_str} ${ref} ${sample_name}.sorted.fastq \
            > ${sample_name}.aln.sam
        rm ${sample_name}.sorted.fastq

        # sort and index
        samtools view -b -o ${sample_name}.aln.bam ${sample_name}.aln.sam
        samtools sort -o ${sample_name}.aln.sorted.bam ${sample_name}.aln.bam
        samtools index ${sample_name}.aln.sorted.bam
    }

    output { 
        File aligned_bam = "${sample_name}.aln.sorted.bam"
        File aligned_bam_index = "${sample_name}.aln.sorted.bai"
    }

    runtime { 
        docker: gitc_docker
        memory: mem_size
        disks: "local-disk " + disk_size + " HDD"
    }
}

task CleanBam { 
    File bam 
    String mem_size
    Int disk_size
    String docker_image
    String sample_name = basename(bam, ".bam") 

    command { 
        gatk PrintReads -I ${bam} \
            -O ${sample_name}.filtered.bam \
            --create-output-bam-index
    } 	     

    runtime { 
        docker: docker_image
        memory: mem_size
        disks: "local-disk " + disk_size + " HDD"
    } 

    output { 
        File cleaned_bam = "${sample_name}.filtered.bam"
        File cleaned_bam_index = "${sample_name}.filtered.bai"
    } 
} 

task CollectFragmentCounts { 
    File bam
    File bam_index
    File intervals_file
    String output_dir
    String mem_size
    Int disk_size
    String docker_image
    String sample_name = basename(bam, ".bam")

    command { 
        gatk CollectFragmentCounts \
        	-I ${bam} \
            -L ${intervals_file} \
            --interval-merging-rule OVERLAPPING_ONLY \
            -O ${sample_name}.counts.hdf5
    } 

    runtime { 
        docker: docker_image
        memory: mem_size
        disks: "local-disk " + disk_size + " HDD"
    }

    output { 
    	File count_file = "${sample_name}.counts.hdf5"
    }  
}

task DeterminePloidyAndCallCNVs { 
    String run_id
    String output_dir 
    File intervals_file
    Array[File] count_files
    File contig_ploidy_priors_table

    Int disk_size
    String mem_size
    String docker_image

    command { 
        # determine contig ploidy 
        gatk DetermineGermlineContigPloidy \
        	--input ${sep=" --input " count_files} \
         	--contig-ploidy-priors ${contig_ploidy_priors_table} \
        	--output ${output_dir} \
        	--output-prefix "contig_ploidy"

        # call CNVs 
        gatk GermlineCNVCaller \
        	--run-mode COHORT \
         	-L ${intervals_file} \
        	--interval-merging-rule OVERLAPPING_ONLY \
          	--contig-ploidy-calls ${output_dir}/contig_ploidy-calls/ \
          	--input ${sep=" --input " count_files} \
          	--output ${output_dir} \
          	--output-prefix "germline_cnvs"

        # zip results
        zip -r gatk4_germline_cnvs.${run_id}.zip ${output_dir}
    } 

    runtime { 
        docker: docker_image
        memory: mem_size
        disks: "local-disk " + disk_size + " HDD"
    } 

    output { 
        File out = "gatk4_germline_cnvs.${run_id}.zip" 
    } 
} 
