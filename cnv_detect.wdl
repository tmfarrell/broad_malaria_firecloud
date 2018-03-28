## 
## cnv_detect.wdl 
## 

## WORKFLOW
workflow CNVDetect { 
    File ref
    File ref_dict
    File ref_index
    String output_dir
    File sample_bam_paths_file

    Boolean align 

    String cnv_detect_docker = "tfarrell/cnv_detect:0.1"
    String cnv_detect_path = "/usr/local/etc/cnv_detect/main.py"

    # get sample_sam_paths 
    Array[File] sample_bams = read_lines(sample_bam_paths_file)

    # align, if applicable
    if (align) { 
        scatter (sample_bam in sample_bams) { 
            call AlignBam { 
                input: 
                ref = ref, 
                ref_dict = ref_dict,
                ref_index = ref_index,
                bam = sample_bam,
                gitc_docker = gitc_docker
            }
        }
    }

    # assumes if not aligning, then sample_bams are already aligned 
    Array[File] aligned_bams = select_first([AlignBam.aligned_bam, sample_bams])

    call ComputeExpectedGC { 
        input: 
        ref = ref, 
        cnv_detect_path = cnv_detect_path,
        cnv_detect_docker = cnv_detect_docker
    }

    scatter (aligned_bam in aligned_bams) { 
        call ComputeGCCorrection {     
            input: 
            ref = ref, 
            bam = bam, 
            cnv_detect_path = cnv_detect_path, 
            cnv_detect_docker = cnv_detect_docker, 
            expected_gc_file = ComputeExpectedGC.expected_gc_file
        }

        call ComputeCorrectedPileup { 
            input: 
            ref = ref, 
            bam = bam, 
            gc_correction_file = ComputeGCCorrection.gc_correction_file,
            cnv_detect_path = cnv_detect_path,
            cnv_detect_docker = cnv_detect_docker
        }

        call ComputeWIGFile { 
            input: 
            pileup = ComputeCorrectedPileup.corrected_pileup,
            cnv_detect_path = cnv_detect_path,
            cnv_detect_docker = cnv_detect_docker,
        }

        call CallReadDepthCNVs { 
            input: 
            ref = ref, 
            wig_file = ComputeWIGFile.wig_file,
            cnv_detect_path = cnv_detect_path,
            cnv_detect_docker = cnv_detect_docker
        }

        call ReadPairDiscordanceCNVs { 
            input: 
            bam = bam,  
            cnv_detect_path = cnv_detect_path, 
            cnv_detect_docker = cnv_detect_docker
        }

        call ComputeFinalCalls { 
            input:  
            read_depth_cnvs = CallReadDepthCNVs.cnvs,
            inverted_read_pair_cnvs = ReadPairDiscordanceCNVs.inverted_cnvs,
            spanning_read_pair_cnvs = ReadPairDiscordanceCNVs.spanning_cnvs,
            cnv_detect_path = cnv_detect_path,
            cnv_detect_docker = cnv_detect_docker
        }
    } 

    output { 
        Array[File] final_cnvs = ComputeFinalCalls.cnvs
    }
}

## TASKS
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

        # filter/ clean bam 
        gatk4/gatk-launch PrintReads \
            -I ${bam} \
            -O ${sample_name}.filtered.bam \
            --create-output-bam-index
    }

    output { 
        File aligned_bam = "${sample_name}.filtered.bam"
        File aligned_bam_index = "${sample_name}.filtered.bai"
    }

    runtime { 
        docker: gitc_docker
        memory: mem_size
        disks: "local-disk " + disk_size + " HDD"
    }
}

task ComputeExpectedGC { 
    File ref

    Int disk_size
    String mem_size
    String cnv_detect_path
    String cnv_detect_docker

    String ref_prefix = basename(ref, ".fasta")

    command { 
        python ${cnv_detect_path} expected prodgc \
            --prefix ${ref_prefix} \
            ${ref}
    }

    runtime { 
        docker: cnv_detect_docker
        memory: mem_size
        disks: "local-disk " + disk_size + " HDD"
    }

    output { 
        File expected_gc_file = "${ref_prefix}_prodGC.e"
    }
}

task ComputeGCCorrection { 
    File ref
    File bam 
    File expected_gc_file

    Int disk_size
    String mem_size
    String cnv_detect_path
    String cnv_detect_docker

    String sample_name = basename(bam, ".bam")

    command { 
        python ${cnv_detect_path} prodgc \
            -p ${sample_name} \
            -f ${ref} -e ${expected_gc_file} \
            ${bam}
    }        

    runtime { 
        docker: cnv_detect_docker
        memory: mem_size
        disks: "local-disk " + disk_size + " HDD"
    } 

    output { 
        File gc_correction_file = "${sample_name}.cfile"
    }
}

task ComputeCorrectedPileup { 
    File ref
    File bam 
    File gc_correction_file

    Int disk_size
    String mem_size
    String cnv_detect_path
    String cnv_detect_docker

    String sample_name = basename(bam, ".bam")

    command { 
        python ${cnv_detect_path} pileup \
            -p ${sample_name} ${sample_name}.cfile \
            ${bam}
    }        

    runtime { 
        docker: cnv_detect_docker
        memory: mem_size
        disks: "local-disk " + disk_size + " HDD"
    } 

    output { 
        File corrected_pileup = "${sample_name}.pileup"
    }
}

task ComputeWIGFile { 
    File pileup

    Int disk_size
    String mem_size
    String cnv_detect_path
    String cnv_detect_docker

    String sample_name = basename(pileup, ".pileup")

    command { 
        python ${cnv_detect_path} wig \
            -p ${sample_name} ${sample_name}.pileup
    }        

    runtime { 
        docker: cnv_detect_docker
        memory: mem_size
        disks: "local-disk " + disk_size + " HDD"
    } 

    output { 
        File wig_file = "${sample_name}.wig"
    }
}

task CallReadDepthCNVs { 
    File ref
    File wig_file

    Int disk_size
    String mem_size
    String cnv_detect_path
    String cnv_detect_docker

    String sample_name = basename(wig_file, ".wig")

    command {
        # gives cnvs for each chromosome 
        python ${cnv_detect_path} segmentation meanshift \
            -p ${sample_name} -f ${ref} ${sample_name}.wig

        # join 
        cat ${sample_name}.*.bed > ${sample_name}.rd_cnvs.bed
    }        

    runtime { 
        docker: cnv_detect_docker
        memory: mem_size
        disks: "local-disk " + disk_size + " HDD"
    } 

    output { 
        File cnvs = "${sample_name}.rd_cnvs.bed"
    }
}

task ReadPairDiscordanceCNVs { 
    File bam 

    Int disk_size
    String mem_size
    String cnv_detect_path
    String cnv_detect_docker

    String sample_name = basename(bam, ".bam")

    command { 
        # inverted discordant read pairs 
        python ${cnv_detect_path} pe inverted \
            -b -n ${sample_name}.inverted.bed \
            ${bam}

        # spanning discordant read pairs 
        python ${cnv_detect_path} pe spanning \
            -b -n ${sample_name}.spanning.bed \
            ${bam}
    }        

    runtime { 
        docker: cnv_detect_docker
        memory: mem_size
        disks: "local-disk " + disk_size + " HDD"
    } 

    output { 
        File spanning_cnvs = "${sample_name}.spanning.bed"
        File inverted_cnvs = "${sample_name}.inverted.bed"
    }
}

task ComputeFinalCalls { 
    File read_depth_cnvs
    File inverted_read_pair_cnvs
    File spanning_read_pair_cnvs

    Int disk_size
    String mem_size
    String cnv_detect_path
    String cnv_detect_docker

    String sample_name = basename(read_depth_cnvs, ".rd_cnvs.bed")

    command { 
        python ${cnv_detect_path} intersect \
            -s ${spanning_read_pair_cnvs} \
            -i ${inverted_read_pair_cnvs} \
            ${read_depth_cnvs} > ${sample_name}.intersect.bed
    }        

    runtime { 
        docker: docker_image
        memory: mem_size
        disks: "local-disk " + disk_size + " HDD"
    } 

    output { 
        File cnvs = "${sample_name}.intersect.bed"
    }
}