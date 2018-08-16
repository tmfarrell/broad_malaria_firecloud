# 
#  gatk3_germline_snps_indels-malaria.wdl: 
#     - standardized pipeline for calling germline snps and indels using HaplotypeCaller. 
#     - developed by Malaria Group, IDMP, Broad Institute. 
#  

## WORKFLOW DEFINITION
workflow GATK3_Germline_Variants {
    ## config params
    # input data
    File ref			       # path to reference file
    File ref_dict
    File ref_index
    String run_name
    File sample_paths_file     # .tsv of (sample_name, sample_sam_path)
    File interval_files_list		   
    
    # gatk/ picard/ genomes-in-the-cloud  
    String gatk_docker         # gatk3 docker (e.g. broadinstitute/gatk3:3.8-0)
    String gatk_path_to_gatk
    String gitc_docker         # genomes-in-the-cloud (gitc) docker
    String gitc_path_to_picard
    String gitc_path_to_gatk
    
    # for base quality score recalibration
    Array[File] known_sites
    Array[File] known_sites_indices
    Array[String] bqsr_intervals_to_exclude

    # variant quality control param
    # either "vqsr" or "hard_filtering"
    String variant_qc          
    
    # vqsr params		
    # if variant_qc == "vqsr"
    # all of these params are required
    Float ts_filter_snp
    Float ts_filter_indel
    Int snp_max_gaussians
    Int indel_max_gaussians
    Int vqsr_mapping_qual_cap
    Array[String] snp_resources
    Array[String] indel_resources
    Array[String] snp_annotations
    Array[String] indel_annotations

    # hard filtering params
    # if variant_qc == "hard_filtering"
    # both of these params are required
    String snp_filter_expr
    String indel_filter_expr
    
    # snpeff 
    File Pf3D7_gff
    Boolean use_snpeff


    ## task calls 
    # run pipeline on each sample, in parallel
    scatter(sample in read_tsv(sample_paths_file)) {
        String sample_name = sample[0]
        String sample_bam_path = sample[1]	

        call MarkDuplicates {
            input:
            sample_name = sample_name,
            sorted_bam = sample_bam_path,
            picard_docker = gitc_docker,
            picard_path = gitc_path_to_picard
        }

        call ReorderBam {
            input:
            ref = ref,
            dict = ref_dict, 
            bam = MarkDuplicates.bam,
            picard_docker = gitc_docker,
            picard_path = gitc_path_to_picard
        }
    	
        # base quality score recalibration 
        call BQSR {
            input:
            ref = ref, 
            ref_dict = ref_dict, 
            ref_index = ref_index,
            bam = ReorderBam.out,
            bam_index = ReorderBam.out_index,
            sample_name = sample_name, 
            known_sites = known_sites,
            known_sites_indices = known_sites_indices,
            intervals_to_exclude = bqsr_intervals_to_exclude,
            output_table_name = sample_name + ".bqsr.table",
            output_bam_name = sample_name + ".bsqr.bam",
            gitc_docker = gitc_docker,
            gatk_path = gitc_path_to_gatk,
            picard_path = gitc_path_to_picard
        }

        call HaplotypeCaller {
            input:
            ref = ref, 
            bam = BQSR.out,
            ref_dict = ref_dict,
            ref_index = ref_index, 
            bam_index = BQSR.out_index, 
            bqsr_table = BQSR.table,
            sample_name = sample_name, 
            intervals = interval_files_list,
            gatk_docker = gatk_docker,
            gatk_path = gatk_path_to_gatk
        }
    } # end scatter

    call GenotypeGVCFs {
        input:
        ref = ref,
        ref_dict = ref_dict,
        ref_index = ref_index, 
        intervals = interval_files_list,
        vcf_files = HaplotypeCaller.vcf,
        gatk_docker = gatk_docker,
        gatk_path = gatk_path_to_gatk
    }
    
    # variant quality control
    if (variant_qc == "vqsr") { 
        # variant quality score recalibration
        # snp vqsr
        call VQSR as SnpVQSR {
            input:
            ref = ref,
            ref_dict = ref_dict,
            ref_index = ref_index, 
            intervals = interval_files_list,
            gvcf = GenotypeGVCFs.out,
            output_filename = "${run_name}.snp_vqsr.g.vcf",
            
            mode = "SNP",
            resources = snp_resources,
      		resource_files = known_sites, 
            resource_file_indices = known_sites_indices,
            annotations = snp_annotations,
            ts_filter = ts_filter_snp,
            max_gaussians = snp_max_gaussians,
            mapping_qual_cap = vqsr_mapping_qual_cap,

            gitc_docker = gitc_docker,
            gatk_path = gitc_path_to_gatk
        }
        # indel vqsr
        call VQSR as IndelVQSR {
            input:
            ref = ref,
            ref_dict = ref_dict,
            ref_index = ref_index, 
            intervals = interval_files_list,
            gvcf = SnpVQSR.out,
            output_filename = "${run_name}.indel_vqsr.g.vcf",

            mode = "INDEL",
            resources = indel_resources,
            resource_files = known_sites, 
            resource_file_indices = known_sites_indices,
            annotations = indel_annotations,
            ts_filter = ts_filter_indel,
            max_gaussians = indel_max_gaussians,
            mapping_qual_cap = vqsr_mapping_qual_cap,

            gitc_docker = gitc_docker,
            gatk_path = gitc_path_to_gatk
        }
    } 
    if (variant_qc == "hard_filtering") { 
        call HardFiltration { 
            input: 
            ref = ref, 
            ref_dict = ref_dict, 
            ref_index = ref_index, 
            vcf = GenotypeGVCFs.out,
            snp_filter_expr = snp_filter_expr,
            indel_filter_expr = indel_filter_expr,
            output_filename = "${run_name}.hard_filtered.g.vcf",
            gatk_docker = gatk_docker, 
            gatk_path = gatk_path_to_gatk
        }
    }

    # add variant annotations using SnpEff
    File vcf = select_first([IndelVQSR.out, HardFiltration.out])
    if (use_snpeff == true) {
        call SnpEff {
            input:
            vcf = vcf,
            ref = ref, 
            Pf3D7_gff = Pf3D7_gff,
            output_filename = "${run_name}.snpeff.g.vcf"
        }
    }
    
    output { 
        File gvcf = select_first([SnpEff.out, IndelVQSR.out, HardFiltration.out])
    }
}


## TASK DEFINITIONS 
# mark duplicate reads in bam 
task MarkDuplicates {
    File sorted_bam
    String sample_name

    Int disk_size
    String mem_size
    String picard_docker
    String picard_path

    command {
        java -Xmx8G -jar ${picard_path} MarkDuplicates \
            I=${sorted_bam} \
            O=${sample_name}.marked_duplicates.bam \
            M=${sample_name}.marked_duplicates.metrics
    }

    output {
        File bam = "${sample_name}.marked_duplicates.bam"
    } 

    runtime {
        docker: picard_docker
        memory: mem_size 
        disks: "local-disk " + disk_size + " HDD"
    }
}

# reorder and index a bam  
task ReorderBam {
    File ref
    File dict
    File bam
    String bam_prefix = basename(bam, '.bam')

    Int disk_size
    String mem_size
    String picard_docker
    String picard_path
    
    command {
        # reorder bam 
        java -Xmx8G -jar ${picard_path} ReorderSam \
            I=${bam} \
            O=${bam_prefix}.reordered.bam \
            R=${ref}

        # then index 
        java -Xmx8G -jar ${picard_path} BuildBamIndex \
            I=${bam_prefix}.reordered.bam
    }
    
    output {
        File out = "${bam_prefix}.reordered.bam"
        File out_index = "${bam_prefix}.reordered.bai"
    } 

    runtime {
        docker: picard_docker
        memory: mem_size 
        disks: "local-disk " + disk_size + " HDD"
    }
}

# base quality score recalibration
task BQSR {
    File ref
    File ref_dict
    File ref_index
    File bam 
    File bam_index
    String sample_name
    String output_table_name
    String output_bam_name
    Array[File] known_sites
    Array[File] known_sites_indices
    Array[String] intervals_to_exclude
    
    Int disk_size
    String mem_size
    String gitc_docker    
    String gatk_path
    String picard_path

    command {  
        # build BQSR table 
        java -Xmx4G -jar ${gatk_path} -T BaseRecalibrator -nt 1 \
            -R ${ref} -I ${bam} \
            -knownSites ${sep=" -knownSites " known_sites} \
            --excludeIntervals ${sep=" --excludeIntervals " intervals_to_exclude} \
            -o ${output_table_name}
            

        # install GATK AnalyzeCovariates R dependencies
        R --vanilla << CODE
        install.packages("gplots", repos="http://cran.us.r-project.org")
        install.packages("gsalib", repos="http://cran.us.r-project.org")
        install.packages("reshape", repos="http://cran.us.r-project.org")
        CODE

        # AnalyzeCovariates
        java -jar ${gatk_path} \
            -T AnalyzeCovariates \
            -R ${ref} \
            --BQSR ${output_table_name} \
            -plots ${sample_name}.bqsr.pdf

        # clean reads, using bqsr if applicable
        java -Xmx4G -jar ${gatk_path} \
            -T PrintReads \
            -nt 1 \
            -R ${ref} \
            -I ${bam} \
            --BQSR ${output_table_name} \
            -o ${output_bam_name}

        # build index
        java -Xmx4G -jar ${picard_path} \
            BuildBamIndex \
            I=${output_bam_name} \
            O=${output_bam_name}.bai
    }
    
    output {
        File out = "${output_bam_name}"
        File out_index = "${output_bam_name}.bai"
        File table = "${output_table_name}"
    } 

    runtime {
        docker: gitc_docker
        memory: mem_size 
        disks: "local-disk " + disk_size + " HDD"
    }
}

# call snp and indel variants 
task HaplotypeCaller {
    File bam
    File ref
    File ref_dict
    File ref_index
    File bam_index
    File bqsr_table
    File intervals
    String sample_name
    String output_name = "${sample_name}.g.vcf"

    Int disk_size
    String mem_size
    String gatk_docker
    String gatk_path 

    command {
        java -Xmx8G -jar ${gatk_path}\
            -T HaplotypeCaller \
            -nt 1 \
            -R ${ref} \
            --input_file ${bam} \
            --intervals ${intervals} \
            --BQSR ${bqsr_table} \
            -ERC GVCF \
            --interval_padding 100 \
            -o ${output_name} \
            -variant_index_type LINEAR \
            -variant_index_parameter 128000
    }

    output {
		File vcf = "${output_name}"
    } 

    runtime {
        docker: gatk_docker
        memory: mem_size 
        disks: "local-disk " + disk_size + " HDD"
    }
  }

# merge and genotype vcfs 
task GenotypeGVCFs {
    File ref
    File ref_dict
    File ref_index
	File intervals
    Boolean all_sites
	Array[File] vcf_files
    String gcvf_out = "genotypedGVCFs.vcf"

    Int disk_size
    String mem_size
    String gatk_docker
    String gatk_path 

    command {
        java -Xmx8G -jar ${gatk_path} \
            -T GenotypeGVCFs \
            -R ${ref} \
            --intervals ${sep=" --intervals " intervals} \
            -o ${gcvf_out} \
            -V ${sep=" -V " vcf_files} \
            ${true="-allSites" false="" all_sites}
    }
    output {
	   File out = gcvf_out
    } 

    runtime {
        docker: gatk_docker
        memory: mem_size 
        disks: "local-disk " + disk_size + " HDD"
    }
}

# variant quality score recalibration
# https://software.broadinstitute.org/gatk/documentation/article.php?id=2805
task VQSR {
	File ref
    File ref_dict
    File ref_index
    File gvcf 
    File intervals
    String output_filename

    String mode
    Float ts_filter
    Array[String] resources
    Array[String] annotations
	Array[File] resource_files
    Array[File] resource_file_indices
    
    Int max_gaussians
    Int mapping_qual_cap

    Int disk_size
    String mem_size
    String gitc_docker
    String gatk_path

    String vqsr_file = "${mode}.recal"
    String rscript_file = "${mode}.plots.R"
    String tranches_file = "${mode}.tranches"

	command {
        # build vqsr file
        java -Xmx8G -jar ${gatk_path} \
            -T VariantRecalibrator \
            -R ${ref} \
            --input ${gvcf} \
            --mode ${mode} \
            --recal_file ${vqsr_file} \
            --tranches_file ${tranches_file} \
            --rscript_file ${rscript_file} \
            --intervals ${sep="--intervals " intervals} \
            --resource:${sep=" --resource:" resources} \
            --use_annotation ${sep=" --use_annotation " annotations} \
            --maxGaussians ${max_gaussians} \
            --MQCapForLogitJitterTransform ${mapping_qual_cap}

        # apply vqsr
        java -Xmx8G -jar ${gatk_path} \
            -T ApplyRecalibration \
            -R ${ref} \
            --input ${gvcf} \
            --ts_filter_level ${ts_filter} \
            --tranches_file ${tranches_file} \
            --recal_file ${vqsr_file} \
            --mode ${mode} \
            -o ${mode}_vqsr.filtered.vcf
    }

    output {
        File vqsr = vqsr_file
        File rscript = rscript_file
		File tranches = tranches_file
        File out = output_filename
    } 

    runtime {
        docker: gitc_docker
        memory: mem_size 
        disks: "local-disk " + disk_size + " HDD"
    }
}

# hard-filter a vcf, if vqsr not available 
# http://gatkforums.broadinstitute.org/gatk/discussion/2806/howto-apply-hard-filters-to-a-call-set
task HardFiltration {
    File vcf
    File ref
    File ref_dict
    File ref_index
    String output_filename

    String snp_filter_expr
    String indel_filter_expr 

    Int disk_size
    String mem_size
    String gatk_docker
    String gatk_path

    command {
        # select snps
        java -Xmx8G -jar ${gatk_path} \
            -T SelectVariants \
            -R ${ref} \
            -V ${vcf} \
            -selectType SNP \
            -o raw_snps.g.vcf

        # filter snps
        java -Xmx8G -jar ${gatk_path} \
            -T VariantFiltration \
            -R ${ref} \
            -V ${vcf} \
            --filterExpression "${snp_filter_expr}" \
            --filterName "snp_filter" \
            -o filtered_snps.g.vcf

        # select indels 
        java -Xmx8G -jar ${gatk_path}\
            -T SelectVariants \
            -R ${ref} \
            -V ${vcf} \
            -selectType INDEL \
            -o raw_indels.g.vcf

        # filter indels
        java -Xmx8G -jar ${gatk_path} \
            -T VariantFiltration \
            -R ${ref} \
            -V ${vcf} \
            --filterExpression "${indel_filter_expr}" \
            --filterName "indel_filter" \
            -o filtered_indels.g.vcf

        # combine variants
        java -Xmx8G -jar ${gatk_path}\
            -T CombineVariants \
            -R ${ref} \
            --variant filtered_snps.g.vcf \
            --variant filtered_indels.g.vcf \
            -o ${output_filename} \
            --genotypemergeoption UNSORTED
    }

    output {
        File out = "${output_filename}"
    } 

    runtime {
        docker: gatk_docker
        memory: mem_size 
        disks: "local-disk " + disk_size + " HDD"
    }
}

# annotate variants
# Based on http://gatkforums.broadinstitute.org/gatk/discussion/50/adding-genomic-annotations-using-snpeff-and-variantannotator
task SnpEff {
    File vcf    
    File ref 
    File Pf3D7_gff
    String output_filename

    Int disk_size
    String mem_size
    String snpeff_path = "/opt/snpEff/snpEff.jar"
    String snpeff_docker = "maxulysse/snpeff:1.3"
    String Pf3D7_db_dir = "/opt/snpEff/data/Pf3D7/"
    String snpeff_config_path = "/opt/snpEff/snpEff.config"

    command {
        # init database
        echo "Pf3D7.genome : Plasmodium_falciparum_3D7" >> ${snpeff_config_path}
        mkdir -p ${Pf3D7_db_dir}
        mv ${ref} ${Pf3D7_db_dir}/sequences.fa
        mv ${Pf3D7_gff} ${Pf3D7_db_dir}/genes.gff

        # build db 
        java -jar ${snpeff_path} build -gff3 -v Pf3D7

        # run snpeff 
        java -Xmx4G -jar ${snpeff_path} \
        	-config ${snpeff_config_path} \
            -formatEff -no-downstream -no-intergenic \
            -no-upstream -no-utr -noStats \
            -treatAllAsProteinCoding false \
            Pf3D7 ${vcf} > ${output_filename}
    }

    output {
        File out = "${output_filename}"
    }

    runtime { 
        docker: snpeff_docker
        memory: mem_size 
        disks: "local-disk " + disk_size + " HDD"
    }
}
