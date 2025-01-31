# Author: Ali Pirani and Joseph Hale
configfile: "config/config.yaml"

import pandas as pd
import os
import numpy as np

samples_df = pd.read_csv(config["samples"])
# BARCODE = list(samples_df['barcode_id'])
# SAMPLE = list(samples_df['sample_id'])
PREFIX = config["prefix"]
# SHORTREADS = list(samples_df['sample_id'])

# don't use unicycler in most up-to-date version

# samples_df['combination'] = samples_df[['barcode_id', 'sample_id', 'sample_id']].agg('/'.join, axis=1)
# COMBINATION = list(samples_df['combination'])

SAMPLE = list(samples_df['sample_id'])

if not os.path.exists("results/"):
    os.system("mkdir %s" % "results/")

if not os.path.exists("results/" + PREFIX):
    os.system("mkdir %s" % "results/" + PREFIX)

rule all:
    input:
        #pycoqc_report = expand("results/{prefix}/pycoqc/{prefix}.html", prefix=PREFIX),
        nanopore_filtlong = expand("results/{prefix}/filtlong/{sample}/{sample}.trimmed.fastq.gz", sample=SAMPLE, prefix=PREFIX),
        #trimmed = expand("results/{prefix}/filtlong/{combination}.trimmed.fastq.gz", barcode=BARCODE, sample=SAMPLE, prefix=PREFIX, combination=COMBINATION),
        nanopore_nanoplot = expand("results/{prefix}/nanoplot/{sample}/{sample}_preqcNanoPlot-report.html", sample=SAMPLE, prefix=PREFIX),
	    #nanoplot = expand("results/{prefix}/nanoplot/{combination}_preqcNanoPlot-report.html", barcode=BARCODE, sample=SAMPLE, prefix=PREFIX, combination=COMBINATION),
        nanopore_flye_assembly = expand("results/{prefix}/flye/{sample}/{sample}_flye.fasta", sample=SAMPLE, prefix=PREFIX),
        #flye_assembly = expand("results/{prefix}/flye/{combination}_flye.fasta", barcode=BARCODE, sample=SAMPLE, prefix=PREFIX, combination=COMBINATION),
        #flye_circ_assembly = expand("results/{prefix}/flye/{combination}_flye_circ.fasta", barcode=BARCODE, sample=SAMPLE, prefix=PREFIX, combination=COMBINATION),
        nanopore_medaka = expand("results/{prefix}/medaka/{sample}/{sample}_medaka.fasta", sample=SAMPLE, prefix=PREFIX),
        #medaka_out = expand("results/{prefix}/medaka/{combination}_medaka.fasta", barcode=BARCODE, sample=SAMPLE, prefix=PREFIX, combination=COMBINATION),
        illumina_trimmomatic_pe = expand("results/{prefix}/trimmomatic/{sample}/{sample}_R1_trim_paired.fastq.gz", sample=SAMPLE, prefix=PREFIX),
        polypolish_bwa_align = expand("results/{prefix}/polypolish/{sample}/{sample}_r1.sam", sample=SAMPLE, prefix=PREFIX),
        polypolish_polish = expand("results/{prefix}/polypolish/{sample}/{sample}_polypolish.fasta", sample=SAMPLE, prefix=PREFIX),
        #polypolish = expand("results/{prefix}/polypolish/{combination}_flye_medaka_polypolish.fasta", barcode=BARCODE, sample=SAMPLE, prefix=PREFIX, combination=COMBINATION),
        #polypolish_unicycler = expand("results/{prefix}/polypolish_unicycler/{combination}_unicycler_polypolish.fasta", barcode=BARCODE, sample=SAMPLE, prefix=PREFIX, combination=COMBINATION),
        #prokka = expand("results/{prefix}/prokka/{combination}_unicycler.gff", barcode=BARCODE, sample=SAMPLE, prefix=PREFIX,combination=COMBINATION),
        quast_out = expand("results/{prefix}/quast/{sample}/{sample}_flye_medaka_polypolish/report.txt", sample=SAMPLE, prefix=PREFIX),
        #busco_out = expand("results/{prefix}/busco/{combination}.unicycler/busco_unicycler.txt", barcode=BARCODE, sample=SAMPLE, prefix=PREFIX, combination=COMBINATION),
        funannotate_sort = expand("results/{prefix}/funannotate/{sample}/{sample}_sorted.fa", sample=SAMPLE, prefix=PREFIX),
        repeatmasker = expand("results/{prefix}/repeatmasker/{sample}/{sample}_masked.fa", sample = SAMPLE, prefix = PREFIX),
        funannotate_train = expand("results/{prefix}/funannotate/{sample}/training/funannotate_train.coordSorted.bam", sample=SAMPLE, prefix=PREFIX),
        funannotate_predict = expand("results/{prefix}/funannotate/{sample}/predict_results/{sample}.proteins.fa", sample=SAMPLE, prefix=PREFIX),
        funannotate_update = expand("results/{prefix}/funannotate/{sample}/update_results/{sample}.proteins.fa", sample=SAMPLE, prefix=PREFIX),        
        interproscan = expand("results/{prefix}/funannotate/{sample}/interproscan/{sample}.proteins.fa.xml", sample=SAMPLE, prefix=PREFIX),
        eggnog = expand("results/{prefix}/funannotate/{sample}/eggnog/{sample}.emapper.annotations", sample=SAMPLE, prefix=PREFIX),
        funannotate_annotate = expand("results/{prefix}/funannotate/{sample}/annotate_results/{sample}.proteins.fa", sample=SAMPLE, prefix=PREFIX),
        busco_final = expand("results/{prefix}/busco/busco_output_prot/batch_summary.txt", prefix = PREFIX),
        multiqc_out = expand("results/{prefix}/multiqc/{prefix}_QC_report.html", prefix=PREFIX),


# rule pycoqc:
#     input:
#         seq_summary = config["long_reads"] + "/sequencing_summary.txt",
#     output:
#        "results/{prefix}/pycoqc/{prefix}.html",
#     log:
#         "logs/{prefix}/pycoqc/{prefix}.log"        
#     conda:
#         "envs/pycoqc.yaml"
#     log:
#         "logs/{prefix}/pycoqc/{barcode}/{sample}/{sample}.log"  
#     shell:
#        'pycoQC -f {input.seq_summary} -o {output} &>{log}'

# Deprecated: Porechop is no longer being supported. 
# rule trim_nano_adaptors:
#     input:
#         longreads = config["long_reads"] + "/{barcode}/",
#         #samplename = lambda wildcards: expand(str(config["short_reads"] + "/" + f"{wildcards.sample}_R1.fastq.gz")),
#     output:
#         trimmed = "results/{prefix}/porechop/{barcode}/{sample}/{sample}.trimmed.fastq"
#     log:
#         "logs/{prefix}/porechop/{barcode}/{sample}/{sample}.log"        
#     conda:
#         "envs/porechop.yaml"
#     log:
#         "logs/{prefix}/porechop/{barcode}/{sample}/{sample}.log"
#     shell:
#         'porechop -i {input.longreads} -o {output.trimmed} -t 8 --discard_middle &>{log}'

rule nanopore_filtlong:
    input:
        longreads = config["long_reads"] + "/{sample}.fastq.gz"
        #longreads = config["long_reads"] + "/{barcode}/",
        #samplename = lambda wildcards: expand(str(config["short_reads"] + "/" + f"{wildcards.sample}_R1.fastq.gz")),
    output:
        trimmed = "results/{prefix}/filtlong/{sample}/{sample}.trimmed.fastq.gz"
        #trimmed = "results/{prefix}/filtlong/{barcode}/{sample}/{sample}.trimmed.fastq.gz"
    # log:
    #     "logs/{prefix}/porechop/{barcode}/{sample}/{sample}.log"        
    # log:
    #     "logs/{prefix}/porechop/{barcode}/{sample}/{sample}.log"
    # conda:
    #     "envs/filtlong.yaml"
    # params:
    #     prefix="{sample}",
    singularity:
        "docker://staphb/filtlong:0.2.1"
    resources:
        mem_mb = 5000,
        runtime = 20
    shell:
        "filtlong --min_length 1000 --keep_percent 95 {input.longreads} --mean_q_weight 10 --target_bases 500000000 | gzip > {output.trimmed}"
        #'filtlong --min_length 1000 --keep_percent 95 {input.longreads}/*.fastq.gz | gzip > {output.trimmed}'
        #'cat {input.longreads}/*.fastq.gz > /tmp/{params.prefix}.gz | filtlong --min_length 1000 --keep_percent 95 /tmp/{params.prefix}.gz | gzip > {output.trimmed} && rm /tmp/{params.prefix}.gz'

rule nanopore_nanoplot:
    input:
        longreads = config["long_reads"] + "/{sample}.fastq.gz",
        trimmed = "results/{prefix}/filtlong/{sample}/{sample}.trimmed.fastq.gz",
    output:
        nanoplot_preqc = "results/{prefix}/nanoplot/{sample}/{sample}_preqcNanoPlot-report.html"
    #log:
    #    "logs/{prefix}/nanoplot/{sample}/{sample}.log"
    params:
        outdir="results/{prefix}/nanoplot/{sample}",
        sample="{sample}",
    singularity:
        "docker://staphb/nanoplot:1.42.0"
    #envmodules:
    #    "Bioinformatics",
    #    "nanoplot"
    resources:
        mem_mb = 5000,
        runtime = 30
    shell:
        """ 
        NanoPlot -o {params.outdir} -p {params.sample}_preqc --tsv_stats --info_in_report --N50 --title {params.sample}_preqc --fastq {input.longreads} && 
        NanoPlot -o {params.outdir} -p {params.sample}_postqc --tsv_stats --info_in_report --N50 --title {params.sample}_postqc --fastq {input.trimmed} 
        """
        # """
        # cat {input.longreads} > /tmp/{params.sample}.gz && 
        # NanoPlot -o {params.outdir} -p {params.sample}_preqc --tsv_stats --info_in_report --N50 --title {params.sample}_preqc --fastq /tmp/{params.sample}.gz && 
        # NanoPlot -o {params.outdir} -p {params.sample}_postqc --tsv_stats --info_in_report --N50 --title {params.sample}_postqc --fastq {input.trimmed} && rm /tmp/{params.sample}.gz 
        # """

# old version
# rule nanoplot:
#     input:
#         longreads = config["long_reads"] + "/{barcode}/",
#         trimmed = lambda wildcards: expand(str("results/" + f"{wildcards.prefix}" + "/filtlong/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}" + ".trimmed.fastq.gz")),
#     output:
#         nanoplot_preqc = "results/{prefix}/nanoplot/{barcode}/{sample}/{sample}_preqcNanoPlot-report.html"
#     log:
#         "logs/{prefix}/nanoplot/{barcode}/{sample}/{sample}.log"
#     params:
#         outdir="results/{prefix}/nanoplot/{barcode}/{sample}",
#         prefix="{sample}",
#     # conda:
#     #     "envs/nanoplot.yaml"
#     singularity:
#         "docker://staphb/nanoplot:1.42.0"
#     shell:
#         "cat {input.longreads}/*.fastq.gz > /tmp/{params.prefix}.gz && NanoPlot -o {params.outdir} -p {params.prefix}_preqc --tsv_stats --info_in_report --N50 --title {params.prefix}_preqc --fastq /tmp/{params.prefix}.gz && NanoPlot -o {params.outdir} -p {params.prefix}_postqc --tsv_stats --info_in_report --N50 --title {params.prefix}_postqc --fastq {input.trimmed} && rm /tmp/{params.prefix}.gz"


rule illumina_trimmomatic_pe:
    input:
        r1 = config["short_reads"] + "/" + "{sample}_R1.fastq.gz",
        r2 = config["short_reads"] + "/" + "{sample}_R2.fastq.gz",
        #r1 = lambda wildcards: expand(str(config["short_reads"] + "/" + f"{wildcards.sample}_R1.fastq.gz")),
        #r2 = lambda wildcards: expand(str(config["short_reads"] + "/" + f"{wildcards.sample}_R2.fastq.gz"))
        #r1 = lambda wildcards: expand(str(config["short_reads"] + "/" + f"{wildcards.sample}_R1_001.fastq.gz")),
        #r2 = lambda wildcards: expand(str(config["short_reads"] + "/" + f"{wildcards.sample}_R2_001.fastq.gz")),  
    output:
        r1 = "results/{prefix}/trimmomatic/{sample}/{sample}_R1_trim_paired.fastq.gz",
        r2 = "results/{prefix}/trimmomatic/{sample}/{sample}_R2_trim_paired.fastq.gz", 
        # reads where trimming entirely removed the mate
        r1_unpaired = "results/{prefix}/trimmomatic/{sample}/{sample}_R1_trim_unpaired.fastq.gz",
        r2_unpaired = "results/{prefix}/trimmomatic/{sample}/{sample}_R2_trim_unpaired.fastq.gz",
    params:
        adapter_filepath=config["adapter_file"],
        seed=config["seed_mismatches"],
        palindrome_clip=config["palindrome_clipthreshold"],
        simple_clip=config["simple_clipthreshold"],
        minadapterlength=config["minadapterlength"],
        keep_both_reads=config["keep_both_reads"],
        window_size=config["window_size"],
        window_size_quality=config["window_size_quality"],
        minlength=config["minlength"],
        headcrop_length=config["headcrop_length"],
        lead_trail_qual=config["lead_trail_qual"],
        #threads = config["ncores"],
    log:
        "logs/{prefix}/trimmomatic/{sample}/{sample}.log"
    threads: 8
    # threads: workflow.cores
    # This has been changed to specifically use 8 cores (the previous default number of cores used)
    # Before, this always allocated the maximum number of cores
    #conda:
    #    "envs/trimmomatic.yaml"
    singularity:
        "docker://staphb/trimmomatic:0.39"
    resources:
        mem_mb = 5000,
        runtime = 30
    #envmodules:
    #    "Bioinformatics",
    #    "trimmomatic"
    shell:
        """
        trimmomatic PE {input.r1} {input.r2} {output.r1} {output.r1_unpaired} {output.r2} {output.r2_unpaired} -threads {threads} \
        ILLUMINACLIP:{params.adapter_filepath}:{params.seed}:{params.palindrome_clip}:{params.simple_clip}:{params.minadapterlength}:{params.keep_both_reads} \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:{params.window_size}:{params.window_size_quality} MINLEN:{params.minlength} HEADCROP:{params.headcrop_length} &>{log}
        """

# old version
# rule trimmomatic_pe:
#     input:
#         r1 = lambda wildcards: expand(str(config["short_reads"] + "/" + f"{wildcards.sample}_R1.fastq.gz")),
#         r2 = lambda wildcards: expand(str(config["short_reads"] + "/" + f"{wildcards.sample}_R2.fastq.gz")),
        
#     output:
#         r1 = f"results/{{prefix}}/trimmomatic/{{barcode}}/{{sample}}/{{sample}}_R1_paired.fastq.gz",
#         r2 = f"results/{{prefix}}/trimmomatic/{{barcode}}/{{sample}}/{{sample}}_R2_paired.fastq.gz", 
#         # reads where trimming entirely removed the mate
#         r1_unpaired = f"results/{{prefix}}/trimmomatic/{{barcode}}/{{sample}}/{{sample}}_R1_unpaired.fastq.gz",
#         r2_unpaired = f"results/{{prefix}}/trimmomatic/{{barcode}}/{{sample}}/{{sample}}_R2_unpaired.fastq.gz",
#     params:
#         adapter_filepath=config["adapter_file"],
#         seed=config["seed_mismatches"],
#         palindrome_clip=config["palindrome_clipthreshold"],
#         simple_clip=config["simple_clipthreshold"],
#         minadapterlength=config["minadapterlength"],
#         keep_both_reads=config["keep_both_reads"],
#         window_size=config["window_size"],
#         window_size_quality=config["window_size_quality"],
#         minlength=config["minlength"],
#         headcrop_length=config["headcrop_length"],
#         threads = config["ncores"],
#     log:
#         "logs/{prefix}/trimmomatic/{barcode}/{sample}/{sample}.log"
#     # conda:
#     #     "envs/trimmomatic.yaml"
#     singularity:
#         "docker://staphb/trimmomatic:0.39"
#     shell:
#         "trimmomatic PE {input.r1} {input.r2} {output.r1} {output.r1_unpaired} {output.r2} {output.r2_unpaired} -threads {params.threads} ILLUMINACLIP:{params.adapter_filepath}:{params.seed}:{params.palindrome_clip}:{params.simple_clip}:{params.minadapterlength}:{params.keep_both_reads} SLIDINGWINDOW:{params.window_size}:{params.window_size_quality} MINLEN:{params.minlength} HEADCROP:{params.headcrop_length} &>{log}"

rule nanopore_flye_assembly:
    input:
        trimmed = "results/{prefix}/filtlong/{sample}/{sample}.trimmed.fastq.gz"
    output:
        assembly = "results/{prefix}/flye/{sample}/{sample}_flye.fasta",
    params:
        assembly_dir = "results/{prefix}/flye/{sample}/",
        size = config["genome_size"],
        #threads = config["threads"],
        #flye_options = config["flye_options"],
        sample = "{sample}",
    #log:
    #    "logs/{prefix}/flye/{sample}/{sample}_flye.log" # Flye has its own log
    singularity:
        "docker://staphb/flye:2.9.4"
    threads: 8
    resources:
        mem_mb = 8000,
        runtime = 30
    #envmodules:
    #    "Bioinformatics",
    #    "flye"
    shell:   
        """
        flye --nano-hq {input.trimmed} -g {params.size} -o {params.assembly_dir} -t {threads} --debug && 
        cp {params.assembly_dir}/assembly.fasta {params.assembly_dir}/{params.sample}_flye.fasta 
        """

# old version
# rule flye:
#     input:
#         trimmed = lambda wildcards: expand(str("results/" + f"{wildcards.prefix}" + "/filtlong/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}" + ".trimmed.fastq.gz")),
#     output:
#         assembly = f"results/{{prefix}}/flye/{{barcode}}/{{sample}}/{{sample}}_flye.fasta",
#     params:
#         assembly_dir = "results/{prefix}/flye/{barcode}/{sample}/",
#         size = config["genome_size"],
#         threads = config["threads"],
#         flye_options = config["flye_options"],
#         prefix = "{sample}",
#     log:
#         "logs/{prefix}/flye/{barcode}/{sample}/{sample}_flye.log"  
#     # conda:
#     #     "envs/flye.yaml"
#     singularity:
#         "docker://staphb/flye:2.9.3"
#     shell:
#         "flye --nano-hq {input.trimmed} -g {params.size} -o {params.assembly_dir} -t {params.threads} {params.flye_options} && cp {params.assembly_dir}/assembly.fasta {params.assembly_dir}/{params.prefix}_flye.fasta &>{log}"

# not needed for c auris
# rule flye_add_circ:
#     input:
#         flye_assembly = lambda wildcards: expand(str("results/" + f"{wildcards.prefix}" + "/flye/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}" + "_flye.fasta")),
#     output:
#         assembly = f"results/{{prefix}}/flye/{{barcode}}/{{sample}}/{{sample}}_flye_circ.fasta",
#     params:
#         assembly_dir = "results/{prefix}/flye/{barcode}/{sample}/",
#         size = config["genome_size"],
#         threads = config["threads"],
#         flye_options = config["flye_options"],
#         prefix = "{sample}",
#     log:
#         "logs/{prefix}/flye/{barcode}/{sample}/{sample}_flye.log"  
#     run:
#         shell("cp {params.assembly_dir}/{params.prefix}_flye.fasta {params.assembly_dir}/{params.prefix}_flye_circ.fasta")
#         assembly_info = pd.read_csv("%s/assembly_info.txt" % params.assembly_dir, sep='\t', header=0)
#         assembly_info["circular"] = np.where(assembly_info["circ."] == "Y", "true", "false")
#         flye_assembly = "%s" % input.flye_assembly
#         flye_assembly_circ = "%s" % output.assembly
#         for index, row in assembly_info.iterrows():
#             circular = "%s;circular=%s" % (row['#seq_name'], row['circular'])
#             shell("sed -i 's/\<%s\>/%s/g' %s" % (row['#seq_name'], circular, flye_assembly_circ))

rule nanopore_medaka:
    input:
        trimmed = "results/{prefix}/filtlong/{sample}/{sample}.trimmed.fastq.gz",
        flye_assembly = "results/{prefix}/flye/{sample}/{sample}_flye.fasta",
    output:
        medaka_out = "results/{prefix}/medaka/{sample}/{sample}_medaka.fasta",
    params:
        medaka_out_dir = "results/{prefix}/medaka/{sample}",
        #threads = config["threads"],
        sample = "{sample}",
    #log:
    #    "logs/{prefix}/medaka/{sample}/{sample}.log"
    singularity:
        "docker://staphb/medaka:2.0.1"
    threads: 8
    resources:
        mem_mb = 15000,
        runtime = 30
    #envmodules:
    #    "Bioinformatics",
    #    "medaka",
    #    "bcftools"
    shell:
        """
        medaka_consensus -i {input.trimmed} -d {input.flye_assembly} -o {params.medaka_out_dir} -t {threads} &&
        cp {params.medaka_out_dir}/consensus.fasta {params.medaka_out_dir}/{params.sample}_medaka.fasta 
        """ 

# old version
# rule medaka:
#     input:
#         trimmed = lambda wildcards: expand(str("results/" + f"{wildcards.prefix}" + "/filtlong/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}" + ".trimmed.fastq.gz")),
#         flye_assembly = lambda wildcards: expand(str("results/" + f"{wildcards.prefix}" + "/flye/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}" + "_flye_circ.fasta")),
#     output:
#         medaka_out = f"results/{{prefix}}/medaka/{{barcode}}/{{sample}}/{{sample}}_medaka.fasta"
#     params:
#         medaka_out_dir = "results/{prefix}/medaka/{barcode}/{sample}",
#         threads = config["threads"],
#         prefix = f"{{sample}}",
#         model = config["medaka_model"],
#     #conda:
#     #    "envs/medaka.yaml"
#     log:
#         "logs/{prefix}/medaka/{barcode}/{sample}/{sample}.log"
#     shell:
#         "module load Bioinformatics medaka bcftools/1.12-g4b275e && medaka_consensus -i {input.trimmed} -d {input.flye_assembly} -o {params.medaka_out_dir} -t {params.threads} -m {params.model} && cp {params.medaka_out_dir}/consensus.fasta {params.medaka_out_dir}/{params.prefix}_medaka.fasta &>{log}"

rule polypolish_bwa_align:
    input:
        r1 = "results/{prefix}/trimmomatic/{sample}/{sample}_R1_trim_paired.fastq.gz",
        r2 = "results/{prefix}/trimmomatic/{sample}/{sample}_R2_trim_paired.fastq.gz", 
        r1_unpaired = "results/{prefix}/trimmomatic/{sample}/{sample}_R1_trim_unpaired.fastq.gz",
        r2_unpaired = "results/{prefix}/trimmomatic/{sample}/{sample}_R2_trim_unpaired.fastq.gz",
        medaka_assembly = "results/{prefix}/medaka/{sample}/{sample}_medaka.fasta"
        #r1 = lambda wildcards: expand(f"results/{wildcards.prefix}/trimmomatic/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}_R1_paired.fastq.gz"),
        #r2 = lambda wildcards: expand(f"results/{wildcards.prefix}/trimmomatic/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}_R2_paired.fastq.gz"),
        #r1_unpaired = lambda wildcards: expand(f"results/{wildcards.prefix}/trimmomatic/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}_R1_unpaired.fastq.gz"),
        #r2_unpaired = lambda wildcards: expand(f"results/{wildcards.prefix}/trimmomatic/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}_R2_unpaired.fastq.gz"),
        #medaka_assembly = f"results/{{prefix}}/medaka/{{barcode}}/{{sample}}/{{sample}}_medaka.fasta"
    output:
        samout_r1 = "results/{prefix}/polypolish/{sample}/{sample}_r1.sam",
        samout_r2 = "results/{prefix}/polypolish/{sample}/{sample}_r2.sam",
        samout_r1unpair = "results/{prefix}/polypolish/{sample}/{sample}_r1unpair.sam",
        samout_r2unpair = "results/{prefix}/polypolish/{sample}/{sample}_r2unpair.sam",
        # samout_1 = "results/{prefix}/polypolish/{barcode}/{sample}/{sample}_1.sam",
        # samout_2 = "results/{prefix}/polypolish/{barcode}/{sample}/{sample}_2.sam",
        # samout_3 = "results/{prefix}/polypolish/{barcode}/{sample}/{sample}_3.sam",
        # samout_4 = "results/{prefix}/polypolish/{barcode}/{sample}/{sample}_4.sam",
    # params:
    #     threads = config["ncores"],
    # conda:
    #     "envs/bwa.yaml"
    singularity:
        "docker://staphb/bwa:0.7.17"
    threads: 8
    resources:
        mem_mb = 15000,
        runtime = 600
    shell:
        """
        bwa index {input.medaka_assembly} && \
        bwa mem -t {threads} -a {input.medaka_assembly} {input.r1} > {output.samout_r1} && \
        bwa mem -t {threads} -a {input.medaka_assembly} {input.r2} > {output.samout_r2} && \
        bwa mem -t {threads} -a {input.medaka_assembly} {input.r1_unpaired} > {output.samout_r1unpair} && \
        bwa mem -t {threads} -a {input.medaka_assembly} {input.r2_unpaired} > {output.samout_r2unpair}
        """
        # "bwa index {input.medaka_assembly} && bwa mem -t12 -a {input.medaka_assembly} {input.r1} > {output.samout_1} && bwa mem -t12 -a {input.medaka_assembly} {input.r2} > {output.samout_2} && bwa mem -t12 -a {input.medaka_assembly} {input.r1_unpaired} > {output.samout_3} && bwa mem -t12 -a {input.medaka_assembly} {input.r2_unpaired} > {output.samout_4}"

# Deprecated: This rule was used for pilon 
#rule bwaalign_polypolish:
#    input:
#        r1 = lambda wildcards: expand(f"results/{wildcards.prefix}/trimmomatic/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}_R1_paired.fastq.gz"),
#        r2 = lambda wildcards: expand(f"results/{wildcards.prefix}/trimmomatic/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}_R2_paired.fastq.gz"),
        #r1_unpaired = lambda wildcards: expand(f"results/{wildcards.prefix}/trimmomatic/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}_R1_unpaired.fastq.gz"),
        #r2_unpaired = lambda wildcards: expand(f"results/{wildcards.prefix}/trimmomatic/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}_R2_unpaired.fastq.gz"),
#        polypolish_assembly = f"results/{{prefix}}/polypolish/{{barcode}}/{{sample}}/{{sample}}_flye_medaka_polypolish.fasta"
#    output:
#       samout = f"results/{{prefix}}/pilon/{{barcode}}/{{sample}}/{{sample}}.sam",
#        bamout = f"results/{{prefix}}/pilon/{{barcode}}/{{sample}}/{{sample}}.bam",
#        bamout_sorted = f"results/{{prefix}}/pilon/{{barcode}}/{{sample}}/{{sample}}_sorted.bam",
#    params:
#        threads = config["ncores"],
    # conda:
    #     "envs/bwa.yaml"
#    singularity:
#        "docker://staphb/bwa:0.7.17"
#    shell:
#        "bwa index {input.polypolish_assembly} && bwa mem -t12 {input.polypolish_assembly} {input.r1} {input.r2} > {output.samout} && samtools view -Sb {output.samout} > {output.bamout} && samtools sort -o {output.bamout_sorted} {output.bamout} && samtools index {output.bamout_sorted}"

# rule bwaalign_unicycler:
#     input:
#         r1 = "results/{prefix}/trimmomatic/{sample}/{sample}_R1_paired.fastq.gz",
#         r2 = "results/{prefix}/trimmomatic/{sample}/{sample}_R2_paired.fastq.gz",
#         r1_unpaired = "results/{prefix}/trimmomatic/{sample}/{sample}_R1_unpaired.fastq.gz",
#         r2_unpaired = "results/{prefix}/trimmomatic/{sample}/{sample}_R2_unpaired.fastq.gz",
#         unicycler_assembly = "results/{prefix}/unicycler/{sample}/{sample}_unicycler.fasta"
#         # r1 = lambda wildcards: expand(f"results/{wildcards.prefix}/trimmomatic/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}_R1_paired.fastq.gz"),
#         # r2 = lambda wildcards: expand(f"results/{wildcards.prefix}/trimmomatic/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}_R2_paired.fastq.gz"),
#         # r1_unpaired = lambda wildcards: expand(f"results/{wildcards.prefix}/trimmomatic/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}_R1_unpaired.fastq.gz"),
#         # r2_unpaired = lambda wildcards: expand(f"results/{wildcards.prefix}/trimmomatic/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}_R2_unpaired.fastq.gz"),
#         # unicycler_assembly = f"results/{{prefix}}/unicycler/{{barcode}}/{{sample}}/{{sample}}_unicycler.fasta"
#     output:
#         samout_r1 = "results/{prefix}/polypolish_unicycler/{sample}/{sample}_1.sam",
#         samout_r2 = "results/{prefix}/polypolish_unicycler/{sample}/{sample}_2.sam",
#         samout_r1unpair = "results/{prefix}/polypolish_unicycler/{sample}/{sample}_3.sam",
#         samout_r2unpair = "results/{prefix}/polypolish_unicycler/{sample}/{sample}_4.sam",
#         # samout_1 = f"results/{{prefix}}/polypolish_unicycler/{{barcode}}/{{sample}}/{{sample}}_1.sam",
#         # samout_2 = f"results/{{prefix}}/polypolish_unicycler/{{barcode}}/{{sample}}/{{sample}}_2.sam",
#         # samout_3 = f"results/{{prefix}}/polypolish_unicycler/{{barcode}}/{{sample}}/{{sample}}_3.sam",
#         # samout_4 = f"results/{{prefix}}/polypolish_unicycler/{{barcode}}/{{sample}}/{{sample}}_4.sam",
#     # params:
#     #     threads = config["ncores"],
#     # conda:
#     #     "envs/bwa.yaml"
#     singularity:
#         "docker://staphb/bwa:0.7.17"
#     threads: 8
#     shell:
#         """
#         bwa index {input.medaka_assembly} && \
#         bwa mem -t {threads} -a {input.medaka_assembly} {input.r1} > {output.samout_r1} && \
#         bwa mem -t {threads} -a {input.medaka_assembly} {input.r2} > {output.samout_r2} && \
#         bwa mem -t {threads} -a {input.medaka_assembly} {input.r1_unpaired} > {output.samout_r1unpair} && \
#         bwa mem -t {threads} -a {input.medaka_assembly} {input.r2_unpaired} > {output.samout_r2unpair}
#         """
#         # "bwa index {input.unicycler_assembly} && bwa mem -t12 -a {input.unicycler_assembly} {input.r1} > {output.samout_1} && bwa mem -t12 -a {input.unicycler_assembly} {input.r2} > {output.samout_2} && bwa mem -t12 -a {input.unicycler_assembly} {input.r1_unpaired} > {output.samout_3} && bwa mem -t12 -a {input.unicycler_assembly} {input.r2_unpaired} > {output.samout_4}"

rule polypolish_polish:
    input:
        r1 = "results/{prefix}/trimmomatic/{sample}/{sample}_R1_trim_paired.fastq.gz",
        r2 = "results/{prefix}/trimmomatic/{sample}/{sample}_R2_trim_paired.fastq.gz", 
        medaka_assembly = "results/{prefix}/medaka/{sample}/{sample}" + "_medaka.fasta",
        samout_1 = "results/{prefix}/polypolish/{sample}/{sample}_r1.sam",
        samout_2 = "results/{prefix}/polypolish/{sample}/{sample}_r2.sam",
        # r1 = lambda wildcards: expand(f"results/{wildcards.prefix}/trimmomatic/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}_R1_paired.fastq.gz"),
        # r2 = lambda wildcards: expand(f"results/{wildcards.prefix}/trimmomatic/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}_R2_paired.fastq.gz"),
        # medaka_assembly = lambda wildcards: expand(f"results/{wildcards.prefix}/medaka/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}" + "_medaka.fasta"),
        # samout_1 = lambda wildcards: expand(f"results/{wildcards.prefix}/polypolish/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}_1.sam"),
        # samout_2 = lambda wildcards: expand(f"results/{wildcards.prefix}/polypolish/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}_2.sam"),
    output:
        filtersam1 = "results/{prefix}/polypolish/{sample}/{sample}_filtered_1.sam",
        filtersam2 = "results/{prefix}/polypolish/{sample}/{sample}_filtered_2.sam",
        flye_medaka_polypolish = "results/{prefix}/polypolish/{sample}/{sample}_polypolish.fasta",
    # params:
    #     threads = config["ncores"],
    # conda:
    #     "envs/polypolish.yaml"
    singularity:
        "docker://staphb/polypolish:0.6.0"
    threads: 4
    resources:
        mem_mb = 15000,
        runtime = 600
    shell:
        """
        polypolish filter --in1 {input.samout_1} --in2 {input.samout_2} --out1 {output.filtersam1} --out2 {output.filtersam2} && \
        polypolish polish {input.medaka_assembly} {output.filtersam1} {output.filtersam2} > {output.flye_medaka_polypolish}
        """
        # "polypolish filter --in1 {input.samout_1} --in2 {input.samout_2} --out1 {output.filtersam1} --out2 {output.filtersam2} && polypolish polish {input.medaka_assembly} {output.filtersam1} {output.filtersam2} > {output.flye_medaka_polypolish}"

# unicycler is for bacterial genomes
# rule unicycler:
#     input:
#         r1 = "results/{prefix}/trimmomatic/{sample}/{sample}_R1_paired.fastq.gz",
#         r2 = "results/{prefix}/trimmomatic/{sample}/{sample}_R2_paired.fastq.gz",
#         trimmed_long = "results/{prefix}/filtlong/{sample}/{sample}.trimmed.fastq.gz",
#         r1_unpaired = "results/{prefix}/trimmomatic/{sample}/{sample}_R1_unpaired.fastq.gz",
#         r2_unpaired = "results/{prefix}/trimmomatic/{sample}/{sample}_R2_unpaired.fastq.gz",
#         # r1 = lambda wildcards: expand(f"results/{wildcards.prefix}/trimmomatic/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}_R1_paired.fastq.gz"),
#         # r2 = lambda wildcards: expand(f"results/{wildcards.prefix}/trimmomatic/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}_R2_paired.fastq.gz"),
#         # trimmed_long = "results/{prefix}/filtlong/{barcode}/{sample}/{sample}.trimmed.fastq.gz",
#         # r1_unpaired = lambda wildcards: expand(f"results/{wildcards.prefix}/trimmomatic/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}_R1_unpaired.fastq.gz"),
#         # r2_unpaired = lambda wildcards: expand(f"results/{wildcards.prefix}/trimmomatic/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}_R2_unpaired.fastq.gz"),
#     output:
#         unicycler_assembly = "results/{prefix}/unicycler/{sample}/{sample}_unicycler.fasta",
#     params:
#         unicycler_assembly_out = "results/{prefix}/unicycler/{sample}",
#         #threads = config["ncores"],
#         sample = "{sample}",
#     # conda:
#     #     "envs/unicycler.yaml"
#     singularity:
#         "docker://staphb/unicycler:0.5.0"
#     threads: 8
#     shell:
#         """
#         unicycler -1 {input.r1} -2 {input.r2} -s {input.r1_unpaired} -s {input.r2_unpaired} -l {input.trimmed_long} -o {params.unicycler_assembly_out} && \
#         cp {params.unicycler_assembly_out}/assembly.fasta {params.unicycler_assembly_out}/{params.sample}_unicycler.fasta && \
#         sed -i 's/length=.*depth=.*circular/circular/g' {params.unicycler_assembly_out}/{params.sample}_unicycler.fasta && \
#         sed -i 's/ circular/;circular/g' {params.unicycler_assembly_out}/{params.sample}_unicycler.fasta
#         """
#         # "unicycler -1 {input.r1} -2 {input.r2} -s {input.r1_unpaired} -s {input.r2_unpaired} -l {input.trimmed_long} -o {params.unicycler_assembly_out} && cp {params.unicycler_assembly_out}/assembly.fasta {params.unicycler_assembly_out}/{params.prefix}_unicycler.fasta && sed -i 's/length=.*depth=.*circular/circular/g' {params.unicycler_assembly_out}/{params.prefix}_unicycler.fasta && sed -i 's/ circular/;circular/g' {params.unicycler_assembly_out}/{params.prefix}_unicycler.fasta"

# unicycler is for bacterial genomes
# rule polypolish_unicycler:
#     input:
#         r1 = "results/{prefix}/trimmomatic/{sample}/{sample}_R1_paired.fastq.gz",
#         r2 = "results/{prefix}/trimmomatic/{sample}/{sample}_R2_paired.fastq.gz",
#         unicycler_assembly = "results/{prefix}/unicycler/{sample}/{sample}_unicycler.fasta",
#         samout_1 = "results/{prefix}/polypolish_unicycler/{sample}/{sample}_1.sam",
#         samout_2 = "results/{prefix}/polypolish_unicycler/{sample}/{sample}_2.sam",
#         # r1 = "results/{wildcards.prefix}/trimmomatic/{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}_R1_paired.fastq.gz",
#         # r2 = "results/{wildcards.prefix}/trimmomatic/{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}_R2_paired.fastq.gz",
#         # unicycler_assembly = "results/{prefix}/unicycler/{sample}/{sample}_unicycler.fasta",
#         # samout_1 = "results/{wildcards.prefix}/polypolish_unicycler/{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}_1.sam",
#         # samout_2 = "results/{wildcards.prefix}/polypolish_unicycler/{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}_2.sam",
#     output:
#         filtersam1 = "results/{prefix}/polypolish_unicycler/{sample}/{sample}_filtered_1.sam",
#         filtersam2 = "results/{prefix}/polypolish_unicycler/{sample}/{sample}_filtered_2.sam",
#         unicycler_polypolish = "results/{prefix}/polypolish_unicycler/{sample}/{sample}_unicycler_polypolish.fasta",
#     # params:
#     #     threads = config["ncores"],
#     # conda:
#     #     "envs/polypolish.yaml"
#     singularity:
#         "docker://staphb/polypolish:0.6.0"
#     threads: 8
#     shell:
#         """
#         polypolish filter --in1 {input.samout_1} --in2 {input.samout_2} --out1 {output.filtersam1} --out2 {output.filtersam2} && \
#         polypolish polish {input.unicycler_assembly} {output.filtersam1} {output.filtersam2} > {output.unicycler_polypolish}
#         """
#         #"polypolish filter --in1 {input.samout_1} --in2 {input.samout_2} --out1 {output.filtersam1} --out2 {output.filtersam2} && polypolish polish {input.unicycler_assembly} {output.filtersam1} {output.filtersam2} > {output.unicycler_polypolish}"

# Deprecated: This rule was used for pilon
#rule bwaalign_polypolish_unicycler:
#    input:
#        r1 = lambda wildcards: expand(f"results/{wildcards.prefix}/trimmomatic/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}_R1_paired.fastq.gz"),
#        r2 = lambda wildcards: expand(f"results/{wildcards.prefix}/trimmomatic/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}_R2_paired.fastq.gz"),
#        polypolish_unicycler_assembly = f"results/{{prefix}}/polypolish_unicycler/{{barcode}}/{{sample}}/{{sample}}_unicycler_polypolish.fasta"
#    output:
#        samout = f"results/{{prefix}}/pilon_unicycler/{{barcode}}/{{sample}}/{{sample}}.sam",
#        bamout = f"results/{{prefix}}/pilon_unicycler/{{barcode}}/{{sample}}/{{sample}}.bam",
#        bamout_sorted = f"results/{{prefix}}/pilon_unicycler/{{barcode}}/{{sample}}/{{sample}}_sorted.bam",
#    params:
#        threads = config["ncores"],
    # conda:
    #     "envs/bwa.yaml"
#    singularity:
#        "docker://staphb/bwa:0.7.17"
#    shell:
#        "bwa index {input.polypolish_unicycler_assembly} && bwa mem -t12 {input.polypolish_unicycler_assembly} {input.r1} {input.r2} > {output.samout} && samtools view -Sb {output.samout} > {output.bamout} && samtools sort -o {output.bamout_sorted} {output.bamout} && samtools index {output.bamout_sorted}"

# rule prokka:
#     input:
#         unicycler_polypolish = f"results/{{prefix}}/polypolish_unicycler/{{barcode}}/{{sample}}/{{sample}}_unicycler_polypolish.fasta",
#         unicycler_assembly = f"results/{{prefix}}/unicycler/{{barcode}}/{{sample}}/{{sample}}_unicycler.fasta",
#         flye_medaka_polypolish = f"results/{{prefix}}/polypolish/{{barcode}}/{{sample}}/{{sample}}_flye_medaka_polypolish.fasta",
#         flye_assembly = f"results/{{prefix}}/flye/{{barcode}}/{{sample}}/{{sample}}_flye_circ.fasta",
#     output:
#         unicycler_annotation = f"results/{{prefix}}/prokka/{{barcode}}/{{sample}}/{{sample}}_unicycler.gff",
#     params:
#         threads = config["ncores"],
#         prefix = "{sample}",
#         options = config["prokka_options"],
#         prokka_dir = directory("results/{prefix}/prokka/{barcode}/{sample}/"),
#     # conda:
#     #     "envs/prokka.yaml"
#     singularity:
#         "docker://staphb/prokka:1.14.6"
#     shell:
#         "prokka {params.options} --strain {params.prefix} -outdir {params.prokka_dir} -prefix {params.prefix}_unicycler {input.unicycler_assembly} && prokka {params.options} --strain {params.prefix} -outdir {params.prokka_dir} -prefix {params.prefix}_flye {input.flye_assembly} && prokka {params.options} --strain {params.prefix} -outdir {params.prokka_dir} -prefix {params.prefix}_flye_medaka_polypolish {input.flye_medaka_polypolish} && prokka {params.options} --strain {params.prefix} -outdir {params.prokka_dir} -prefix {params.prefix}_unicycler_polypolish {input.unicycler_polypolish}"

rule quast:
    input:
        #unicycler_polypolish = "results/{prefix}/polypolish_unicycler/{sample}/{sample}_unicycler_polypolish.fasta",
        #unicycler_assembly = "results/{prefix}/unicycler/{sample}/{sample}_unicycler.fasta",
        flye_medaka_polypolish = "results/{prefix}/polypolish/{sample}/{sample}_polypolish.fasta",
        flye_assembly = "results/{prefix}/flye/{sample}/{sample}_flye.fasta",
    output:
        quast_out = "results/{prefix}/quast/{sample}/{sample}_flye_medaka_polypolish/report.txt",
    params:
        # threads = config["ncores"],
        quast_dir = directory("results/{prefix}/quast/{sample}/{sample}"),
    # conda:
    #     "envs/quast.yaml"
    singularity:
        "docker://staphb/quast:5.0.2"
    #threads: 8
    shell:
        """
        quast.py {input.flye_medaka_polypolish} -o {params.quast_dir}_flye_medaka_polypolish && \
        quast.py {input.flye_assembly} -o {params.quast_dir}_flye
        """
        #"quast.py {input.unicycler_assembly} -o {params.quast_dir}_unicycler && quast.py {input.unicycler_polypolish} -o {params.quast_dir}_unicycler_polypolish && quast.py {input.flye_medaka_polypolish} -o {params.quast_dir}_flye_medaka_polypolish && quast.py {input.flye_assembly} -o {params.quast_dir}_flye"

# rule busco:
#     input:
#         unicycler_polypolish = "results/{prefix}/polypolish_unicycler/{sample}/{sample}_unicycler_polypolish.fasta",
#         unicycler = "results/{prefix}/unicycler/{sample}/{sample}_unicycler.fasta",
#         flye_medaka_polypolish = "results/{prefix}/polypolish/{sample}/{sample}_flye_medaka_polypolish.fasta",
#         flye_assembly = "results/{prefix}/flye/{sample}/{sample}_flye.fasta",
#     output:
#         busco_out = "results/{prefix}/busco/{sample}/{sample}.unicycler/busco_unicycler.txt",
#     params:
#         busco_outpath = "results/{prefix}/busco/{sample}/{sample}",
#         unicycler_busco_out = "short_summary.specific.bacteria_odb10.{sample}.unicycler.txt",
#         unicycler_polypolish_busco_out = "short_summary.specific.bacteria_odb10.{sample}.unicycler_polypolish.txt",
#         flye_medaka_polypolish_busco_out = "short_summary.specific.bacteria_odb10.{sample}.flye_medaka_polypolish.txt",
#         medaka_assembly_busco_out = "short_summary.specific.bacteria_odb10.{sample}.medaka_assembly.txt",
#         flye_assembly_busco_out = "short_summary.specific.bacteria_odb10.{sample}.flye_assembly.txt",
#         flye_medaka_polypolish_pilon_busco_out = "short_summary.specific.bacteria_odb10.{sample}.flye_medaka_polypolish_pilon.txt",
#         unicycler_polypolish_pilon_busco_out = "short_summary.specific.bacteria_odb10.{sample}.unicycler_polypolish_pilon.txt",
#         #threads = config["ncores"],
#     # conda:
#     #     "envs/busco.yaml"
#     singularity:
#         "docker://ezlabgva/busco:v5.7.0_cv1"
#     threads: 8
#     shell:
#         """
#         busco -f -i {input.unicycler} -m genome -l bacteria_odb10 -o {params.busco_outpath}.unicycler && \
#         cp {params.busco_outpath}.unicycler/{params.unicycler_busco_out} {params.busco_outpath}.unicycler/busco_unicycler.txt && \
#         busco -f -i {input.unicycler_polypolish} -m genome -l bacteria_odb10 -o {params.busco_outpath}.unicycler_polypolish && \
#         cp {params.busco_outpath}.unicycler_polypolish/{params.unicycler_polypolish_busco_out} {params.busco_outpath}.unicycler_polypolish/busco_unicycler_polypolish.txt && \
#         busco -f -i {input.flye_medaka_polypolish} -m genome -l bacteria_odb10 -o {params.busco_outpath}.flye_medaka_polypolish && \
#         cp {params.busco_outpath}.flye_medaka_polypolish/{params.flye_medaka_polypolish_busco_out} {params.busco_outpath}.flye_medaka_polypolish/busco_flye_medaka_polypolish.txt && \
#         busco -f -i {input.flye_assembly} -m genome -l bacteria_odb10 -o {params.busco_outpath}.flye_assembly && \
#         cp {params.busco_outpath}.flye_assembly/{params.flye_assembly_busco_out} {params.busco_outpath}.flye_assembly/busco_flye_assembly.txt
#         """
#         #"busco -f -i {input.unicycler} -m genome -l bacteria_odb10 -o {params.busco_outpath}.unicycler && cp {params.busco_outpath}.unicycler/{params.unicycler_busco_out} {params.busco_outpath}.unicycler/busco_unicycler.txt && busco -f -i {input.unicycler_polypolish} -m genome -l bacteria_odb10 -o {params.busco_outpath}.unicycler_polypolish && cp {params.busco_outpath}.unicycler_polypolish/{params.unicycler_polypolish_busco_out} {params.busco_outpath}.unicycler_polypolish/busco_unicycler_polypolish.txt && busco -f -i {input.flye_medaka_polypolish} -m genome -l bacteria_odb10 -o {params.busco_outpath}.flye_medaka_polypolish && cp {params.busco_outpath}.flye_medaka_polypolish/{params.flye_medaka_polypolish_busco_out} {params.busco_outpath}.flye_medaka_polypolish/busco_flye_medaka_polypolish.txt && busco -f -i {input.flye_assembly} -m genome -l bacteria_odb10 -o {params.busco_outpath}.flye_assembly && cp {params.busco_outpath}.flye_assembly/{params.flye_assembly_busco_out} {params.busco_outpath}.flye_assembly/busco_flye_assembly.txt"




# copied from funQCD

rule funannotate_sort:
    input:
        flye_medaka_polypolish = "results/{prefix}/polypolish/{sample}/{sample}_polypolish.fasta",
    output:
        sorted_assembly = "results/{prefix}/funannotate/{sample}/{sample}_sorted.fa"
    params:
        out_dir = "results/{prefix}/funannotate/{sample}/",
        sample = "{sample}"
    singularity:
        "docker://nextgenusfs/funannotate:v1.8.17"
    resources:
        mem_mb = 5000,
        runtime = 60
    shell:
        """
        funannotate clean -i {input.flye_medaka_polypolish} -o {params.out_dir}{params.sample}_cleaned.fa
        funannotate sort -i {params.out_dir}{params.sample}_cleaned.fa -o {params.out_dir}{params.sample}_sorted.fa --minlen 0
        """

rule repeatmasker:
    input:
        sorted_assembly = "results/{prefix}/funannotate/{sample}/{sample}_sorted.fa"
    output:
        masked_assembly = "results/{prefix}/repeatmasker/{sample}/{sample}_masked.fa"
    params:
        out_dir = "results/{prefix}/repeatmasker/{sample}/",
        repeat_lib = config["funqcd_lib"] + "repeat_libraries/fungi_b8441/b8441_fungi_repeatlib.fa",
        prefix = "{prefix}",
        sample = "{sample}",
    threads: 8
    resources:
        mem_mb = 5000,
        runtime = 120,
    singularity:
        "docker://dfam/tetools:1.89.2"
    shell:
        """
        RepeatMasker -xsmall -dir {params.out_dir} -lib {params.repeat_lib} \
        results/{params.prefix}/funannotate/{params.sample}/{params.sample}_sorted.fa -pa {threads}
        mv {params.out_dir}/{params.sample}_sorted.fa.masked {params.out_dir}/{params.sample}_masked.fa
        """

# this requires the RNA-seq data from Teresa, and assumes that the files are in a specific format
# specifically, that read 1 contains 'R1' in the file name, that all files are in fastq.gz format and in RF order for stranded RNA-seq, and that they can be ordered alphabetically
# This is also currently hard-coded to use 30G of memory
rule funannotate_train:
    input:
        masked_assembly = "results/{prefix}/repeatmasker/{sample}/{sample}_masked.fa"
    output:
        funannotate_training_dir = directory("results/{prefix}/funannotate/{sample}/training/"),
        funannotate_training_rna_bam = "results/{prefix}/funannotate/{sample}/training/funannotate_train.coordSorted.bam",
        funannotate_training_pasa_gff = "results/{prefix}/funannotate/{sample}/training/funannotate_train.pasa.gff3",
        funannotate_training_stringtie = "results/{prefix}/funannotate/{sample}/training/funannotate_train.stringtie.gtf",
        funannotate_training_transc_align = "results/{prefix}/funannotate/{sample}/training/funannotate_train.transcripts.gff3",
    params:
        out_dir = "results/{prefix}/funannotate/{sample}/",
        sample = "{sample}",
        rna_data_r1 = ' '.join(sorted([config["funqcd_lib"] + 'rna_seq_data/' + x for x in os.listdir(config["funqcd_lib"] + 'rna_seq_data/') if '_R1_' in x and 'fastq.gz' in x])),
        rna_data_r2 = ' '.join(sorted([config["funqcd_lib"] + 'rna_seq_data/' + x for x in os.listdir(config["funqcd_lib"] + 'rna_seq_data/') if '_R2_' in x and 'fastq.gz' in x])),
        mem_g = "30G",
        test = vars(workflow.resources)
    # threads: workflow.cores
    threads: 8
    resources:
        mem_mb = 32000,
        runtime = 700
    singularity:
        "docker://nextgenusfs/funannotate:v1.8.17"
    shell:
        """
        funannotate train --input {input.masked_assembly} --out {params.out_dir} \
        --left {params.rna_data_r1} --right {params.rna_data_r2} --stranded RF \
        --jaccard_clip --species "Candida auris" --isolate {params.sample} --cpus {threads} --memory {params.mem_g}
        """

# This should automatically detect the four training files generated previously, even without explicit input
# All steps should run with 'pasa' or 'rna-bam' under Training-Methods. Nothing should run with 'busco'.
rule funannotate_predict:
    input:
        masked_assembly = "results/{prefix}/repeatmasker/{sample}/{sample}_masked.fa",
        funannotate_training_rna_bam = "results/{prefix}/funannotate/{sample}/training/funannotate_train.coordSorted.bam",
        busco_db = config["funqcd_lib"] + "busco/lineages/saccharomycetes_odb10/dataset.cfg"
    output:
        funannotate_predict_proteins = "results/{prefix}/funannotate/{sample}/predict_results/{sample}.proteins.fa",
        funannotate_predict_out = directory("results/{prefix}/funannotate/{sample}/predict_results/")
    params:
        out_dir = "results/{prefix}/funannotate/{sample}/",
        sample = "{sample}",
        genemark_path = config["funqcd_lib"] + "genemark/gmes_linux_64_4/",
    threads: 8
    resources:
        mem_mb = 15000,
        runtime = 500
    singularity:
        "docker://nextgenusfs/funannotate:v1.8.17"
    shell:
        # removed --species "Candida auris" --isolate {params.sample}
        # consider changing output names and re-adding
        """
        funannotate predict --input {input.masked_assembly} --out {params.out_dir} \
        --species {params.sample} --force \
        --busco_seed_species candida_albicans --busco_db saccharomycetes_odb10 --cpus {threads} \
        --GENEMARK_PATH {params.genemark_path}
        """

rule funannotate_update:
    input:
        funannotate_predict_proteins = "results/{prefix}/funannotate/{sample}/predict_results/{sample}.proteins.fa"
    output:
        funannotate_update_proteins = "results/{prefix}/funannotate/{sample}/update_results/{sample}.proteins.fa",
        funannotate_update_out = directory("results/{prefix}/funannotate/{sample}/update_results/")
    params:
        out_dir = "results/{prefix}/funannotate/{sample}/",
        sample = "{sample}"
    threads: 8
    resources:
        mem_mb = 15000,
        runtime = 600
    singularity:
        "docker://nextgenusfs/funannotate:v1.8.17"
    shell:
        """
        funannotate update --input {params.out_dir} --cpus {threads}
        """


rule interproscan:
    input:
        interproscan_data = config["funqcd_lib"] + "interproscan_data/data/antifam/",
        test_log = config["funqcd_lib"] + "interproscan_data/test_log.txt",
        funannotate_predict_proteins = "results/{prefix}/funannotate/{sample}/predict_results/{sample}.proteins.fa"
    output:
        interproscan_out = "results/{prefix}/funannotate/{sample}/interproscan/{sample}.proteins.fa.xml"
    singularity:
        "docker://interpro/interproscan:5.71-102.0"
    params:
        out_dir = "results/{prefix}/funannotate/{sample}/interproscan/",
        sample = "{sample}",
        # cpus = config["ncores"]
    threads: 8
    resources:
        mem_mb = 15000,
        runtime = 600
    shell:
        """
        bash /opt/interproscan/interproscan.sh --input {input.funannotate_predict_proteins} --output-dir {params.out_dir} \
        --disable-precalc --cpu {threads}
        """


rule eggnog:
    input:
        eggnog_data = config["funqcd_lib"] + "eggnog_data/eggnog.db",
        funannotate_predict_proteins = "results/{prefix}/funannotate/{sample}/predict_results/{sample}.proteins.fa"
    output:
        eggnog_out = "results/{prefix}/funannotate/{sample}/eggnog/{sample}.emapper.annotations"
    singularity:
        "docker://nanozoo/eggnog-mapper:2.1.9--4f2b6c0"
    params:
        eggnog_data_dir = config["funqcd_lib"] + "eggnog_data/",
        out_dir = "results/{prefix}/funannotate/{sample}/eggnog/",
        sample = "{sample}"
    threads: 8
    resources:
        mem_mb = 15000,
        runtime = 600
    shell:
        """
        emapper.py -i {input.funannotate_predict_proteins} --itype proteins --data_dir {params.eggnog_data_dir} -m diamond \
        --output {params.sample} --output_dir {params.out_dir} --cpu {threads} --override
        """

rule funannotate_annotate:
    input:
        funannotate_predict_out = "results/{prefix}/funannotate/{sample}/update_results/",
        interproscan_out = "results/{prefix}/funannotate/{sample}/interproscan/{sample}.proteins.fa.xml",
        eggnog_out = "results/{prefix}/funannotate/{sample}/eggnog/{sample}.emapper.annotations",
        busco_db = config["funqcd_lib"] + "busco/lineages/saccharomycetes_odb10/dataset.cfg"
    output:
        funannotate_annotate_proteins = "results/{prefix}/funannotate/{sample}/annotate_results/{sample}.proteins.fa",
        funannotate_annotate_assembly = "results/{prefix}/funannotate/{sample}/annotate_results/{sample}.scaffolds.fa",
    params:
        out_dir = "results/{prefix}/funannotate/{sample}/",
        sample = "{sample}",
        # cpus = config["ncores"]
    threads: 8
    resources:
        mem_mb = 3000,
        runtime = 80
    singularity:
        "docker://nextgenusfs/funannotate:v1.8.17"
    shell:
        """
        funannotate annotate -i {input.funannotate_predict_out} -o {params.out_dir} --cpus {threads} \
        --iprscan {input.interproscan_out} --eggnog {input.eggnog_out} --busco_db saccharomycetes_odb10
        """

# The line 'rm -rf RM_*' removes the directories that RepeatMasker generates in the working directory
rule busco_final:
    input:
        funannotate_annotate_proteins = expand("results/{prefix}/funannotate/{sample}/annotate_results/{sample}.proteins.fa", prefix = PREFIX, sample = SAMPLE),
        funannotate_annotate_nucleotides = expand("results/{prefix}/funannotate/{sample}/annotate_results/{sample}.scaffolds.fa", prefix = PREFIX, sample = SAMPLE),       
    output:
        busco_out_p = "results/{prefix}/busco/busco_output_prot/batch_summary.txt",
        busco_out_n = "results/{prefix}/busco/busco_output_nucl/batch_summary.txt",
    params:
        prefix = "{prefix}",
        busco_db = config["funqcd_lib"] + "busco/"
    threads: 8
    resources:
        mem_mb = 15000,
        runtime = 45
    singularity:
        "docker://ezlabgva/busco:v5.7.0_cv1"
    shell:
        """
        mkdir -p results/{params.prefix}/busco/input/prot
	mkdir results/{params.prefix}/busco/input/nucl
        cp results/{params.prefix}/funannotate/*/annotate_results/*.proteins.fa results/{params.prefix}/busco/input/prot
        cp results/{params.prefix}/funannotate/*/annotate_results/*.scaffolds.fa results/{params.prefix}/busco/input/nucl
        busco -f --in results/{params.prefix}/busco/input/prot --mode protein --lineage_dataset saccharomycetes_odb10 --out_path results/{params.prefix}/busco/ -c {threads} --out busco_output_prot --offline --download_path {params.busco_db}
        busco -f --in results/{params.prefix}/busco/input/nucl --mode genome --lineage_dataset saccharomycetes_odb10 --out_path results/{params.prefix}/busco/ -c {threads} --out busco_output_nucl --offline --download_path {params.busco_db}
        rm -rf RM_*
        """



# this needs to be updated at the end (after busco_final)
rule multiqc:
    input:
        quast_out = expand("results/{prefix}/quast/{sample}/{sample}_flye_medaka_polypolish/report.txt",prefix=PREFIX, sample=SAMPLE),
        nanoplot_preqc = expand("results/{prefix}/nanoplot/{sample}/{sample}_preqcNanoPlot-report.html",prefix=PREFIX, sample=SAMPLE),
        busco_out_p = "results/{prefix}/busco/busco_output_prot/batch_summary.txt",
        busco_out_n = "results/{prefix}/busco/busco_output_nucl/batch_summary.txt",
        #unicycler_annotation = "results/{prefix}/prokka/",
    output:
        multiqc_report = "results/{prefix}/multiqc/{prefix}_QC_report.html",
    params:
        outdir = "results/{prefix}/multiqc",
        prefix = "{prefix}",
        quast_dir = "results/{prefix}/quast/",
        nanoplot_dir = "results/{prefix}/nanoplot/",
        busco_n_dir = "results/{prefix}/busco/busco_output_nucl/",
        busco_p_dir = "results/{prefix}/busco/busco_output_prot/",
        #multiqc_out_dir = directory("results/{prefix}/multiqc/"),
        #prokka_dir_out = directory("results/{prefix}/prokka"),
        #quast_dir_out = directory("results/{prefix}/quast"),
        #pycoqc_dir_out = directory("results/{prefix}/pycoqc"),
        #busco_dir_out = directory("results/{prefix}/busco"),
        #threads = config["ncores"],
    #priority: 100000    
    # conda:
    #     "envs/multiqc.yaml"
    resources:
        mem_mb = 1000,
        runtime = 20
    threads: 1
    singularity:
        "docker://multiqc/multiqc:v1.25.1"
    shell:
        """
        multiqc -f --outdir {params.outdir} -n {params.prefix}_QC_report -i {params.prefix}_QC_report \
        {params.quast_dir} {params.busco_n_dir} {params.busco_p_dir} {params.nanoplot_dir}
        """

