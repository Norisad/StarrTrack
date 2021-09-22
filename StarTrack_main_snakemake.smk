# COMMAND TO LAUNCH IT : `snakemake -s snake_capstarr_rep_by_rep.smk -c 'qsub -V -q tagc -l nodes=1:ppn=15' -j 6`
# If threads is specified in config.yaml `snakemake -s snake_capstarr_rep_by_rep.smk -c 'qsub -V -q tagc -l nodes=1:ppn={threads}' -j 6`
# in local mode with job in parallel  `snakemake -s snake_capstarr_rep_by_rep.smk -c 'qsub -V -q tagc -l nodes=1:ppn=15' -j 6`
# Command to unlock dir if needed : `snakemake --unlock -s snake_file.py`
# DAG of jobs : snakemake -s snake_capstarr_rep_by_rep.smk --dag | dot -Tsvg > dag.svg

#snakemake -s snake_capstarr_rep_by_rep.smk -c 'qsub -V -q tagc -l nodes=1:ppn=15' -j 6 --use-conda

from function_snakemake import *
import sys
import os
import re
import pandas as pd
import numpy as np


#configfile: "config.yaml"
#instead of set env -> setup config file to use appropriate conda env

#Working directory
workdir:'/gpfs/tagc/home/sadouni/silencer_project/core_silencer_analysis/scp1'
# workdir:'/gpfs/tagc/home/sadouni/silencer_project/core_silencer_analysis/scp1' # on mac

# Dataset : For now here then I have to setup config.yaml to don't touch this file
CDNA,=glob_wildcards('/gpfs/tagc/home/sadouni/silencer_project/core_silencer_analysis/scp1/data/input/cdna/{cdna}.bam')
LIBRARY,=glob_wildcards('/gpfs/tagc/home/sadouni/silencer_project/core_silencer_analysis/scp1/data/input/library/{lib}.bam')
COND,=glob_wildcards('/gpfs/tagc/home/sadouni/silencer_project/core_silencer_analysis/scp1/data/input/{cond}.bam')

# localrule all:
#     input: data/annotated_clones/rep2_paste_input.annotated.bed, data/annotated_clones/rep1_paste_input.annotated.bed, data/wc_bed/cdna/rep2.wc.txt, data/wc_bed/cdna/rep1.wc.txt, data/wc_bed/library/input.wc.txt



rule all :
    input:
        expand("data/bed_rgb/by_clones/{cdna}_paste_{lib}.rgb.bed",cdna=CDNA,lib=LIBRARY),
        expand("data/bed_rgb/by_clones/merged_replicates_paste_{lib}.rgb.bed",lib=LIBRARY),
        expand("data/bed_rgb/by_region_subregion/{cdna}_paste_{lib}.considered_region.rgb.bed",cdna=CDNA,lib=LIBRARY),
        expand("data/bed_rgb/by_region_subregion/merged_replicates_paste_{lib}.considered_region.rgb.bed",lib=LIBRARY),
        expand("data/summary/{cdna}_paste_{lib}.summary.txt",cdna=CDNA,lib=LIBRARY),
        expand("data/summary/merged_replicates_paste_{lib}.summary.txt",lib=LIBRARY),
        expand("data/R_envs_save/{cdna}_paste_{lib}.env.RData",cdna=CDNA,lib=LIBRARY),
        expand("data/R_envs_save/merged_replicates_paste_{lib}.env.RData",lib=LIBRARY),
        expand("data/core_silencer/{cdna}_paste_{lib}.region.silencer.bed",cdna=CDNA,lib=LIBRARY),
        expand("data/core_silencer/{cdna}_paste_{lib}.core_silencer.bed",cdna=CDNA,lib=LIBRARY),
        expand("data/core_silencer/merged_replicates_paste_{lib}.region.silencer.bed",lib=LIBRARY),
        expand("data/core_silencer/merged_replicates_paste_{lib}.core_silencer.bed",lib=LIBRARY)


# rule get_normalization_value:

# rule make_frag_great_again:
#     input:
#         clone_list="clone_list/clone_list.1bp.{cond}.count.bed"
#     params:
#         elong=estimate_frag_size("/home/sadouni/silencer_project/silencer_project/capstarr-pipeline/t5_peaks.xls"),genome="/home/sadouni/genome/mm9/mm9_chr_size.txt"
#     output:
#         "clone_list/clone_list.1bp.{cond}.slop.count.bed"
#     shell:'''
#         bedtools slop -b {params.elong} -i {input.clone_list} -g {params.genome}
#     '''

# rule get_estimate_frag_size:
#     input:
#         "temp/macs_est.xls"
#     output:
#         "temp/estimate_size_frag.txt"
#     run:
#         estimate_frag_size("/home/sadouni/silencer_project/silencer_project/capstarr-pipeline/t5_peaks.xls","temp/estimate_size_frag.txt")

# rule macs_estimate_frag_size:
#     input:
#         expand("bed_files/{cond}.bed",cond=CONDITION)
#     output:
#         "temp/macs_est.xls"
#     shell: '''
#         cat {input} | sort -k1,1 -k2,2n -k6,6 | awk -v OFS='\t' '/^chr/{{$4=".";print $0}}' > temp/concat_bed.bed
#         macs2 callpeak -t temp/concat_bed.bed -f BED -g mm -n macs_est
#     '''

# rule concatenate_replicate:
#     input:
#         expand("data/annotated_clones/{cdna}_paste_{lib}.annotated.bed", cdna=CDNA,lib=LIBRARY)
#     output:
#         "data/merge_replicates/merge.annotated.bed"
#     shell:'''
#         cat {input} > {output}
#     '''

#Compute RGB color code proportional to FC
rule compute_rgb:
    input:
        "data/tmp/fold_change/{cdna}_paste_{lib}.fc.bed"
    output:
        "data/bed_rgb/by_clones/{cdna}_paste_{lib}.rgb.bed"
    params:
        bed_header = "{cdna}"
    shell: '''
        awk -v maxlogratio=5 -v FS='\t' -v OFS='\t' '/^chr/{{r=0;g=0;if($9<1){{r=-255*(log($9)/log(2))/maxlogratio;if(r>255){{r=255}}}};if($9>1){{g=255*(log($9)/log(2))/maxlogratio;if(g>255){{g=255}}}};print $0,r","g",0"}}' {input} | awk -v FS='\t' -v OFS='\t' 'BEGIN{{print"track","type=bed","name={params.bed_header}","itemRgb=on"}} /^chr/{{print $1,$2,$3,".",$5,$6,$2,$3,$10}}' > {output}
    '''

#compute FC -> normalisation using number of read
rule compute_fc:
    input:
        annotate= "data/tmp/annotated_clones/{cdna}_paste_{lib}.annotated.bed",
        wc_cdna= "data/tmp/wc_bed/cdna/{cdna}.wc.txt",
        wc_lib= "data/tmp/wc_bed/library/{lib}.wc.txt"
    output:
        "data/tmp/fold_change/{cdna}_paste_{lib}.fc.bed"
    shell: '''
        SIZE_CDNA=`cut -f1 -d ' ' {input.wc_cdna}`
        SIZE_LIB=`cut -f1 -d ' ' {input.wc_lib}`
        awk -v cdna1=$SIZE_CDNA -v inp=$SIZE_LIB -v pseudocount=1 -v FS='\t' -v OFS='\t' '/^chr/{{ratio=(($7+pseudocount)/($8+pseudocount))*(inp/cdna1);print $0,ratio}}' {input.annotate} | cut -f 1-8,11 > {output}
    '''

# R analysis
rule find_core_edge_silencer:
    input:
        rep="data/tmp/annotated_clones_subregion/{cdna}_paste_{lib}.enlarged_dhs.considered_region.subregion.bed"
    output:
        summary = "data/summary/{cdna}_paste_{lib}.summary.txt",
        region_silencer_out = "data/silencer_files/{cdna}_paste_{lib}.region.silencer.bed",
        region_enhancer_out = "data/enhancer_files/{cdna}_paste_{lib}.region.silencer.bed",
        region_rgb = "data/bed_rgb/by_region_subregion/{cdna}_paste_{lib}.considered_region.rgb.bed",

        subregion_silencer_out = "data/silencer_files/{cdna}_paste_{lib}.core_silencer.bed",
        subregion_enhancer_out = "data/enhancer_files/{cdna}_paste_{lib}.core_enhancer.bed",
        out_env = "data/R_envs_save/{cdna}_paste_{lib}.env.RData"
    conda:
        "envs/r_dev.yaml"
    shell: '''
        Rscript --vanilla subregion_analysis.R
        {input.rep}
        {output.summary}
        {output.region_silencer_out}
        {output.region_enhancer_out}
        {output.region_rgb}
        {output.subregion_silencer_out}
        {output.subregion_enhancer_out}
        {output.out_env}
        '''




rule annotated_considered_region_subregion_bins:
    input:
        clones="data/tmp/slop_clone_list/{cdna}_paste_{lib}.bed", region_annotated="data/tmp/enlarged_coordinates/{cdna}/{cdna}_paste_{lib}.enlarged_dhs.considered_region.subregion.bed"
    output:
        "data/tmp/annotated_clones_subregion/{cdna}_paste_{lib}.enlarged_dhs.considered_region.subregion.bed"
    shell: '''
        set +u; source /gpfs/tagc/home/sadouni/Apps/anaconda3/bin/activate dev; set -u
        awk -v FS='\t' -v OFS='\t' '/^chr/{{$4=".";print $0}}' {input.clones} | intersectBed -wao -a stdin -b {input.region_annotated} | cut -f 1-8,10-14 > {output}
    '''

rule make_bins_considered_enlarged_region:
    input:
        bed_in="data/tmp/enlarged_coordinates/{cdna}/{cdna}_paste_{lib}.enlarged_dhs.considered_region.bed"
    output:
        bed_out="data/tmp/enlarged_coordinates/{cdna}/{cdna}_paste_{lib}.enlarged_dhs.considered_region.subregion.bed"
    params:
        bins=50
    run:
        make_bins(params.bins,input.bed_in,output.bed_out)

rule get_considered_enlarged_region:
    input:
        "data/tmp/annotated_clones_enlarged_dhs/{cdna}/{cdna}_paste_{lib}.annotated.enlarged_dhs.bed"
    output:
        "data/tmp/enlarged_coordinates/{cdna}/{cdna}_paste_{lib}.enlarged_dhs.considered_region.bed"
    run:
        enlarge_overlapping_region(input[0],output[0])

#Annotate the clones with extended DHS
rule annotated_captured_regions_enlarged_dhs:
    input:
        clones="data/tmp/slop_clone_list/{cdna}_paste_{lib}.bed",dhs="data/tmp/enlarged_coordinates/{cdna}/{cdna}_paste_{lib}.enlarged_dhs.bed"
    output:
        "data/tmp/annotated_clones_enlarged_dhs/{cdna}/{cdna}_paste_{lib}.annotated.enlarged_dhs.bed"
    shell: '''
        set +u; source /gpfs/tagc/home/sadouni/Apps/anaconda3/bin/activate dev; set -u
        awk -v FS='\t' -v OFS='\t' '/^chr/{{$4=".";print $0}}' {input.clones} | intersectBed -wao -a stdin -b {input.dhs} | cut -f 1-8,12-13 > {output}
    '''

#get coordinates of extend DHS
rule get_enlarged_dhs:
    input:
        "data/tmp/annotated_clones/{cdna}_paste_{lib}.annotated.bed"
    output:
        "data/tmp/enlarged_coordinates/{cdna}/{cdna}_paste_{lib}.enlarged_dhs.bed"
    run:
        # lambda wildcards: enlarge_overlapping_region(f"{wildcards.input}",f"{wildcards.output}")
        enlarge_overlapping_region(input[0],output[0])

#Annotate each clone by overlap region of interest
rule annotated_captured_regions:
    input:
        "data/tmp/slop_clone_list/{files}.bed"
    output:
        "data/tmp/annotated_clones/{files}.annotated.bed"
    params:
        region_type="/gpfs/tagc/home/sadouni/genome/mm9/DHS_capture_regions_mm9.bed"
    shell: '''
        set +u; source /gpfs/tagc/home/sadouni/Apps/anaconda3/bin/activate dev; set -u
        awk -v FS='\t' -v OFS='\t' '/^chr/{{$4=".";print $0}}' {input} | intersectBed -wao -a stdin -b {params.region_type} | cut -f 1-8,12-13 > {output}
    '''

#slop to original frag size for now fix length by hand
rule slop_bed:
    input:
        "data/tmp/paste/{files}.bed"
    output:
        "data/tmp/slop_clone_list/{files}.bed"
    params:
        frag_size=314 ,genome="/gpfs/tagc/home/sadouni/genome/mm9/mm9_chr_size.txt"
    shell:'''
        set +u; source /gpfs/tagc/home/sadouni/Apps/anaconda3/bin/activate dev; set -u
        bedtools slop -s -l 0 -r {params.frag_size} -i {input} -g {params.genome}  > {output}
    '''

#Paste cDNA with respective Library input
rule paste_cdna_input:
    input:
        dna="data/tmp/count/{cdna}/{cdna}_vs_{lib}.{cdna}_clone_{lib}.count.bed",
        lib="data/tmp/count/{cdna}/{lib}_vs_{cdna}.{cdna}_clone_{lib}.count.bed"
    output:
       "data/tmp/paste/{cdna}_paste_{lib}.bed"
    shell:'''
        paste {input.dna} {input.lib} | cut -f 1-7,14  > {output}
    '''

#Count the frequency of clone. Have to be done at each time for cdna and library
rule bedtools_intersect_to_count:
    input:
        cdna="data/bed_files/cdna/{cdna}.bed",lib="data/bed_files/library/{lib}.bed",clone_list="data/tmp/clone_list_1bp/{cdna}/{clone_list}.bed"
    output:
        cdna_out="data/tmp/count/{cdna}/{cdna}_vs_{lib}.{clone_list}.count.bed",lib_out="data/tmp/count/{cdna}/{lib}_vs_{cdna}.{clone_list}.count.bed"
    shell:'''
        set +u; source /gpfs/tagc/home/sadouni/Apps/anaconda3/bin/activate dev; set -u
        awk -v FS='\t' -v OFS='\t' '((/^chr/)&&($2>=0)){{$3=$2+1;print $0}}' {input.cdna} | intersectBed -c -a {input.clone_list} -b stdin > {output.cdna_out}
        awk -v FS='\t' -v OFS='\t' '((/^chr/)&&($2>=0)){{$3=$2+1;print $0}}' {input.lib} | intersectBed -c -a {input.clone_list} -b stdin > {output.lib_out}
    '''

#create a list of clone where all cDNA and library conditions are concatenated
rule concat_clone_list:
    input:
        cdna="data/bed_files/cdna/{cdna}.bed",lib="data/bed_files/library/{lib}.bed"
    output:
        "data/tmp/clone_list_1bp/{cdna}/{cdna}_clone_{lib}.bed"
    shell:'''
        cat {input.cdna} {input.lib} | sort -u -k1,1 -k2,2n -k6,6 | awk -v OFS='\t' '/^chr/{{$4=".";print $0}}' | awk -v FS='\t' -v OFS='\t' '((/^chr/)&&($2>=0)){{$3=$2+1;print $0}}'  > {output}
    '''

rule wc_bed_files:
    input:
        "data/bed_files/{files}.bed",
    output:
        "data/tmp/wc_bed/{files}.wc.txt"
    shell:'''
        wc -l {input} > {output}
    '''
#convert bam file to bed
rule bam_to_bed:
    input:
        rep="data/input/{files}.bam"
    output:
        rep_out="data/bed_files/{files}.bed"
    shell:'''
        set +u; source /gpfs/tagc/home/sadouni/Apps/anaconda3/bin/activate dev; set -u
        bedtools bamtobed -i {input.rep} | sort -k1,1 -k2,2n -k3,3n  > {output.rep_out}
    '''

rule merge_bam:
    input:
        expand("data/input/cdna/{cdna}.bam",cdna=CDNA)
    output:
        "data/input/cdna/merged_replicates.bam"
    shell:'''
    set +u; source /gpfs/tagc/home/sadouni/Apps/anaconda3/bin/activate dev; set -u
    samtools merge {output} {input}
    '''
