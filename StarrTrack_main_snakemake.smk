##### Running commands ######
# COMMAND TO LAUNCH IT :
# snakemake -s StarrTrack/StarrTrack_main_snakemake.smk -c 'qsub -V -q tagc -l nodes=1:ppn={threads}' -j 6 --use-conda
# Command to unlock dir if needed : `snakemake --unlock -s snake_file.py`
# DAG of jobs : snakemake -s StarrTrack/StarrTrack_main_snakemake.smk --dag | dot -Tsvg > dag.svg

#Load attached script
configfile: "StarrTrack/config.yaml"
from function_snakemake import *


#Working directory
# workdir:'/gpfs/tagc/home/sadouni/integrative_reporter_assay/run423'
workdir: config['outdir']

#path to input files
#cdna and library have to be separate by folder cdna / library in input folder
# example : data/input/cdna
#           data/input/library
# path_to_files = config['path_to_files']

# Dataset
CDNA,=glob_wildcards(config['path_to_files']+'/cdna/{cdna}.bam')
LIBRARY,=glob_wildcards(config['path_to_files']+'/library/{lib}.bam')
COND,=glob_wildcards(config['path_to_files']+'/{cond}.bam')




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
        expand("data/bed_activity/{cdna}_paste_{lib}.region.silencer.bed",cdna=CDNA,lib=LIBRARY),
        expand("data/bed_activity/{cdna}_paste_{lib}.core_silencer.bed",cdna=CDNA,lib=LIBRARY),
        expand("data/bed_activity/{cdna}_paste_{lib}.region.enhancer.bed",cdna=CDNA,lib=LIBRARY),
        expand("data/bed_activity/{cdna}_paste_{lib}.core_enhancer.bed",cdna=CDNA,lib=LIBRARY),
        expand("data/bed_activity/merged_replicates_paste_{lib}.region.silencer.bed",lib=LIBRARY),
        expand("data/bed_activity/merged_replicates_paste_{lib}.core_silencer.bed",lib=LIBRARY),
        expand("data/bed_activity/merged_replicates_paste_{lib}.region.enhancer.bed",lib=LIBRARY),
        expand("data/bed_activity/merged_replicates_paste_{lib}.core_enhancer.bed",lib=LIBRARY)


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
        region_silencer_out = "data/bed_activity/{cdna}_paste_{lib}.region.silencer.bed",
        region_enhancer_out = "data/bed_activity/{cdna}_paste_{lib}.region.enhancer.bed",
        region_rgb = "data/bed_rgb/by_region_subregion/{cdna}_paste_{lib}.considered_region.rgb.bed",
        subregion_silencer_out = "data/bed_activity/{cdna}_paste_{lib}.core_silencer.bed",
        subregion_enhancer_out = "data/bed_activity/{cdna}_paste_{lib}.core_enhancer.bed",
        out_env = "data/R_envs_save/{cdna}_paste_{lib}.env.RData"
    conda:
        "envs/starrtrak_envs.yaml"
    shell: '''
        Rscript --vanilla StarrTrack/subregion_analysis.R {input.rep} {output.summary} {output.region_rgb} {output.region_silencer_out} {output.region_enhancer_out} {output.subregion_silencer_out} {output.subregion_enhancer_out} {output.out_env}
        '''

rule annotated_considered_region_subregion_bins:
    input:
        clones="data/tmp/slop_clone_list/{cdna}_paste_{lib}.bed", region_annotated="data/tmp/enlarged_coordinates/{cdna}/{cdna}_paste_{lib}.enlarged_dhs.considered_region.subregion.bed"
    output:
        "data/tmp/annotated_clones_subregion/{cdna}_paste_{lib}.enlarged_dhs.considered_region.subregion.bed"
    conda:
        "envs/starrtrak_envs.yaml"
    shell: '''
        awk -v FS='\t' -v OFS='\t' '/^chr/{{$4=".";print $0}}' {input.clones} | intersectBed -wao -a stdin -b {input.region_annotated} | cut -f 1-8,10-14 > {output}
    '''

rule make_bins_considered_enlarged_region:
    input:
        bed_in="data/tmp/enlarged_coordinates/{cdna}/{cdna}_paste_{lib}.enlarged_dhs.considered_region.bed"
    output:
        bed_out="data/tmp/enlarged_coordinates/{cdna}/{cdna}_paste_{lib}.enlarged_dhs.considered_region.subregion.bed"
    params:
        bins=50
    # conda:
    #     "StarrTrack/envs/starrtrak_envs.yaml"
    run:
        make_bins(params.bins,input.bed_in,output.bed_out)

rule get_considered_enlarged_region:
    input:
        "data/tmp/annotated_clones_enlarged_dhs/{cdna}/{cdna}_paste_{lib}.annotated.enlarged_dhs.bed"
    output:
        "data/tmp/enlarged_coordinates/{cdna}/{cdna}_paste_{lib}.enlarged_dhs.considered_region.bed"
    # conda:
    #     "StarrTrack/envs/starrtrak_envs.yaml"
    run:
        enlarge_overlapping_region(input[0],output[0])

#Annotate the clones with extended DHS
rule annotated_captured_regions_enlarged_dhs:
    input:
        clones="data/tmp/slop_clone_list/{cdna}_paste_{lib}.bed",dhs="data/tmp/enlarged_coordinates/{cdna}/{cdna}_paste_{lib}.enlarged_dhs.bed"
    output:
        "data/tmp/annotated_clones_enlarged_dhs/{cdna}/{cdna}_paste_{lib}.annotated.enlarged_dhs.bed"
    conda:
        "envs/starrtrak_envs.yaml"
    shell: '''
        awk -v FS='\t' -v OFS='\t' '/^chr/{{$4=".";print $0}}' {input.clones} | intersectBed -wao -a stdin -b {input.dhs} | cut -f 1-8,12-13 > {output}
    '''

#get coordinates of extend DHS
rule get_enlarged_dhs:
    input:
        "data/tmp/annotated_clones/{cdna}_paste_{lib}.annotated.bed"
    output:
        "data/tmp/enlarged_coordinates/{cdna}/{cdna}_paste_{lib}.enlarged_dhs.bed"
    # conda:
    #     "StarrTrack/envs/starrtrak_envs.yaml"
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
        # region_type="/gpfs/tagc/home/sadouni/genome/captured_regions/hProm.hg19.annotated_id.bed"
        region_type=config['captured_regions_ref']
    conda:
        "envs/starrtrak_envs.yaml"
    shell: '''
        awk -v FS='\t' -v OFS='\t' '/^chr/{{$4=".";print $0}}' {input} | intersectBed -wao -a stdin -b {params.region_type} | cut -f 1-8,12-13 > {output}
    '''

#slop to original frag size for now fix length by hand
rule slop_bed:
    input:
        "data/tmp/paste/{files}.bed"
    output:
        "data/tmp/slop_clone_list/{files}.bed"
    params:
        # frag_size=314 ,genome="/gpfs/tagc/home/sadouni/genome/hg19/hg19.chr_size.txt"
        frag_size=314 ,genome=config['chr_size_file']

    conda:
        "envs/starrtrak_envs.yaml"
    shell:'''
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
        cdna="data/tmp/ori_bed_files/cdna/{cdna}.original_frag_length.bed",lib="data/tmp/ori_bed_files/library/{lib}.original_frag_length.bed",clone_list="data/tmp/clone_list_1bp/{cdna}/{clone_list}.bed"
    output:
        cdna_out="data/tmp/count/{cdna}/{cdna}_vs_{lib}.{clone_list}.count.bed",lib_out="data/tmp/count/{cdna}/{lib}_vs_{cdna}.{clone_list}.count.bed"
    conda:
        "envs/starrtrak_envs.yaml"
    shell:'''
        awk -v FS='\t' -v OFS='\t' '((/^chr/)&&($2>=0)){{$3=$2+1;print $0}}' {input.cdna} | intersectBed -c -a {input.clone_list} -b stdin > {output.cdna_out}
        awk -v FS='\t' -v OFS='\t' '((/^chr/)&&($2>=0)){{$3=$2+1;print $0}}' {input.lib} | intersectBed -c -a {input.clone_list} -b stdin > {output.lib_out}
    '''

#create a list of clone where all cDNA and library conditions are concatenated
rule concat_clone_list:
    input:
        cdna="data/tmp/ori_bed_files/cdna/{cdna}.original_frag_length.bed",lib="data/tmp/ori_bed_files/library/{lib}.original_frag_length.bed"
    output:
        "data/tmp/clone_list_1bp/{cdna}/{cdna}_clone_{lib}.bed"
    shell:'''
        cat {input.cdna} {input.lib} | sort -u -k1,1 -k2,2n -k6,6 | awk -v OFS='\t' '/^chr/{{$4=".";print $0}}' | awk -v FS='\t' -v OFS='\t' '((/^chr/)&&($2>=0)){{$3=$2+1;print $0}}' > {output}
    '''

#restore original frag lenth
rule original_frag_length:
    input:
        clone_list_ind="data/bed_files/{files}.bed"
    output:
        original_frag_length="data/tmp/ori_bed_files/{files}.original_frag_length.bed"
    params:
        frag_size=314
    shell:'''
        awk -v OFS='\t' '/^chr/{{if($6=="+"){{$3=$2+{params.frag_size}-1}}else{{$2=$3-{params.frag_size}+1}}print$0}}' {input.clone_list_ind} > {output.original_frag_length}
    '''

#Get normalisation number
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
    conda:
        "envs/starrtrak_envs.yaml"
    shell:'''
        bedtools bamtobed -i {input.rep} | sort -k1,1 -k2,2n -k3,3n  > {output.rep_out}
    '''
#merge all merged_replicates
#analysis is done in parallele for separate and merged replicates
rule merge_bam:
    input:
        expand("data/input/cdna/{cdna}.bam",cdna=CDNA)
    output:
        "data/input/cdna/merged_replicates.bam"
    conda:
        "envs/starrtrak_envs.yaml"
    shell:'''
    samtools merge {output} {input}
    '''
