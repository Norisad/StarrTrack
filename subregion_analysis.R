#!/usr/bin/env Rscript
##############################################
##### CapSTARR-seq Statistical Analysis ######
############ Subregions ######################
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "out.txt"
}
#Library
library("ggplot2")
library("dplyr")
library("ggpubr")
library("data.table")

#Functions
#Compute activity using a list of clones and dataframe
# get_activities=function(list_of_clones,df){
#   summary_activities= data.frame(
#     region_totals_cdna = sapply( list_of_clones, function(x) sum( df[x, "cdna_count"] )),
#     region_totals_lib = sapply( list_of_clones, function(x) sum( df[x, "lib_count"] )))
#   # region_activity_repX is calculated using replicates separatly
#   summary_activities$region_activity = apply(summary_activities,1,function(x){(x[1]/x[2])* ((sum_count_full_df["coverage_lib"]+sum_count_bg["lib_count"])/ (sum_count_full_df["coverage_cdna"]+sum_count_bg["cdna_count"]))})
#   return(summary_activities)}

# get_activities=function(list_of_clones,df){
#   summary_activities= data.frame(
#     region_totals_cdna = sapply( list_of_clones, function(x) sum( df[x, "cdna_count"] )),
#     region_totals_lib = sapply( list_of_clones, function(x) sum( df[x, "lib_count"] )))
#   # region_activity_repX is calculated using replicates separatly
#   summary_activities$region_activity = apply(summary_activities,1,function(x){(x[1]/x[2])* ((sum_count["lib_count"]+sum_count_bg["lib_count"])/ (sum_count["cdna_count"]+sum_count_bg["cdna_count"]))}) 
#   return(summary_activities)}

#Compute RPKM
# Normalize each count region by the sequencing depth /1.10^6 bp
# Then normalize by the size of the region 
get_activities_normalized_coverage=function(list_of_clones,df){
  summary_activities= data.frame(
    fpkm_cdna = sapply( list_of_clones, function(x) 
      ((sum( df[x,"coverage_cdna_region"]))/
         ((sum_count_full_df["coverage_cdna"]+sum_count_bg["cdna_count"])/1000000))/(
           mean(df[x,"end_region"]-df[x,"start_region"])/1000)),
    fpkm_lib = sapply( list_of_clones, function(x) 
      ((sum( df[x,"coverage_lib_region"]))/
         ((sum_count_full_df["coverage_lib"]+sum_count_bg["lib_count"])/1000000))/(
           mean(df[x,"end_region"]-df[x,"start_region"])/1000))
  )
  summary_activities$region_activity = apply(summary_activities,1,function(x){log2((x[1]/x[2]))})
  return(summary_activities)}

#Function to have RGB code from red to green / palette
scale_to_rgb <- function(val, bounds = NA,
                         colors = c("#ff0000", "#000000", "#00ff00"),
                         format = c("rgb", "comma"), ...) {
  if (anyNA(bounds)) bounds <- range(val, na.rm = TRUE)
  format = match.arg(format)
  
  isna <- is.na(val)
  ispos <- !isna & val >= 0
  isneg <- !isna & !ispos
  cols <- matrix(NA, nrow = length(val), ncol = 3)
  valneg <- pmax(bounds[1], val[isneg]) / bounds[1]
  valpos <- pmin(bounds[2], val[ispos]) / bounds[2]
  
  cols[isneg,] <- colorRamp(colors[2:1], ...)(valneg)
  cols[ispos,] <- colorRamp(colors[2:3], ...)(valpos)
  
  if (format == "rgb") {
    cols <- cols / 255
    rgb(cols[,1], cols[,2], cols[,3])
  } else {
    cols <- round(cols, 0)
    paste(cols[,1], cols[,2], cols[,3], sep = ",")
  }
}

#Compute the inflexion point
calculate_cutoff <- function(inputVector, drawPlot=TRUE,...){
  inputVector <- sort(inputVector)
  inputVector[inputVector<0] <- 0 #set those regions with more control than ranking equal to zero
  slope <- (max(inputVector)-min(inputVector))/length(inputVector) #This is the slope of the line we want to slide. This is the diagonal.
  xPt <- floor(
    optimize(
      numPts_below_line,
      lower = 1,
      upper = length(inputVector),
      myVector = inputVector,
      slope = slope
    )$minimum
  ) #Find the x-axis point where a line passing through that point has the minimum number of points below it. (ie. tangent)
  y_cutoff <- inputVector[xPt] #The y-value at this x point. This is our cutoff.
  
  if(drawPlot){  #if TRUE, draw the plot
    plot(
      1:length(inputVector),
      inputVector,
      type = "l",
      ...
    )
    b <- y_cutoff-(slope* xPt)
    abline(
      v = xPt,
      h = y_cutoff,
      lty = 2,
      col = 8
    )
    points(
      xPt,
      y_cutoff,
      pch = 16,
      cex = 0.9,
      col = 2
    )
    abline(
      coef = c(b,slope),
      col = 2
    )
    title(
      paste(
        "x=",
        xPt,
        "\ny=",
        signif(
          y_cutoff,
          3
        ),
        "\nFold over Median=",
        signif(
          y_cutoff/median(inputVector),
          3
        ),
        "x\nFold over Mean=",
        signif(
          y_cutoff/mean(inputVector),
          3
        ),
        "x",
        sep=""
      )
    )
    axis(
      1,
      sum(inputVector==0),
      sum(inputVector==0),
      col.axis = "pink",
      col = "pink"
    ) #Number of regions with zero signal
  }
  return(
    list(
      absolute = y_cutoff,
      overMedian = y_cutoff/median(inputVector),
      overMean = y_cutoff/mean(inputVector)
    )
  )
}

numPts_below_line <- function(myVector,slope,x){
  yPt <- myVector[x]
  b <- yPt-(slope*x)
  xPts <- 1:length(myVector)
  return(sum(myVector<=(xPts*slope+b)))
}

# get_coverage=function(list_of_clones,df){
#   summary_coverage= data.frame(
#     coverage_cdna = sapply( list_of_clones, function(x) sum( df[!duplicated(df[,c('chr','start_clone','end_clone')]),][x, "cdna_count"] )),
#     coverage_lib = sapply( list_of_clones, function(x) sum( df[!duplicated(df[,c('chr','start_clone','end_clone')]),][x, "lib_count"] )))
#   # region_activity_repX is calculated using replicates separatly
#   return(summary_coverage)}


############################################################################
########### Threshold used to determinate enhancer / silencer ##############
############################################################################

enhancer_threshold=1 # need to be change to > than inflexion point
silencer_threshold=-1

###################
#### Load Data ####
###################

# input_clone_dataf=as.data.frame(read.table("/Users/nori/sacapus_remote/silencer_project/silencer_pipeline_output_first_part/pgk/data/annotated_clones_subregion/chr11.pPGK_r1.bed",sep="\t", dec=".",fill=TRUE,header = FALSE,col.names=(c("chr","start_clone","end_clone","point","score","strands","cdna_count","lib_count","start_subregion","end_subregion","region_type","region_id","subregion_id"))))

# input_clone_dataf=as.data.frame(read.table("sacapus_remote/silencer_project/silencer_pipeline_output_first_part/pgk/data/annotated_clones_subregion/chr11.pPGK.r3.bed",sep="\t", dec=".",fill=TRUE,header = FALSE,col.names=(c("chr","start_clone","end_clone","point","score","strands","cdna_count","lib_count","start_subregion","end_subregion","region_type","region_id","subregion_id"))))

#Input in command line
input_clone_dataf = as.data.frame(read.table(args[1],sep="\t", dec=".",fill=TRUE,header = FALSE,col.names=(c("chr","start_clone","end_clone","point","score","strands","cdna_count","lib_count","start_subregion","end_subregion","region_type","region_id","subregion_id"))))
#Use the file name to extract the condition name -> use on bed RGB file
header_rgb=paste(c("track","type=bed",paste("name",unlist(strsplit(unlist(strsplit(args[1],"\\/"))[3],"\\."))[1],sep = "="),"itemRgb=on\n"),collapse  = "\t")

#Annotation file
# Need to change path to arg[x] tu be used on command line
DHS_annotated=as.data.frame(read.table("/gpfs/tagc/home/sadouni/genome/mm9/DHS_region_id_annotated_gene.txt",sep="\t", dec=".",fill=TRUE,header = FALSE,col.names = c("chr","start_dhs","end_dhs","region_type","region_id","genes")))

#Add one to each column to avoid 0 division
input_clone_dataf[, 7:8]=input_clone_dataf[, 7:8] + 1

#Remove non annotated clone : decrease time computation
background_df<-input_clone_dataf[(input_clone_dataf$region_id=="."),]
input_clone_dataf<-input_clone_dataf[!(input_clone_dataf$region_id=="."),]

##########################
#### Regions activity ####
##########################
#Overall count : Make the sum of different column : Normalization values
# sum_count=apply(input_clone_dataf[,c("cdna_count","lib_count")],2,sum)
sum_count_bg=apply(background_df[,c("cdna_count","lib_count")],2,sum)
sum_count=apply(input_clone_dataf[,c("cdna_count","lib_count")],2,sum)

#Make a List of region ID / subregion ID
region_list <- levels(input_clone_dataf$region_id)
subregion_list <- levels(input_clone_dataf$subregion_id)

#List of all unique region ID / subregion ID which contain all the clone that belong this region/subregion /!\ Long step
# region_list_clones <- lapply(region_list, function(x) which(input_clone_dataf$region_id == x))
# subregion_list_clones <- lapply(subregion_list, function(x) which(input_clone_dataf$subregion_id == x))


#Get coverage by region
# To compute subregion activity some clones are present several times per region
# Need to count each clone one times per region
#If 1 clone overlaps 2 subregions then need to divid by 2 the counts
# Will be used to get the coverage not normalized of the region
# This step is not requirend to compute coverage of subregion because
# you cannot have 2 sames clones overlapping the same subregion
sum_input_clone_dataf=input_clone_dataf %>%
  group_by(chr,start_clone,end_clone,region_id,strands) %>%
  summarise_at(c("cdna_count","lib_count"),mean)
# This step is not necessary for subregion because we don't have duplicates clones per subregions

#Make the sums of the clones that belong a region
cov_input_region_clone_dataf=sum_input_clone_dataf %>%
  group_by(region_id) %>%
  summarise_at(c("cdna_count","lib_count"),sum) %>%
  rename("coverage_cdna"="cdna_count","coverage_lib"="lib_count")

#Get coverage by subregion
cov_input_subregion_clone_dataf=input_clone_dataf %>%
  group_by(subregion_id) %>%
  summarise_at(c("cdna_count","lib_count"),sum) %>%
  rename("coverage_cdna"="cdna_count","coverage_lib"="lib_count")

#Get normalization value
# Better to use this one instead of sum_count
# Because some clones are count several times
sum_count_full_df=apply(cov_input_region_clone_dataf[,c("coverage_cdna","coverage_lib")],2,sum)

#Get region / subregion activity
# region_activity=get_activities(region_list_clones,input_clone_dataf)$region_activity
# subregion_activity=get_activities(subregion_list_clones,input_clone_dataf)$region_activity

#Copy of original dataframe -> will be used as output
# output_clone_df=input_clone_dataf

#Subset of original dataframe : merge the original dataframe with the coverage per region
# -> used to have extend region dhs + information of coverage lib / input
region_df=merge(input_clone_dataf[,c("chr","start_subregion","end_subregion","region_id")],cov_input_region_clone_dataf,by = "region_id")

#Get by subregion_id the coverage and activity
subregion_df=merge(input_clone_dataf[,c("chr","start_subregion","end_subregion","subregion_id","region_id","region_type")],cov_input_subregion_clone_dataf,by = "subregion_id")

#remove duplicates due to the presence of clones per subregion 
# I need 1 row per region_id / subregion_id with the coverage
region_df = region_df[!duplicated(region_df[,c('chr','start_subregion','end_subregion','region_id','coverage_cdna','coverage_lib')]),]

subregion_df = subregion_df[!duplicated(subregion_df[,c('chr','start_subregion','end_subregion','subregion_id','coverage_cdna','coverage_lib')]),]

# Compute subregion activity and center it
# No need to use the function get_activities_normalized_coverage 
# I don't need the FPKM per subregion
# So I compute activity using the coverage_cdna/coverage_lib*(sum_lib/sum_cdna)
# check that the activity is correct : done (same as previous activity computed)
subregion_df<-subregion_df%>%
  mutate(subregion_activity=log2((coverage_cdna/coverage_lib)*(sum_count_full_df["coverage_lib"]+sum_count_bg["lib_count"])/(sum_count_full_df["coverage_cdna"]+sum_count_bg["cdna_count"]))) %>%
  mutate(subregion_activity_center=scale(subregion_activity,center = TRUE,scale = FALSE))


#Get the region size
# Start of the region is the start of the first clone overlapping the region
# End is the end of the last clone overlapping the region
final_region_df=setDF(setDT(region_df)[, .(chr=chr[1],start_region=min(start_subregion),end_region=max(end_subregion),coverage_cdna_region=coverage_cdna[1],coverage_lib_region=coverage_lib[1]) ,by=region_id])

#Normalize the coverage using FPKM : Done by get_activities_normalized_coverage
# This function compute FPKM -> log2(activity) = FPKM_cdna/FPKM_library

# First get the list of clones that belong each region
region_list_clones_fpkm <- lapply(region_list, function(x) which(final_region_df$region_id == x))

#Call function to have the fpkm cdna / fpkm input
fpkm_final_region_df= get_activities_normalized_coverage(region_list_clones_fpkm,final_region_df)

#ID to put the value calculated on fpkm_final_region_df in -> final_region_df
final_region_df.region_indexes = match(final_region_df[,"region_id"], region_list)
#Add FPKM value to the dataframe
final_region_df$fpkm_cdna_region=fpkm_final_region_df$fpkm_cdna[final_region_df.region_indexes]
final_region_df$fpkm_lib_region=fpkm_final_region_df$fpkm_lib[final_region_df.region_indexes]
#Activity is already in log2 scale -> function get_activities_normalized_coverage
final_region_df$region_activity=fpkm_final_region_df$region_activity[final_region_df.region_indexes]

#Exclude region which have fpkm < 1
final_region_df_exclude=subset(final_region_df, final_region_df[,"fpkm_lib_region"] >= 1)  
#Center the activity
final_region_df_exclude$region_activity_center=scale(final_region_df_exclude$region_activity,scale = FALSE,center = TRUE)



#Index for region and subregion
# input_clone_dataf.region_indexes = match(output_clone_df[,"region_id"], region_list)
# input_clone_dataf.subregion_indexes = match(output_clone_df[,"subregion_id"], subregion_list)

#Put activity and reads count on df
# output_clone_df$subregion_activity=log2(subregion_activity[input_clone_dataf.subregion_indexes])
# output_clone_df$region_activity=log2(region_activity[input_clone_dataf.region_indexes])
# output_clone_df$coverage_cdna=get_activities(subregion_list_clones,output_clone_df)$region_totals_cdna[input_clone_dataf.subregion_indexes]
# output_clone_df$coverage_lib=get_activities(subregion_list_clones,output_clone_df)$region_totals_lib[input_clone_dataf.subregion_indexes]

#Regroup by subregion ID
# final_subregion_df=setDF(setDT(output_clone_df)[, .(chr=chr[1],
#                                                     region_type=region_type[1],
#                                                     region_id=region_id[1],
#                                                     region_activity=region_activity[1],
#                                                     start_subregion = min(start_subregion), end_subregion = max(end_subregion),
#                                                     # coverage_cdna=coverage_cdna[1],
#                                                     # coverage_lib=coverage_lib[1],
#                                                     subregion_activity=subregion_activity[1]),by=subregion_id])

#Get RPKM  instead of basic coverage
# normalized_cov_subregion=get_activities_normalized_coverage(subregion_list_clones,final_subregion_df)

#Merge the extend region with subregion
final_df_exclude=merge(final_region_df_exclude, subregion_df, by=c("chr","region_id"))%>%
  select(chr,start_subregion,end_subregion,fpkm_cdna_region,fpkm_lib_region,subregion_id,subregion_activity_center,start_region,end_region,region_type,region_id,region_activity_center)# NA's match
#Order by region_id and subregion 
final_df_exclude=final_df_exclude[order(final_df_exclude$region_id,final_df_exclude$start_subregion),]
# this dataframe contains per row the subregion coordinates with subregion activity
# Can be used as bed files to have the coordinates of the region
# Need to add color RGB code by subregion activity
# Maybed need to add the closest genes

#Compute RGB code directly here without function
# final_region_subregion_rgb=final_df_exclude%>%
#   select(chr,start_subregion,end_subregion,subregion_id,subregion_activity_center,start_region,end_region,region_id,region_activity_center)%>%
#   mutate(rgb=ifelse(subregion_activity_center<0, paste(ifelse(-255*(subregion_activity_center)/5 > 255, 255, -255*(subregion_activity_center)/5),0,0,sep = ","),
#                     paste(0, ifelse(255*(subregion_activity_center)/5 > 255, 255,255*(subregion_activity_center)/5),0,sep=",")))

#Use ColorRamp and function : Can be usefull to got also color code palette for this just replace comma by rgb -> look function below
final_region_subregion_rgb=final_df_exclude%>%
  select(chr,start_subregion,end_subregion,subregion_id,subregion_activity_center,start_region,end_region,region_id,region_activity_center)%>%
  mutate(
    rgb = scale_to_rgb(subregion_activity_center, bounds = c(-2, 2), format = "comma"))

#Change dataframe to have bed 9 rgb file that can be upload on UCSC
final_region_subregion_rgb=final_region_subregion_rgb[order(final_region_subregion_rgb$region_id,final_region_subregion_rgb$start_subregion),] %>%
  mutate(strand="+",
         start_sub=start_subregion,
         end_sub=end_subregion,
         subregion_id=1) %>%
  select(chr,start_subregion,end_subregion,region_id,subregion_id,strand,start_sub,end_sub,rgb)


#Reorder dataframe and sort by region_id
# col_order=c("chr","start_subregion","end_subregion","coverage_cdna_region","coverage_lib_region","subregion_id","subregion_activity","start_region","end_region","region_type","region_id","region_activity")
# final_df <- final_df[, col_order]
# final_df=final_df[order(final_df$region_id,final_df$start_subregion),]


# subset of the dataframe with all the subregion
# Needed to conserve a dataframe with all the region
# Next step is to find region with core silencer
# Here we keep a df with all the value per region
full_region_activity=final_df_exclude %>%
  select(chr,start_region, end_region,region_id,region_activity_center,fpkm_cdna_region,fpkm_lib_region,region_type)
#Remove duplicates
full_region_activity = full_region_activity[!duplicated(full_region_activity[,c('chr','start_region', 'end_region','region_id','region_activity_center','fpkm_cdna_region','fpkm_lib_region','region_type')]),]

#Find core silencer region
#On all data
#Make groups. If value >1 increase of 1 grp ID
#Merge data with same grp ID in order to get the core silencer and the edge silencer
final_df.core_silencer=final_df_exclude %>%
  group_by(chr,region_id,start_region,end_region,region_type,region_activity_center,fpkm_cdna_region,fpkm_lib_region,gr = data.table::rleid(subregion_activity_center <= silencer_threshold)) %>%
  filter(subregion_activity_center <= silencer_threshold) %>%
  summarise(start_core_silencer = first(start_subregion),
            end_core_silencer = last(end_subregion),
            activity_core_silencer = mean(subregion_activity_center),
            silencer_edge_activity = min(subregion_activity_center),
            start_edge_silencer = min(start_subregion[subregion_activity_center == silencer_edge_activity]),
            end_edge_silencer = max(end_subregion[subregion_activity_center == silencer_edge_activity])) %>%
  select(-gr)%>%
  select(chr,start_region,end_region,fpkm_cdna_region,fpkm_lib_region,region_type,region_id,region_activity_center,start_core_silencer,end_core_silencer,activity_core_silencer,start_edge_silencer,end_edge_silencer,silencer_edge_activity)

#Find core enhancer region
#On all data
#Make groups. If value >1 increase of 1 grp ID
#Merge data with same grp ID in order to get the core silencer and the edge silencer
final_df.core_enhancer=final_df_exclude %>%
  group_by(chr,region_id,start_region,end_region,region_type,region_activity_center,fpkm_cdna_region,fpkm_lib_region,gr = data.table::rleid(subregion_activity_center >= enhancer_threshold)) %>%
  filter(subregion_activity_center >= enhancer_threshold) %>%
  summarise(start_core_enhancer = first(start_subregion),
            end_core_enhancer = last(end_subregion),
            activity_core_enhancer = mean(subregion_activity_center),
            enhancer_edge_activity = min(subregion_activity_center),
            start_edge_enhancer = min(start_subregion[subregion_activity_center == enhancer_edge_activity]),
            end_edge_enhancer = max(end_subregion[subregion_activity_center == enhancer_edge_activity])) %>%
  select(-gr)%>%
  select(chr,start_region,end_region,fpkm_cdna_region,fpkm_lib_region,region_type,region_id,region_activity_center,start_core_enhancer,end_core_enhancer,activity_core_enhancer,start_edge_enhancer,end_edge_enhancer,enhancer_edge_activity)

#merge core enhancer with full df

full_df_enhancer=merge(final_df.core_enhancer, full_region_activity, by=c("chr","start_region", "end_region","region_id","region_activity_center","fpkm_cdna_region","fpkm_lib_region","region_type"),all = TRUE) # NA's match

#Merge the dataframe of region which have at least one core silencer with all the region df
full_df=merge(final_df.core_silencer, full_df_enhancer, by=c("chr","start_region", "end_region","region_id","region_activity_center","fpkm_cdna_region","fpkm_lib_region","region_type"),all = TRUE) # NA's match

#In very few cases (1 or 2) because of normalization some region which have activity very close to -1 (-1.002) -> Then they don't have core silencer (very close to -1)
full_df=full_df %>% 
  mutate(start_core_silencer = ifelse(region_activity_center <= silencer_threshold & is.na(start_core_silencer),start_region,start_core_silencer),
         end_core_silencer = ifelse(region_activity_center <= silencer_threshold & is.na(end_core_silencer),end_region,end_core_silencer),
         activity_core_silencer = ifelse(region_activity_center <= silencer_threshold & is.na(activity_core_silencer),region_activity_center,activity_core_silencer),
         start_edge_silencer = ifelse(region_activity_center <= silencer_threshold & is.na(start_edge_silencer),start_region,start_edge_silencer),
         end_edge_silencer = ifelse(region_activity_center <= silencer_threshold & is.na(end_edge_silencer),end_region,end_edge_silencer),
         silencer_edge_activity = ifelse(region_activity_center <= silencer_threshold & is.na(silencer_edge_activity),region_activity_center,silencer_edge_activity)
  )


#Normalize the reads count cdna/lib using rpkm
# sum_count_full_df=apply(full_df[,c("coverage_cdna_region","coverage_lib_region")],2,sum)
# region_list_clones_full_df <- lapply(region_list, function(x) which(full_df$region_id == x))
# 
# #Get RPKM 
# normalized_cov_df_activity_full_df=get_activities_normalized_coverage(region_list_clones_full_df,full_df)
# 
# #Clean df
# region_list_full_df <- levels(full_df$region_id)
# 
# full_df.region_indexes = match(full_df[,"region_id"], region_list_full_df)
# full_df$region_activity <- NULL
# full_df$coverage_cdna_region <- NULL
# full_df$coverage_lib_region <- NULL
# 
# #Add activity / RPKM cDNA / RPKM library to the final df
# full_df$fpkm_cdna=normalized_cov_df_activity_full_df$fpkm_cdna[full_df.region_indexes]
# full_df$fpkm_lib=normalized_cov_df_activity_full_df$fpkm_lib[full_df.region_indexes]
# full_df$region_activity=log2(normalized_cov_df_activity_full_df$region_activity[full_df.region_indexes])

# Annotate and Reorder df
# full_df=full_df %>%
# select(chr,start_region,end_region,fpkm_cdna,fpkm_lib,region_type,region_id,region_activity,start_core_silencer,end_core_silencer,activity_core_silencer,start_edge,end_edge,activity_edge)

# Annotate the df with gene
# all is FALSE because we exclude some region (fpkm<1) : to avoid NA value
full_df_annotated=merge(full_df,DHS_annotated, by=c("region_id","chr","region_type"),all=FALSE)%>%
  select(chr,start_region,end_region,region_id,start_dhs,end_dhs,
         fpkm_cdna_region,fpkm_lib_region,region_type,region_activity_center,
         start_core_silencer,end_core_silencer,activity_core_silencer,
         start_edge_silencer,end_edge_silencer,silencer_edge_activity,
         start_core_enhancer,end_core_enhancer,activity_core_enhancer,
         start_edge_enhancer,end_edge_enhancer,enhancer_edge_activity,
         genes)

#Split dataframe to got only region with activity lower than -1/ core silencer
# Not need if I create file with region and subregion in rgb

full_df_annotated_region_silencer=full_df_annotated %>%
  filter(region_activity_center <= silencer_threshold)  %>%
  select(c(chr,start_region,end_region,region_id,start_dhs,end_dhs,fpkm_cdna_region,fpkm_lib_region,region_type,region_activity_center,genes))
  full_df_annotated_region_silencer = full_df_annotated_region_silencer[!duplicated(full_df_annotated_region_silencer[,c('region_id')]),]

#Compute the inflexion point : threshold to determine the enhancer

th_inflexion_point_enhancer <- log2(calculate_cutoff(2^full_df_annotated$region_activity_center, drawPlot=F)$absolute)

calculate_cutoff(dat[idx,5], drawPlot=F)$absolute
#Get the enhancer region
full_df_annotated_region_enhancer=full_df_annotated %>%
    filter(region_activity_center >= th_inflexion_point_enhancer)  %>%
    select(c(chr,start_region,end_region,region_id,start_dhs,end_dhs,fpkm_cdna_region,fpkm_lib_region,region_type,region_activity_center,genes))
full_df_annotated_region_enhancer = full_df_annotated_region_enhancer[!duplicated(full_df_annotated_region_enhancer[,c('region_id')]),]
  
#Get the core silencer which are in silencer region
final_df.core_silencer_annotated=merge(final_df.core_silencer,DHS_annotated, by=c("region_id","chr","region_type"),all=FALSE)%>%
  filter(region_activity_center <= silencer_threshold)  %>%
  select(chr,start_core_silencer,end_core_silencer,activity_core_silencer,
         start_edge_silencer,end_edge_silencer,silencer_edge_activity,
         region_type,region_id,region_activity_center,genes,start_dhs,end_dhs)

#Get the core enhancer which are in enhancer region
final_df.core_enhancer_annotated=merge(final_df.core_enhancer,DHS_annotated, by=c("region_id","chr","region_type"),all=FALSE)%>%
  filter(region_activity_center >= th_inflexion_point_enhancer)  %>%
  select(chr,start_core_enhancer,end_core_enhancer,activity_core_enhancer,
         start_edge_enhancer,end_edge_enhancer,enhancer_edge_activity,
         region_type,region_id,region_activity_center,genes,start_dhs,end_dhs)


#Get weak and strong silencer
# weak_silencer=filter(final_df.core_silencer_annotated, activity_core_silencer > -2)
# strong_silencer=filter(final_df.core_silencer_annotated, activity_core_silencer <= -2)


#Get df with only edge
# final_df.edge_annotated=merge(final_df.core_silencer,DHS_annotated, by=c("region_id","chr","region_type"),all=FALSE)%>%
#   filter(region_activity_center <= -1)  %>%
#   select(chr,start_edge,end_edge,activity_edge,region_type,region_id,genes,start_dhs,end_dhs)


#####################################
############### OUTPUT############### 
#####################################

# Output for the summary df : contain all the informations
write.table(full_df_annotated, file=args[3],col.names = TRUE,row.names=FALSE,sep="\t",quote=FALSE)

#Output of region 
write.table(full_df_annotated_region_silencer, file=args[4],col.names = TRUE,row.names=FALSE,sep="\t",quote=FALSE)

write.table(full_df_annotated_region_enhancer, file=args[5],col.names = TRUE,row.names=FALSE,sep="\t",quote=FALSE)

#RGB file
# rgb_file<-file(args[6])
cat(header_rgb, file =args[5])
fwrite(x = final_region_subregion_rgb,
       file = args[4],
       sep = "\t",
       col.names=F,
       append=T)


# Bed file of core silencer
write.table(final_df.core_silencer_annotated, file=args[6],col.names = TRUE,row.names=FALSE,sep="\t",quote=FALSE)
# Bed file of weak core silencer
# write.table(weak_silencer, file=args[6],col.names = TRUE,row.names=FALSE,sep="\t",quote=FALSE)
# # Bed file of strong core silencer
# write.table(strong_silencer, file=args[7],col.names = TRUE,row.names=FALSE,sep="\t",quote=FALSE)
# 
# # Bed file of edge silencer
# write.table(final_df.edge_annotated, file=args[8],col.names = TRUE,row.names=FALSE,sep="\t",quote=FALSE)

# Save environment : if need to check something on the data on the R processing
save.image(file=args[7])

#####################
#Test on some regions
#####################
# test=data.frame(
#   start=c(1,50,100,150,200,250,300,350,400,450,500,100,150,200),
#   end=c(50,100,150,200,250,300,350,400,450,500,550,150,200,250),
#   sub_id=c("subid_1_1","subid_1_2","subid_1_3","subid_1_4","subid_1_5","subid_1_6","subid_1_7","subid_1_8","subid_1_9","subid_1_10","subid_1_11","subid_2_1","subid_2_2","subid_2_3"),
#   id=c(rep("id_1",11),rep("id_2",3)),
#   activity=c(-0.2,-0.6,-1,-1.2,-1.6,-1.6,-0.5,0.2,-1.2,-0.8,0.1,0.1,1.2,0.5))
# 
# library(dplyr)
# 
# threshold <- -1
# 
# test.group <- test %>%
#   mutate(grp = cumsum(activity > threshold)) 
# 
# print(test.group, row.names = F)
# 
# test.result <- 
#   test.group %>%
#   subset(activity <= -1) %>%
#   group_by(id, grp) %>%
#   arrange(activity) %>%
#   summarise(
#     start.min = first(start),
#     end.min = first(end),
#     sub_id.min = first(sub_id),
#     activity.min = first(activity),
#     start = min(start), 
#     end = max(end), 
#     activity = mean(activity)
#   ) %>% 
#   ungroup() %>%
#   select(start, end, id, activity, start.min, end.min, sub_id.min, activity.min)
# 
# 
# region_3313=final_df[final_df$region_id=='region_3313',]

# region_3313.test=region_3313 %>%
#   group_by(chr,start_region,end_region,region_type,region_activity,coverage_cdna_region,coverage_lib_region,gr = data.table::rleid(subregion_activity <= -1)) %>%
#   filter(subregion_activity <= -1) %>%
#   summarise(start_core_silencer = first(start_subregion),
#             end_core_silencer = last(end_subregion),
#             region_id = first(region_id),
#             activity_core_silencer = mean(subregion_activity),
#             activity_edge = min(subregion_activity),
#             start_edge = min(start_subregion[subregion_activity == activity_edge]),
#             end_edge = max(end_subregion[subregion_activity == activity_edge])) %>%
#   select(-gr)%>%
#   select(chr,start_region,end_region,coverage_cdna_region,coverage_lib_region,region_type,region_id,region_activity,start_core_silencer,end_core_silencer,activity_core_silencer,start_edge,end_edge,activity_edge)

# 
# region_3313_min=subset(region_3313, subregion_activity==min(subregion_activity, na.rm=TRUE))
# region_5412_sil(test, subregion_activity==min(subregion_activity, na.rm=TRUE)) # Get the lower activity rows

#Compute coverage by region
# region_3330=input_clone_dataf[input_clone_dataf$region_id=='region_3330',]
# region_3330_3313 <- input_clone_dataf[input_clone_dataf$region_id %in% c("region_3330", "region_3313"), ]
# 
# coverage_cdna_region_3330_3313 = sapply( region_3330_3313, function(x) sum( df[!duplicated(df[,c('chr','start_clone','end_clone')]),][x, "cdna_count"] ))
# 
# region_3330_3313 %>% 
#   group_by(chr,start_clone,end_clone,region_id) %>% 
#   summarise_all(sum(cdna_count))
# 
# uniq_region_3330=region_3330[!duplicated(region_3330[,c('chr','start_clone','end_clone')]),]
# 
# coverage_region_3330 = sum( region_3330[!duplicated(region_3330[,c('chr','start_clone','end_clone')]),][, "cdna_count"] )
# 
# region_activity_df=get_activities(region_list_clones,input_clone_dataf)
# 
# https://stackoverflow.com/questions/43077633/sum-values-in-different-rows-sharing-the-same-value-in-a-column

#Identify core silencer
# region_5412=final_df[final_df$region_id=='region_5412',]
# 
# region_5412.test=region_5412 %>%
#   group_by(gr = data.table::rleid(subregion_activity <= -1)) %>%
#   filter(subregion_activity <= -1) %>%
#   summarise(start_core_silencer = first(start_subregion), 
#             end_core_silencer = last(end_subregion), 
#             region_id = first(region_id),
#             activity_core_silencer = mean(subregion_activity), 
#             activity_edge = min(subregion_activity), 
#             start_edge = min(start_subregion[subregion_activity == activity_edge]), 
#             end_edge = max(end_subregion[subregion_activity == activity_edge])) %>%
#   select(-gr)
# region_5412=region_5412[order(region_5412$start_subregion),]
# 
# region_5412.group <- region_5412 %>%
#   mutate(grp = cumsum(subregion_activity > -1))
# 
# region_5412.result <-
#   region_5412.group %>%
#   subset(subregion_activity <= -1) %>%
#   group_by(chr,region_id, grp,start_region,end_region,region_type,region_activity) %>%
#   summarise(start_subregion = min(start_subregion), end_subregion = max(end_subregion),
#             core_silencer_activity = mean(subregion_activity)) %>% ungroup() %>%
#   select(chr,start_subregion,end_subregion,core_silencer_activity,start_region,end_region,region_type,region_id,region_activity)

#Example of region with 2 sub regions silencer
# region_3373=final_df[final_df$region_id=='region_3373',]
# 
# region_3373_test=region_3373 %>%
#   group_by(chr,start_region,end_region,region_type,region_activity,coverage_cdna_region,coverage_lib_region,gr = data.table::rleid(subregion_activity <= -1)) %>%
#   filter(subregion_activity <= -1) %>%
#   summarise(start_core_silencer = first(start_subregion),
#             end_core_silencer = last(end_subregion),
#             region_id = first(region_id),
#             activity_core_silencer = mean(subregion_activity),
#             activity_edge = min(subregion_activity),
#             start_edge = min(start_subregion[subregion_activity == activity_edge]),
#             end_edge = max(end_subregion[subregion_activity == activity_edge])) %>%
#   select(-gr)%>%
#   select(chr,start_region,end_region,coverage_cdna_region,coverage_lib_region,region_type,region_id,region_activity,start_core_silencer,end_core_silencer,activity_core_silencer,start_edge,end_edge,activity_edge)
# 
# # 
# region_3373.group <- region_3373 %>%
#   mutate(grp = cumsum(subregion_activity > -1))
# 
# region_3373.result <- 
#   region_3373.group %>%
#   subset(subregion_activity <= -1) %>%
#   group_by(chr,region_id, grp,start_region,end_region,region_type,region_activity) %>%
#   summarise(start_core_silencer = min(start_subregion), end_core_subregion = max(end_subregion),
#             core_silencer_activity = mean(subregion_activity),start_edge=start_subregion[,min(subregion_activity)],end_edge=end_subregion[,min(subregion_activity)],activity_edge=min(subregion_activity)) %>% ungroup() %>%
#   select(chr,start_core_silencer,end_core_subregion,core_silencer_activity,start_edge,end_edge,activity_edge,start_region,end_region,region_type,region_id,region_activity)

# region_3373.result <-
#   region_3373.group %>%
#   subset(subregion_activity <= -1) %>%
#   group_by(chr,region_id, grp,start_region,end_region,region_type,region_activity) %>%
#   arrange(subregion_activity) %>%
#   summarise(
#     start_edge = first(start_subregion),
#     end_edge = first(end_subregion),
#     edge_silencer_id = first(subregion_id),
#     activity_edge = first(subregion_activity),
#     start_core_silencer = min(start_subregion),
#     end_core_silencer = max(end_subregion),
#     activity_core_silencer = mean(subregion_activity)
#   ) %>%
#   ungroup() %>%
#   select(chr,start_region,end_region,region_type,region_id,region_activity,start_core_silencer,end_core_silencer,activity_core_silencer,start_edge,end_edge,activity_edge,edge_silencer_id)

# local_min=which(diff(sign(diff(final_df$subregion_activity[final_df$region_id=='region_5412'])))==-2)+1
# local_min=which(diff(sign(diff(region_5412$subregion_activity)))>0)+1

# region_5412_sil=region_5412[region_5412$subregion_activity <= -1, ] 
# 
# region_5412_sil(test, subregion_activity==min(subregion_activity, na.rm=TRUE)) # Get the lower activity rows

# write.table(region_5412_sil, file="sacapus_remote/silencer_project/subregion_analysis/region_5412_sil.bed",col.names = TRUE,row.names=FALSE,sep="\t",quote=FALSE)
