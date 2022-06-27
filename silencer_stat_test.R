##############################################
##### CapSTARR-seq Statistical Analysis ######
##############################################


### Upload Source and Library ###

# library("SciViews")
library("ggplot2")
# library("quantreg")
library("dplyr")
# library("reshape")
library(ggpubr)
library(data.table)
source(file = "/home/sadouni/script/R/capstarr_functions.R")

###################
#### Load Data ####
###################

##### PGK ######
pPGK_clone=as.data.frame(read.table("/home/sadouni/TAGC/IRCAN/PGK/clone_list.1bp.PGK.all_counts.DHS.Peaks.Final.bed",sep="\t", dec=".",fill=TRUE,header = FALSE,col.names=(c("chr","start","end","point","one","strand","r1_count","r2_count","r3_count","lib_count","region_type","region_id"))))

#Input_clone_list will be used to compute params(mean) of negative binomial
#input_clone_list=as.data.frame(read.table(file="/home/sadouni/TAGC/IRCAN/PGK/negativ_binom_test/input_clone_list.divided_in_half.bed",sep='\t',header=FALSE, skip=1))

#Activity threeshold used to compute p-value
activity_cutoff=0.5 


#Overall count : Make the sum of different column
sum_count=apply(pPGK_clone[,c("r1_count","r2_count","r3_count","lib_count")],2,sum)
total_reads_all_reps_together<-sum(pPGK_clone$r1_count,pPGK_clone$r2_count,pPGK_clone$r3_count)

#Make a List of region ID
#region_list contain all region ID
region_list <- levels(pPGK_clone$region_id)

#List of all unique region ID which contain all the clone that belong this region
region_list_clones <- lapply(region_list, function(x) which(pPGK_clone$region_id == x))



# Parse bash result
ids = read.table("/home/sadouni/TAGC/IRCAN/PGK/tmp/region_id_line_number_dict.txt")




#Create background model which will be the random region + region not containing inside a peak or DHS
background_region_list=lapply(1:length(region_list),function(x){
sample(get_random_regions(pPGK_clone),size=length(region_list_clones[[x]]),replace=FALSE)
})

##########################
#### Regions activity ####
##########################
region_activity=get_activities(region_list_clones,pPGK_clone)$region_activity
background_activities=get_activities(background_region_list,pPGK_clone)$region_activity

region_activity_fdr=sapply(region_activity, function(x){
  mean(background_activities<x,na.rm = TRUE)/mean(region_activity<x,na.rm = TRUE)
})

#Used to make violin plot of DHS activity VS Random regions
# clone_activities_df=data.frame(clone=region_activity,region_type="DHS")
# background_activities_df=data.frame(clone=background_activities,region_type="Background")
# DHS_vs_Background=rbind(clone_activities_df,background_activities_df)
# DHS_vs_Background$log2_clone=log2(DHS_vs_Background$clone)
# 
# my_comparisons <- list(c("DHS", "Background"))
# (violin_dhs_background<-
#   ggplot(DHS_vs_Background, aes(x=region_type, y=log2_clone))+
#   geom_violin(scale = "width",adjust = .5,fill='#A4A4A4', color="darkred")+
#   geom_boxplot(width=0.1,outlier.shape = NA) + theme_gray()+
#   stat_summary(fun.y=mean, geom="point",size=1,color="red",aes(shape="Mean")) +
#   stat_compare_means(comparisons = my_comparisons,method = "wilcox.test")+ # Add stars comparisons p-value
#   theme_grey()+
#   labs(title="DHS vs Random region",x="region type",y="Log2(regions activity)"))


#################################
### Statistical Analysis using ##
###### Negative Binomiam  #######
##########  pvals  ##############
#################################

#Create a dataframe for each library count the correct mean and var based on the second half of the library
#df_half_data_mean_with_NA=create_means_var_table(input_df = input_clone_list,repcolA = 7)

#Input count col of pPGK dataframe = 10
mean_var_table=create_means_var_table_rpois(pPGK_clone,datacolA = 10)

# i. First we calcul the p-values for each replicates of each clone
#    based on hypothesis that the region have no activity (=1)
#	   so the expected number of 'hits' in the replicate = activity x clone_library_reads = 1 x normalization value x library reads

#For the value which have NA fill with the params of linear regression
clone_pvals_null_hyp.rep1=calc_pvals(region_list_clones,test_activity=1,replicate=1)
clone_pvals_null_hyp.rep2=calc_pvals(region_list_clones,test_activity=1,replicate=2)
clone_pvals_null_hyp.rep3=calc_pvals(region_list_clones,test_activity=1,replicate=3)

background_pvals_null_hyp.rep1=calc_pvals(background_region_list,test_activity=1,replicate=1)
background_pvals_null_hyp.rep2=calc_pvals(background_region_list,test_activity=1,replicate=2)
background_pvals_null_hyp.rep3=calc_pvals(background_region_list,test_activity=1,replicate=3)


# ii. Second hypothesis : clone activity is equal to observed region activity
# 	  so the lambda = region_activity x clone_library_reads

# Vector of region indexes, the same length as the vector of clones
#pPGK_clone.region_indexes = match(pPGK_clone[,"region_id"], region_list)
#Calcul of pval based on observed activity (here activities were calculated with all replicates combined)

# clone_pvals_obs.rep1=calc_pvals(region_list_clones,test_activity=region_activity,replicate=1)
# clone_pvals_obs.rep2=calc_pvals(region_list_clones,test_activity=region_activity,replicate=2)
# clone_pvals_obs.rep3=calc_pvals(region_list_clones,test_activity=region_activity,replicate=3)
# 
# background_pvals_obs.rep1=calc_pvals(background_region_list,test_activity=background_activities,replicate=1)
# background_pvals_obs.rep2=calc_pvals(background_region_list,test_activity=background_activities,replicate=2)
# background_pvals_obs.rep3=calc_pvals(background_region_list,test_activity=background_activities,replicate=3)


#Calcul of pval with cut off activity of 0.5

clone_pvals_cutoff.rep1=calc_pvals(region_list_clones,test_activity=activity_cutoff,replicate=1)
clone_pvals_cutoff.rep2=calc_pvals(region_list_clones,test_activity=activity_cutoff,replicate=2)
clone_pvals_cutoff.rep3=calc_pvals(region_list_clones,test_activity=activity_cutoff,replicate=3)

background_pvals_cutoff.rep1=calc_pvals(background_region_list,test_activity=activity_cutoff,replicate=1)
background_pvals_cutoff.rep2=calc_pvals(background_region_list,test_activity=activity_cutoff,replicate=2)
background_pvals_cutoff.rep3=calc_pvals(background_region_list,test_activity=activity_cutoff,replicate=3)



#### Calcul of Likelihood  ####
# Calcul of the 'likelihood' of each hypothesis at each region 
#	= product of all the pvalues for each hypothesis for each clone in the region

#Null hypothesis
###
region_likelihoods_null_hyp.rep1 = sapply(clone_pvals_null_hyp.rep1, function(x)  sum(log(x)))
region_likelihoods_null_hyp.rep2 = sapply(clone_pvals_null_hyp.rep2, function(x)  sum(log(x)))
region_likelihoods_null_hyp.rep3 = sapply(clone_pvals_null_hyp.rep3, function(x)  sum(log(x)))

background_region_likelihoods_null_hyp.rep1 = sapply(background_pvals_null_hyp.rep1, function(x)  sum(log(x)))
background_region_likelihoods_null_hyp.rep2 = sapply(background_pvals_null_hyp.rep2, function(x)  sum(log(x)))
background_region_likelihoods_null_hyp.rep3 = sapply(background_pvals_null_hyp.rep3, function(x)  sum(log(x)))

#obs hypothesis
# region_likelihoods_obs.rep1 = sapply(clone_pvals_obs.rep1, function(x)  sum(log(x)))
# region_likelihoods_obs.rep2 = sapply(clone_pvals_obs.rep2, function(x)  sum(log(x)))
# region_likelihoods_obs.rep3 = sapply(clone_pvals_obs.rep3, function(x)  sum(log(x)))

# background_region_likelihoods_obs.rep1 = sapply(background_pvals_obs.rep1, function(x)  sum(log(x)))
# background_region_likelihoods_obs.rep2 = sapply(background_pvals_obs.rep2, function(x)  sum(log(x)))
# background_region_likelihoods_obs.rep3 = sapply(background_pvals_obs.rep3, function(x)  sum(log(x)))

# #obs hypothesis based on the activities per replicates
# region_likelihoods_obs.rep1.actr1 = sapply( region_list_clones, function(x)  10^(sum(log10( clone_pvals_obs.rep1.actr1[ x ]), na.rm=TRUE ) ))
# region_likelihoods_obs.rep2.actr2 = sapply( region_list_clones, function(x)  10^(sum(log10( clone_pvals_obs.rep2.actr2[ x ]), na.rm=TRUE ) ))
# region_likelihoods_obs.rep3.actr3 = sapply( region_list_clones, function(x)  10^(sum(log10( clone_pvals_obs.rep3.actr3[ x ]), na.rm=TRUE ) ))

#based on activity cutoff 0.5
region_likelihoods_cutoff.rep1 = sapply(clone_pvals_cutoff.rep1, function(x)  sum(log(x)))
region_likelihoods_cutoff.rep2 = sapply(clone_pvals_cutoff.rep2, function(x)  sum(log(x)))
region_likelihoods_cutoff.rep3 = sapply(clone_pvals_cutoff.rep3, function(x)  sum(log(x)))

background_region_likelihoods_cutoff.rep1 = sapply(background_pvals_cutoff.rep1, function(x)  sum(log(x)))
background_region_likelihoods_cutoff.rep2 = sapply(background_pvals_cutoff.rep2, function(x)  sum(log(x)))
background_region_likelihoods_cutoff.rep3 = sapply(background_pvals_cutoff.rep3, function(x)  sum(log(x)))

## Based on strand 

region_likelihoods_null_hyp.rep1.plus = compute_likelihood_strand(clone_pvals_null_hyp.rep1,df = pPGK_clone,strand = "+" )
region_likelihoods_null_hyp.rep1.minus = compute_likelihood_strand(clone_pvals_null_hyp.rep1,df = pPGK_clone,strand = "-" )

region_likelihoods_null_hyp.rep2.plus = compute_likelihood_strand(clone_pvals_null_hyp.rep2,df = pPGK_clone,strand = "+" )
region_likelihoods_null_hyp.rep2.minus = compute_likelihood_strand(clone_pvals_null_hyp.rep2,df = pPGK_clone,strand = "-" )

region_likelihoods_null_hyp.rep3.plus = compute_likelihood_strand(clone_pvals_null_hyp.rep3,df = pPGK_clone,strand = "+" )
region_likelihoods_null_hyp.rep3.minus = compute_likelihood_strand(clone_pvals_null_hyp.rep3,df = pPGK_clone,strand = "-" )

background_region_likelihoods_null_hyp.rep1.plus = compute_likelihood_strand(background_pvals_null_hyp.rep1,df = pPGK_clone,strand = "+" )
background_region_likelihoods_null_hyp.rep1.minus = compute_likelihood_strand(background_pvals_null_hyp.rep1,df = pPGK_clone,strand = "-" )

background_region_likelihoods_null_hyp.rep2.plus = compute_likelihood_strand(background_pvals_null_hyp.rep2,df = pPGK_clone,strand = "+" )
background_region_likelihoods_null_hyp.rep2.minus = compute_likelihood_strand(background_pvals_null_hyp.rep2,df = pPGK_clone,strand = "-" )

background_region_likelihoods_null_hyp.rep3.plus = compute_likelihood_strand(background_pvals_null_hyp.rep3,df = pPGK_clone,strand = "+" )
background_region_likelihoods_null_hyp.rep3.minus = compute_likelihood_strand(background_pvals_null_hyp.rep3,df = pPGK_clone,strand = "-" )


#Cut-off

region_likelihoods_cutoff.rep1.plus = compute_likelihood_strand(clone_pvals_cutoff.rep1,df = pPGK_clone,strand = "+" )
region_likelihoods_cutoff.rep1.minus = compute_likelihood_strand(clone_pvals_cutoff.rep1,df = pPGK_clone,strand = "-" )

region_likelihoods_cutoff.rep2.plus = compute_likelihood_strand(clone_pvals_cutoff.rep2,df = pPGK_clone,strand = "+" )
region_likelihoods_cutoff.rep2.minus = compute_likelihood_strand(clone_pvals_cutoff.rep2,df = pPGK_clone,strand = "-" )

region_likelihoods_cutoff.rep3.plus = compute_likelihood_strand(clone_pvals_cutoff.rep3,df = pPGK_clone,strand = "+" )
region_likelihoods_cutoff.rep3.minus = compute_likelihood_strand(clone_pvals_cutoff.rep3,df = pPGK_clone,strand = "-" )

background_region_likelihoods_cutoff.rep1.plus = compute_likelihood_strand(background_pvals_cutoff.rep1,df = pPGK_clone,strand = "+" )
background_region_likelihoods_cutoff.rep1.minus = compute_likelihood_strand(background_pvals_cutoff.rep1,df = pPGK_clone,strand = "-" )

background_region_likelihoods_cutoff.rep2.plus = compute_likelihood_strand(background_pvals_cutoff.rep2,df = pPGK_clone,strand = "+" )
background_region_likelihoods_cutoff.rep2.minus = compute_likelihood_strand(background_pvals_cutoff.rep2,df = pPGK_clone,strand = "-" )

background_region_likelihoods_cutoff.rep3.plus = compute_likelihood_strand(background_pvals_cutoff.rep3,df = pPGK_clone,strand = "+" )
background_region_likelihoods_cutoff.rep3.minus = compute_likelihood_strand(background_pvals_cutoff.rep3,df = pPGK_clone,strand = "-" )


##############################################
## Chi-squared test on log likelihood ratio ##
###############################################
#	Final p-value for the silencer activity of each region

pchisq_rep1=pchisq(2*((region_likelihoods_obs.rep1)-(region_likelihoods_null_hyp.rep1)),1,lower.tail = FALSE)
pchisq_rep2=pchisq(2*((region_likelihoods_obs.rep2)-(region_likelihoods_null_hyp.rep2)),1,lower.tail = FALSE)
pchisq_rep3=pchisq(2*((region_likelihoods_obs.rep3)-(region_likelihoods_null_hyp.rep3)),1,lower.tail = FALSE)

#pchisq_rep1_less_than=ifelse(region_activity<1,pchisq_rep1,1)

pchisq_background_rep1=pchisq(2*((background_region_likelihoods_obs.rep1)-(background_region_likelihoods_null_hyp.rep1)),1,lower.tail = FALSE)
pchisq_background_rep2=pchisq(2*((background_region_likelihoods_obs.rep2)-(background_region_likelihoods_null_hyp.rep2)),1,lower.tail = FALSE)
pchisq_background_rep3=pchisq(2*((background_region_likelihoods_obs.rep3)-(background_region_likelihoods_null_hyp.rep3)),1,lower.tail = FALSE)

pchisq_fdr_rep1=sapply(pchisq_rep1, function(x){
  mean(pchisq_background_rep1<=x,na.rm = TRUE)/mean(pchisq_rep1<=x,na.rm = TRUE)
})
#pchisquare with activity cut off 0.5 and activity region calculated per replicates

pchisq_rep1_silencer= compute_pchisq_silencer(region_likelihoods_obs.rep1,region_likelihoods_cutoff.rep1,region_type = region_activity) 
pchisq_rep2_silencer= compute_pchisq_silencer(region_likelihoods_obs.rep2,region_likelihoods_cutoff.rep2,region_type = region_activity) 
pchisq_rep3_silencer= compute_pchisq_silencer(region_likelihoods_obs.rep3,region_likelihoods_cutoff.rep3,region_type = region_activity) 

pchisq_rep1_silencer_bg= compute_pchisq_silencer(background_region_likelihoods_obs.rep1,background_region_likelihoods_cutoff.rep1,region_type = background_activities) 
pchisq_rep2_silencer_bg= compute_pchisq_silencer(background_region_likelihoods_obs.rep2,background_region_likelihoods_cutoff.rep2,region_type = background_activities) 
pchisq_rep3_silencer_bg= compute_pchisq_silencer(background_region_likelihoods_obs.rep3,background_region_likelihoods_cutoff.rep3,region_type = background_activities) 

pchisq_fdr_rep1_silencer=compute_fdr(pchisq_rep1_silencer,pchisq_rep1_silencer_bg)
pchisq_fdr_rep2_silencer=compute_fdr(pchisq_rep2_silencer,pchisq_rep2_silencer_bg)
pchisq_fdr_rep3_silencer=compute_fdr(pchisq_rep3_silencer,pchisq_rep3_silencer_bg)


# pchisq_rep1_cutoff=pchisq(2*((region_likelihoods_obs.rep1)-(region_likelihoods_cutoff.rep1)),1,lower.tail = FALSE)
# pchisq_rep2_cutoff=pchisq(2*((region_likelihoods_obs.rep2)-(region_likelihoods_cutoff.rep2)),1,lower.tail = FALSE)
# pchisq_rep3_cutoff=pchisq(2*((region_likelihoods_obs.rep3)-(region_likelihoods_cutoff.rep3)),1,lower.tail = FALSE)
# 
# 
# pchisq_background_rep1_cutoff=pchisq(2*((background_region_likelihoods_obs.rep1)-(background_region_likelihoods_cutoff.rep1)),1,lower.tail = FALSE)
# pchisq_background_rep2_cutoff=pchisq(2*((background_region_likelihoods_obs.rep2)-(background_region_likelihoods_cutoff.rep2)),1,lower.tail = FALSE)
# pchisq_background_rep3_cutoff=pchisq(2*((background_region_likelihoods_obs.rep3)-(background_region_likelihoods_cutoff.rep3)),1,lower.tail = FALSE)

# pchisq_rep1_cutoff_df=data.frame(pchisq=-log(pchisq_rep1_cutoff),region_type="DHS")
# pchisq_background_rep1_cutoff_df=data.frame(pchisq=-log(pchisq_background_rep1_cutoff),region_type="Background")
# pchisq_DHS_vs_Background=rbind(pchisq_rep1_cutoff_df,pchisq_background_rep1_cutoff_df)
# 
# my_comparisons <- list(c("DHS", "Background"))
# (violin_dhs_background<-
#     ggplot(pchisq_DHS_vs_Background, aes(x=region_type, y=pchisq))+
#     geom_violin(scale = "width",adjust = .5,fill='#A4A4A4', color="darkred")+
#     #geom_boxplot(width=0.1,outlier.shape = NA) + theme_gray()+
#     #stat_summary(fun.y=mean, geom="point",size=1,color="red",aes(shape="Mean")) +
#     #stat_compare_means(comparisons = my_comparisons,method = "wilcox.test")+ # Add stars comparisons p-value
#     theme_grey()+
#     labs(title="pchisq : DHS vs Random region",x="region type",y="-log(pchisq))"))

### not used
# pchisq_rep1_cutoff_df_less_than=data.frame(pchisq=-log(pchisq_rep1_cutoff_less_than),region_type="DHS")
# pchisq_background_rep1_cutoff_df=data.frame(pchisq=-log(pchisq_background_rep1_cutoff_less_than),region_type="Background")
# pchisq_DHS_vs_Background=rbind(pchisq_rep1_cutoff_df,pchisq_background_rep1_cutoff_df)
# 
# my_comparisons <- list(c("DHS", "Background"))
# (violin_dhs_background<-
#     ggplot(pchisq_DHS_vs_Background, aes(x=region_type, y=pchisq))+
#     geom_violin(scale = "width",adjust = .5,fill='#A4A4A4', color="darkred")+
#     #geom_boxplot(width=0.1,outlier.shape = NA) + theme_gray()+
#     #stat_summary(fun.y=mean, geom="point",size=1,color="red",aes(shape="Mean")) +
#     #stat_compare_means(comparisons = my_comparisons,method = "wilcox.test")+ # Add stars comparisons p-value
#     theme_grey()+
#     labs(title="pchisq : DHS vs Random region",x="region type",y="-log(pchisq))"))
# 


# pchisq_rep2_cutoff_less_than=ifelse(region_activity<activity_cutoff,pchisq_rep2_cutoff,1)
# pchisq_background_rep2_cutoff_less_than=ifelse(background_activities<activity_cutoff,pchisq_background_rep2_cutoff,1)
# 
# pchisq_rep3_cutoff_less_than=ifelse(region_activity<activity_cutoff,pchisq_rep3_cutoff,1)
# pchisq_background_rep3_cutoff_less_than=ifelse(background_activities<activity_cutoff,pchisq_background_rep3_cutoff,1)
# 

# pchisq_fdr_rep1_cutoff_less_than=compute_fdr(pchisq_rep1_cutoff_less_than,pchisq_background_rep1_cutoff_less_than)
# pchisq_fdr_rep2_cutoff_less_than=compute_fdr(pchisq_rep2_cutoff_less_than,pchisq_background_rep2_cutoff_less_than)
# pchisq_fdr_rep2_cutoff_less_than=compute_fdr(pchisq_rep3_cutoff_less_than,pchisq_background_rep3_cutoff_less_than)



# plot(x = log2(region_activity),y = -log(pchisq_rep1_cutoff_less_than),ylim = c(0,50),xlim = c(-6,0))
# plot(x = log2(background_activities),y = -log(pchisq_background_rep1_cutoff_less_than),ylim = c(0,50),xlim = c(-6,0))
# 
# #Used to make violin plot of DHS activity VS Random regions
# pchisq_rep1_df=data.frame(pchisq=-log(pchisq_rep1),region_type="DHS")
# pchisq_background_rep1_df=data.frame(pchisq=-log(pchisq_background_rep1),region_type="Background")
# pchisq_DHS_vs_Background=rbind(pchisq_rep1_df,pchisq_background_rep1_df)
# 
# my_comparisons <- list(c("DHS", "Background"))
# (violin_dhs_background<-
#     ggplot(pchisq_DHS_vs_Background, aes(x=region_type, y=pchisq))+
#     geom_violin(scale = "width",adjust = .5,fill='#A4A4A4', color="darkred")+
#     #geom_boxplot(width=0.1,outlier.shape = NA) + theme_gray()+
#     #stat_summary(fun.y=mean, geom="point",size=1,color="red",aes(shape="Mean")) +
#     #stat_compare_means(comparisons = my_comparisons,method = "wilcox.test")+ # Add stars comparisons p-value
#     theme_grey()+
#     labs(title="pchisq : DHS vs Random region",x="region type",y="log(pchisq))"))
# 
# pchisq_rep2_df=data.frame(pchisq=pchisq_rep2,region_type="DHS")
# pchisq_background_rep2_df=data.frame(pchisq=pchisq_background_rep2,region_type="Background")
# pchisq_DHS_vs_Background_rep2=rbind(pchisq_rep2_df,pchisq_background_rep2_df)
# 
# my_comparisons <- list(c("DHS", "Background"))
# (violin_dhs_background<-
#     ggplot(pchisq_DHS_vs_Background_rep2, aes(x=region_type, y=pchisq))+
#     geom_violin(scale = "width",adjust = .5,fill='#A4A4A4', color="darkred")+
#     #geom_boxplot(width=0.1,outlier.shape = NA) + theme_gray()+
#     #stat_summary(fun.y=mean, geom="point",size=1,color="red",aes(shape="Mean")) +
#     stat_compare_means(comparisons = my_comparisons,method = "wilcox.test")+ # Add stars comparisons p-value
#     theme_grey()+
#     labs(title="pchisq : DHS vs Random region",x="region type",y="log(pchisq))"))

# This part just complete the dataframe with new data
# This will be usefull to create new webfiles to visualize data in UCSC genome browser

# New dataframe to not modify the original (can be swtich off and just put new information in the original dataframe)
pPGK=pPGK_clone


pPGK_clone.region_indexes = match(pPGK_clone[,"region_id"], region_list)


# Complete the dataframe with region activities and pchisq per replicates
pPGK$region_activity=log2(region_activity[pPGK_clone.region_indexes])
pPGK$sum_r1_count=get_activities(region_list_clones,pPGK)$region_totals_rep1[pPGK_clone.region_indexes]
pPGK$sum_r2_count=get_activities(region_list_clones,pPGK)$region_totals_rep2[pPGK_clone.region_indexes]
pPGK$sum_r3_count=get_activities(region_list_clones,pPGK)$region_totals_rep3[pPGK_clone.region_indexes]
pPGK$sum_lib_count=get_activities(region_list_clones,pPGK)$region_totals_lib[pPGK_clone.region_indexes]
pPGK$p_value_r1=pchisq_rep1_cutoff_less_than[pPGK_clone.region_indexes]
pPGK$fdr_r1=pchisq_fdr_rep1_cutoff_less_than[pPGK_clone.region_indexes]
pPGK$p_value_r2=pchisq_rep2_cutoff_less_than[pPGK_clone.region_indexes]
pPGK$fdr_r2=pchisq_fdr_rep2_cutoff_less_than[pPGK_clone.region_indexes]
pPGK$p_value_r3=pchisq_rep3_cutoff_less_than[pPGK_clone.region_indexes]
pPGK$fdr_r3=pchisq_fdr_rep3_cutoff_less_than[pPGK_clone.region_indexes]


pPGK_final=setDT(pPGK)[, .(chr=chr[1],start=start[1], end=end[.N],sum_r1_count=sum_r1_count[1],sum_r2_count=sum_r2_count[1],
                       sum_r3_count=sum_r3_count[1],
                       sum_lib_count=sum_lib_count[1],
                       region_type=region_type[1],
                       region_activity=region_activity[1],
                       p_value_r1=p_value_r1[1],fdr_r1=fdr_r1[1],
                       p_value_r2=p_value_r2[1],fdr_r2=fdr_r2[1],
                       p_value_r3=p_value_r3[1],fdr_r3=fdr_r3[1]),by=region_id]


##############################################################
##############################################################
##########       Select some specific candidates   ###########
##############################################################
##############################################################

top500_DHS_r1_pval=pPGK_final %>% filter((region_type == 'DNaseI'| region_type == 'CRM')&fdr_r1<=0.05)%>%
                                  top_n(-500, wt = p_value_r1)
top500_DHS_r2_pval=pPGK_final %>% filter((region_type == 'DNaseI'| region_type == 'CRM')&fdr_r2<=0.05)%>%
                                  top_n(-500, wt = p_value_r2)
top500_DHS_r3_pval=pPGK_final %>% filter((region_type == 'DNaseI'| region_type == 'CRM') &fdr_r3<=0.05)%>%
                                  top_n(-500, wt = p_value_r3)


top500_DHS_r1_pval=pPGK_final %>% filter(region_type == 'DNaseI'| region_type == 'CRM')%>%
                                          filter(fdr_r1<=0.05)

One_percent_fdr=pPGK_final %>% filter(fdr_r1 <=0.01)

Ten_percent_fdr=pPGK_final %>% filter(fdr_r1 <=0.1)

write.table(One_percent_fdr, file="/home/sadouni/TAGC/IRCAN/PGK/negativ_binom_test/pPGK.0.01fdr_r1.bed",col.names = FALSE,row.names=FALSE,sep="\t",quote=FALSE)

write.table(Ten_percent_fdr, file="/home/sadouni/TAGC/IRCAN/PGK/negativ_binom_test/pPGK.0.1fdr_r1.bed",col.names = FALSE,row.names=FALSE,sep="\t",quote=FALSE)

region_length=(One_percent_fdr$end)-(One_percent_fdr$start)
#Complete the dataframe with p chisq and the -log(pchisq) of each replicates
pPGK$pchisq_rep1=pchisq_rep1[ pPGK_clone.region_indexes ]
pPGK$minus_log_pchisq_rep1=-log(pPGK$pchisq_rep1)
pPGK$pchisq_rep2=pchisq_rep2[ pPGK_clone.region_indexes ]
pPGK$minus_log_pchisq_rep2=-log(pPGK$pchisq_rep2)
pPGK$pchisq_rep3=pchisq_rep3[ pPGK_clone.region_indexes ]
pPGK$minus_log_pchisq_rep3=-log(pPGK$pchisq_rep3)



#Complete the dataframe with pval of activity_cutoff
pPGK$pchisq_rep1_cutoff=pchisq_rep1_cutoff[ pPGK_clone.region_indexes ]
pPGK$minus_log_pchisq_cutoff_rep1=-log(pPGK$pchisq_rep1_cutoff)
pPGK$pchisq_rep2_cutoff=pchisq_rep2_cutoff[ pPGK_clone.region_indexes ]
pPGK$minus_log_pchisq_cutoff_rep2=-log(pPGK$pchisq_rep2_cutoff)
pPGK$pchisq_rep3_cutoff=pchisq_rep3_cutoff[ pPGK_clone.region_indexes ]
pPGK$minus_log_pchisq_cutoff_rep3=-log(pPGK$pchisq_rep3_cutoff)





# Browser track of probability

#Create file of -log(pchisq)
write.table(pPGK[,c("chr","start","end","minus_log_pchisq_rep1")], file="pPGK_silencer_probability.rep1.logp.bedGraph",col.names = FALSE,row.names=FALSE,sep="\t",quote=FALSE)
write.table(pPGK[,c("chr","start","end","minus_log_pchisq_rep2")], file="pPGK_silencer_probability.rep2.logp.bedGraph",col.names = FALSE,row.names=FALSE,sep="\t",quote=FALSE)
write.table(pPGK[,c("chr","start","end","minus_log_pchisq_rep3")], file="pPGK_silencer_probability.rep3.logp.bedGraph",col.names = FALSE,row.names=FALSE,sep="\t",quote=FALSE)

#Create file of -log(pchisq) with activity cut off
write.table(pPGK[,c("chr","start","end","minus_log_pchisq_cutoff_rep1")], file="pPGK_silencer_probability.cutoff_rep1.logp.bedGraph",col.names = FALSE,row.names=FALSE,sep="\t",quote=FALSE)
write.table(pPGK[,c("chr","start","end","minus_log_pchisq_cutoff_rep2")], file="pPGK_silencer_probability.cutoff_rep2.logp.bedGraph",col.names = FALSE,row.names=FALSE,sep="\t",quote=FALSE)
write.table(pPGK[,c("chr","start","end","minus_log_pchisq_cutoff_rep3")], file="pPGK_silencer_probability.cutoff_rep3.logp.bedGraph",col.names = FALSE,row.names=FALSE,sep="\t",quote=FALSE)

#Create a file of log(region_activity)
write.table(pPGK[,c("chr","start","end","log_region_activity")], file="pPGK_silencer_activity_rep1.logp.bedGraph",col.names = FALSE,row.names=FALSE,sep="\t",quote=FALSE)
write.table(pPGK[,c("chr","start","end","log_region_activity")], file="pPGK_silencer_activity_rep2.logp.bedGraph",col.names = FALSE,row.names=FALSE,sep="\t",quote=FALSE)
write.table(pPGK[,c("chr","start","end","log_region_activity")], file="pPGK_silencer_activity_rep3.logp.bedGraph",col.names = FALSE,row.names=FALSE,sep="\t",quote=FALSE)
