############################################
## Function used in capstarr-seq pipeline ##
############################################
# 07/10/2019 #
##############



#Create a dataframe which take in input the summary clone count dataframe
#repcolA is the first half of cDNA
#repcolB if is not specifiy will be repcolA +1 so second half of cDNA and you can also specify a specific number
create_means_var_table<-function(input_df,repcolA,repcolB=repcolA+1){
  table=data.frame(data_means = 
                     sapply(0:max(input_df[,repcolA],input_df[,repcolB]),function(y) mean(input_df[(input_df[,repcolA]==y),repcolB])),
                   range_mean=(0:max(input_df[,repcolA],input_df[,repcolB])))
  table$estimated_means=estimate_means(table)
  table$estimated_vars=estimate_vars(table)
  return(table)}

## more general version using entire data & synthesizing second replicate in silico
##datacolA -> input count
create_means_var_table_rpois<-function(input_df,datacolA) {
  dataA=input_df[,datacolA]
  dataB=sapply(dataA,function(x) mean(rpois(n=1,lambda=x))) #Maybe increase n and get the mean to have a better fit??
  table=data.frame(data_means = sapply(0:max(dataA,dataB),function(x) mean(dataB[dataA==x])),
                   range_mean=(0:max(dataA,dataB)))
  table$estimated_means=estimate_means(table)
  table$estimated_vars=estimate_vars(table)
  
  return(table)
}

#Fill NA in dataframe for each clone which have NA in the range_mean using linear regression
estimate_means<-function(x) {
  lm_result=lm(x$range_mean ~x$data_means)
  estimation_mean=ifelse(is.na(x$data_means), coef(lm_result)["(Intercept)"] +coef(lm_result)["x$data_means"]*x$range_mean,x$data_means)
  return(estimation_mean)}
#Compute Variancen for now var = 2* mean
estimate_vars<-function(x){
  return(2*x$estimated_mean)}

#Extract index of random region -> Will be used as background
get_random_regions<- function(df){
  return(which(pPGK_clone[ ,"region_type" ] == "Random" | pPGK_clone[ ,"region_type" ] == "."))}

#Return a dataframe containing the number of element per regions and their activities
get_activities=function(list_of_clones,df){
  summary_activities= data.frame(
    region_totals_rep1 = sapply( list_of_clones, function(x) sum( df[x, "r1_count"] )),
    region_totals_rep2 = sapply( list_of_clones, function(x) sum( df[x, "r2_count"] )),
    region_totals_rep3 = sapply( list_of_clones, function(x) sum( df[x, "r3_count"] )),
    region_totals_lib = sapply( list_of_clones, function(x) sum( df[x, "lib_count"] )))
  ## region_activity_repX is calculated using replicates separatly
  summary_activities$region_activity_rep1 = apply(summary_activities,1,function(x){(x[1]/x[4])* (sum_count["lib_count"] / sum_count["r1_count"])}) 
  summary_activities$region_activity_rep2 = apply(summary_activities,1,function(x){(x[2]/x[4])* (sum_count["lib_count"] / sum_count["r2_count"])}) 
  summary_activities$region_activity_rep3 = apply(summary_activities,1,function(x){(x[3]/x[4])* (sum_count["lib_count"] / sum_count["r3_count"])}) 
  # #calculate the ratio of the sums. Activity based on region_totals
  summary_activities$region_activity = apply(summary_activities,1, function(x){(sum(x[1:3])/x[4])*(sum_count["lib_count"]/total_reads_all_reps_together)})
  return(summary_activities)}

###Simplify the original function of negative binomial
mydnbinom=function(x,mean,var){dnbinom(x,mu=mean,size=mean^2/(var-mean))}
my_pval= function(x,mean,var){mydnbinom(x,mean,var)}

#Compute negative binomial p-values
#Take in input list of clones that belong the regions
#activity can be null =1 for null hypothesis / = to obs activty or decided cutoff
calc_pvals = function(index_list, test_activity=1,replicate) {
  if(length(test_activity)==1){test_activity=rep(test_activity,times=length(index_list))}
  sapply(1:length(index_list), function(x) {
    as.vector(sapply(index_list[x], function(y) {
      activity= test_activity[x]
      #get the activity of clones
      pval=my_pval( pPGK_clone[y,paste("r",replicate,"_count",sep = "")], 
                    mean=df_half_data_mean_with_NA$estimated_means[ pPGK_clone[y,"lib_count"]+1]*
                      (sum_count[paste("r",replicate,"_count",sep = "")]/sum_count["lib_count"])*activity,
                    var=df_half_data_mean_with_NA$estimated_vars[ pPGK_clone[y,"lib_count"]+1]*
                      (sum_count[paste("r",replicate,"_count",sep = "")]/sum_count["lib_count"])*activity)
      return(pval)
    }))
  })
}

compute_likelihood_strand=function(clone_p_value_list,df,strands=c("+","-")){
  sapply(1:length(clone_p_value_list), function(x)  sum(log(clone_p_value_list[[x]][df$strand[x]==strands])))
}

compute_pchisq_silencer=function(likelihood_obs,likelihood_cutoff,region_type=c(region_activity,background_activities)){
  ifelse(region_type<activity_cutoff,pchisq(2*((likelihood_obs)-(likelihood_cutoff)),1,lower.tail = FALSE),1)
}

compute_fdr=function(pchisq_list_active_clone,pchisq_list_background){
  sapply(pchisq_list_active_clone, function(x){
  mean(pchisq_list_background<=x,na.rm = TRUE)/mean(pchisq_list_active_clone<=x,na.rm = TRUE)
  })}
