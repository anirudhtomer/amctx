library(ggplot2)
library(doParallel)
registerDoParallel(cores = detectCores())

##############################################
# Load the data set and check missing data....no missing for this one
##############################################
amctx = read.csv2(file.choose(), header = T)
apply(amctx, MARGIN = 2, FUN= function(x){any(is.na(x))})

#remove subject 346 because according to the data the number of subjects per response is 
#log(pcr) = 238, log(creatinine) = 239. So they gotta be equal for the analysis

amctx = amctx[amctx$amctx!=346, ]

##############################################
# Do basic clean up
##############################################

#The zis number has one to one correspondance with amctx
amctx$zis = NULL
amctx$amctx = as.factor(amctx$amctx)

#The cols not of interest are the ones with medicinal usage 
colsOfInterest = c("amctx", "type", "measure", "value", "tx_s_days", "rec_age",
                   "gl_failure", "days_tx_gl", "rec_gender", "tx_previoustx", 
                   "d_age", "d_gender", "d_bmi", "rec_bmi", "tx_hla", "tx_pra",
                   "tx_dgf", "tx_cit", "tx_dm", "tx_hvdis", "rr_sys", "rr_dias", 
                   "rr_map", "ah_nr", "tx_dial_days", "d_type", "d_cadaveric")

amctx = amctx[, colsOfInterest]

amctx$years_tx_gl = amctx$days_tx_gl/365
amctx$tx_s_years = amctx$tx_s_days/365

#We have been asked to make tx_hla as continuous
#amctx$tx_hla = factor(amctx$tx_hla, ordered = T)

#############
# Subject 316 has negative rec_bmi. make it positive
###############
amctx$rec_bmi[amctx$amctx == 316] = abs(amctx$rec_bmi[amctx$amctx == 316])

##############################################
# create a data set for survival analysis. 
# Always take first measurement, 
# because rec_age otherwise can be rec_age at time of loss of follow up
# now it is rec_age at first follow up
# removing cols related to longitudinal measurements
##############################################
amctx_cumsum = cumsum(table(amctx$amctx))
first_row_index_eachsub = c(0, amctx_cumsum[-length(amctx_cumsum)]) + 1
amctx.id = amctx[first_row_index_eachsub,]
colnames(amctx.id)[which(colnames(amctx.id) == "rec_age")] = "rec_age_fwp1"

amctx$rec_age_fwp1 = rep(amctx.id$rec_age, table(amctx$amctx))

##############################################
# Create two data sets, each for pcr and creatinine
##############################################
amctx_pcr = amctx[amctx$measure=="pcr",]
amctx_creatinine = amctx[amctx$measure=="creatinine",]

#as.numeric is used otherwise for the missing person number 346, NA's are added.
amctx_pcr$visit_num = factor(unlist(sapply(table(as.numeric(amctx_pcr$amctx)), function(len){1:len})))
amctx_creatinine$visit_num = factor(unlist(sapply(table(amctx_creatinine$amctx), function(len){1:len})))

##############################################
# For certain subjects such as subject 3, 
# there are multiple measurements of same type at same time
# Checking if they exist for PCR
#############################################
idList = unique(amctx_pcr$amctx)
pcr_rep=foreach(i=1:length(idList),.combine='c') %dopar%{
  amctx_pcr_i = amctx_pcr[amctx_pcr$amctx == idList[i],]
  
  freq_time = table(amctx_pcr_i$tx_s_days)
  any(freq_time>1)
}
any(pcr_rep)

###############
# Seems no repetitions at the same time for PCR
# Creatinine table has multiple measurements per person at the same time
################
idList = unique(amctx_creatinine$amctx)
creatinine_rep=foreach(i=1:length(idList),.combine='c') %dopar%{
  amctx_creatinine_i = amctx_creatinine[amctx_creatinine$amctx == idList[i],]
  
  freq_time = table(amctx_creatinine_i$tx_s_days)
  freq_time>1
}

any(creatinine_rep)
View(amctx_creatinine[amctx_creatinine$amctx %in% idList[creatinine_rep],])

##########################################
# Take the average of the creatinine measurements which are multiple 
# at the same time for a patient
##########################################
idList = unique(amctx_creatinine$amctx)
amctx_creatinine=foreach(i=1:length(idList),.combine='rbind') %dopar%{
  amctx_creatinine_i = amctx_creatinine[amctx_creatinine$amctx == idList[i],]
  
  freq_time = table(amctx_creatinine_i$tx_s_days)
  
  retObj = amctx_creatinine_i[cumsum(freq_time) - freq_time + 1,]
  #Take the average of the observed creatinine
  timesOfInterest = as.numeric(attributes(which(freq_time>1))$names)
  for(time in timesOfInterest){
    retObj$value[which(retObj$tx_s_days == time)] = mean(amctx_creatinine_i$value[which(amctx_creatinine_i$tx_s_days == time)])
  }
  
  retObj
}

##########################################
# merge creatinine and pcr
##########################################
idList = unique(amctx$amctx)
amctx_merged=foreach(i=1:length(idList),.combine='rbind') %dopar%{
  amctx_pcr_i = amctx_pcr[amctx_pcr$amctx == idList[i],]
  amctx_creatinine_i =amctx_creatinine[amctx_creatinine$amctx == idList[i],]
  
  #create empty data frame
  amctx_i = rbind(amctx_pcr_i, amctx_creatinine_i)
  amctx_i = amctx_i[order(amctx_i$tx_s_days),]
  
  amctx_i[,c("pcr", "creatinine")] = t(sapply(1:nrow(amctx_i), function(rownum){
    #0 is purposefully used instead of NA here. We later make them NA
    if(amctx_i$measure[rownum]=="pcr"){
      c(amctx_i$value[rownum], 0)
    }else{
      c(0, amctx_i$value[rownum])
    }
  }))
  
  freq_time = table(amctx_i$tx_s_days)
  
  amctx_i_1 = amctx_i[amctx_i$tx_s_days %in% names(freq_time[freq_time==1]),]
  amctx_i_2 = amctx_i[amctx_i$tx_s_days %in% names(freq_time[freq_time==2]),]
  
  if(nrow(amctx_i_2)>0){
    for(i in 1:(nrow(amctx_i_2)/2)){
      amctx_i_2[i*2,]$pcr = amctx_i_2[i*2,]$pcr + amctx_i_2[i*2-1,]$pcr
      amctx_i_2[i*2,]$creatinine = amctx_i_2[i*2,]$creatinine + amctx_i_2[i*2-1,]$creatinine
      amctx_i_1 = rbind(amctx_i_1, amctx_i_2[i*2,])
    }
  }
  
  #removing cols which are related to value, and measurement type
  amctx_i_1[order(amctx_i_1$tx_s_days),-c(2,3,4,5)]
}

amctx_merged$pcr = sapply(amctx_merged$pcr, function(pcr){
  if(pcr==0){
    NA
  }
  else{
    pcr
  }
}, simplify = T)

amctx_merged$creatinine = sapply(amctx_merged$creatinine, function(creatinine){
  if(creatinine==0){
    NA
  }else{
    creatinine
  }
}, simplify = T)

# ###################################################
# # Another type of merging
# ###################################################
# idList = unique(amctx$amctx)
# amctx_merged=foreach(i=1:length(idList),.combine='rbind') %dopar%{
#   amctx_pcr_i = amctx_pcr[amctx_pcr$amctx == idList[i],]
#   amctx_creatinine_i =amctx_creatinine[amctx_creatinine$amctx == idList[i],]
#   amctx.id_i = amctx.id[amctx.id$amctx==idList[i], -c(2,3,4,5,6, ncol(amctx.id))]
#   
#   l_pcr = nrow(amctx_pcr_i)
#   l_creatinine = nrow(amctx_creatinine_i)
#   
#   l = max(l_pcr, l_creatinine)
#   
#   
#   creatinine = c(amctx_creatinine_i$value, rep(NA, l-l_creatinine))
#   pcr = c(amctx_pcr_i$value,  rep(NA, l-l_pcr))
#   tx_s_days_creatinine = c(amctx_creatinine_i$tx_s_days, rep(NA, l-l_creatinine))
#   tx_s_days_pcr = c(amctx_pcr_i$tx_s_days, rep(NA, l-l_pcr))
#   
#   cbind(do.call("rbind", replicate(l, amctx.id_i, simplify = FALSE)), 
#         creatinine, tx_s_days_creatinine, pcr, tx_s_days_pcr)  
# }