#set working directory
setwd("/home/....")

# Creating a sample dataframe
df <- data.frame(
  Gene=c("A_2kb","B_4kb","C_1kb","D_10kb"),
  Rep1=as.numeric(c(15,30,5,0)),
  Rep2=as.numeric(c(14,20,6,0)),
  Rep3=as.numeric(c(30,40,15,1))
)

#get the length of respective genes in kb
length<-c(2,4,1,10)

#read counts from data
counts<-df[,-1]

# Calculate the total number of reads in each sample
total_reads <- colSums(counts)

# Calculate the scaling factor for each sample

##########################################################
#############---------IMPORTANT-------####################
#############--CHANGE 10 to 10e6 below while dealing with real datasets--(here 10 is used as data is less)
##########################################################
scaling_factors <- total_reads/10  #CHANGE 10 TO 10e6

# Convert  to RPM
rpm<-counts
for (i in 1:3) {
  rpm[, i] <- counts[, i]/scaling_factors[i]
}

# Convert  to RPKM
rpkm<-rpm
for (i in 1:4) {
  rpkm[i,] <- rpkm[i,]/length[i]
}

#########################################################
#TO SAVE FILES
#########################################################

#change column names of rpm values
colnames(rpm)<-c("Rep1_rpm","Rep2_rpm","Rep3_rpm")

#merge rpm dataframe with genes
rpm_df<-cbind(df,rpm)

#save the dataframe
write.csv(rpm_df,"rpm.csv",row.names = F)

#change column names of rpkm values
colnames(rpkm)<-c("Rep1_rpkm","Rep2_rpkm","Rep3_rpkm")

#merge rpm dataframe with genes
rpkm_df<-cbind(df,rpkm)

#save the dataframe
write.csv(rpkm_df,"rpkm.csv",row.names = F)
