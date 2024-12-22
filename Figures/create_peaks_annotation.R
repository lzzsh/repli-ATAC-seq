library(MASS)
library(moments)
library(ggplot2)

# normalization the reads of each RT
peaks_reads<-read.table("~/Desktop/Rfiles/idr_peaks/EdU-idr_reads_ZH11.txt", sep = "\t")
colnames(peaks_reads)<-c("chr","start","end","ES","MS","LS")
peaks_reads$ES<-peaks_reads$ES/2482331*1000000 #3987051
peaks_reads$MS<-peaks_reads$MS/42792694*1000000 #10164368
peaks_reads$LS<-peaks_reads$LS/15965054*1000000 #10328346

#peaks_reads<-peaks_reads[which( !(peaks_reads[,4]>peaks_reads$ES & peaks_reads[,4]>peaks_reads$MS & peaks_reads[,4]>peaks_reads$LS)),]
peaks_reads$max<-apply(peaks_reads[,4:6],1,max)
peaks_reads[,8:10]<-0
peaks_reads$RT<-0

# normalization and then remove outliers
if (min(peaks_reads$max, na.rm = TRUE) <= 0) {
  adjustment <- 1 - min(peaks_reads$max, na.rm = TRUE)
  peaks_reads$max_adjusted <- peaks_reads$max + adjustment
} else {
  peaks_reads$max_adjusted <- peaks_reads$max
}

# Perform Box-Cox transformation
bc_lambda <- boxcox(peaks_reads$max_adjusted ~ 1, data = peaks_reads, lambda = seq(-2, 2, len = 50))
lambda_opt <- bc_lambda$x[which.max(bc_lambda$y)]

# Apply the optimal λ value for transformation
data_bc <- if (lambda_opt == 0) {
  log(peaks_reads$max_adjusted)
} else {
  ((peaks_reads$max_adjusted^lambda_opt - 1) / lambda_opt)
}

# Calculate skewness and output the result
skew_bc <- skewness(data_bc)
cat("Skewness after Box-Cox Transformation:", skew_bc, "\n")

# Calculate mean and standard deviation
mean_value <- mean(data_bc, na.rm = TRUE)
std_dev <- sd(data_bc, na.rm = TRUE)

# Define the bounds for outliers
lower_bound <- mean_value - 2 * std_dev
upper_bound <- mean_value + 2 * std_dev

# Add a column to label outliers and normal values
peaks_reads$is_outlier <- ifelse(
  data_bc < lower_bound | data_bc > upper_bound,
  "Outlier",
  "Normal"
)

# Count the number of outliers and normals
outliers_count <- sum(peaks_reads$is_outlier == "Outlier", na.rm = TRUE)
normals_count <- sum(peaks_reads$is_outlier == "Normal", na.rm = TRUE)

cat("Number of outliers:", outliers_count, "\n")
cat("Number of normals:", normals_count, "\n")

# Plot the density plot of original data with outlier detection
ggplot(data.frame(x = peaks_reads$max_adjusted, is_outlier = peaks_reads$is_outlier), aes(x = x, fill = is_outlier)) +
  geom_density(alpha = 0.6) +
  theme_minimal() +
  labs(title = "Density Plot of Data with Outlier Detection", x = "Value", y = "Density") +
  scale_fill_manual(values = c("Outlier" = "#FF9999", "Normal" = "#99CCFF"))

# Plot the density plot of data after Box-Cox transformation with outlier detection
ggplot(data.frame(x = data_bc, is_outlier = peaks_reads$is_outlier), aes(x = x, fill = is_outlier)) +
  geom_density(alpha = 0.6) +
  theme_minimal() +
  labs(title = "Density Plot after Box-Cox Transformation with Outlier Detection", x = "Transformed Value", y = "Density") +
  scale_fill_manual(values = c("Outlier" = "#FF9999", "Normal" = "#99CCFF"))

for(i in 1:nrow(peaks_reads))
{
  if(peaks_reads[i,4] > 0.9*peaks_reads[i,7])
  {peaks_reads[i,8]="E"}
  if(peaks_reads[i,5] > 0.9*peaks_reads[i,7])
  {peaks_reads[i,9]="M"}
  if(peaks_reads[i,6] > 0.9*peaks_reads[i,7])
  {peaks_reads[i,10]="L"}
  peaks_reads[i,11]<-gsub("0","",paste(peaks_reads[i,8],peaks_reads[i,9],peaks_reads[i,10],sep = ""))
}

peaks_RT<-peaks_reads[,c(1,2,3,11)]
peaks_RT$ES<-0
for(i in 1:nrow(peaks_RT))
{
  if(peaks_RT[i,2]==0)
  {peaks_RT[i,2]=peaks_RT[i,2]+1}
  if(peaks_RT[i,4]=="E")
  {peaks_RT[i,5]="Name=ES;color=#2250F1;"}
  if(peaks_RT[i,4]=="EM")
  {peaks_RT[i,5]="Name=ESMS;color=#28C5CC;"}
  if(peaks_RT[i,4]=="M")
  {peaks_RT[i,5]="Name=MS;color=#1A8A12;"}
  if(peaks_RT[i,4]=="ML")
  {peaks_RT[i,5]="Name=MSLS;color=#FFFD33;"}
  if(peaks_RT[i,4]=="L")
  {peaks_RT[i,5]="Name=LS;color=#FB0018;"}
  if(peaks_RT[i,4]=="EL")
  {peaks_RT[i,5]="Name=ESLS;color=#EA3CF2;"}
  if(peaks_RT[i,4]=="EML")
  {peaks_RT[i,5]="Name=ESMSLS;color=#FAB427;"}
}
peaks_RT[,6]<-"."
peaks_RT[,7]<-"."
peaks_RT[,8]<-"."
peaks_RT[,9]<-"."
peaks_RT[,10]<-"peaks"
peaks_RT<-peaks_RT[,c(1,6,10,2,3,7,8,9,5)]
# write.table(peaks_reads,"~/Desktop/Rfiles/idr_peaks/peaks_reads.txt",row.names = F,quote=F,sep = "\t",col.names = T)
# write.table(peaks_RT,"~/Desktop/Rfiles/idr_peaks/rbq_gff.gff3",row.names = F,quote=F,sep = "\t",col.names = F)




#####ZH11
peaks_reads_ZH11<-read.table("~/Desktop/Rfiles/idrpeaks_depth_ZH11.txt")
peaks_reads_ZH11[,4]<-peaks_reads_ZH11[,4]/43880448*1000000
peaks_reads_ZH11[,5]<-peaks_reads_ZH11[,5]/2474361*1000000
peaks_reads_ZH11[,6]<-peaks_reads_ZH11[,6]/42767849*1000000
peaks_reads_ZH11[,7]<-peaks_reads_ZH11[,7]/15943293*1000000
colnames(peaks_reads_ZH11)<-c("chr","start","end","G1","ES","MS","LS")
#peaks_reads_ZH11<-peaks_reads_ZH11[which( !(peaks_reads_ZH11[,4]>peaks_reads_ZH11[,5] & peaks_reads_ZH11[,4]>peaks_reads_ZH11[,6] & peaks_reads_ZH11[,4]>peaks_reads_ZH11[,7])),]
peaks_reads_ZH11[,8]<-apply(peaks_reads_ZH11[,5:7],1,max)
peaks_reads_ZH11[,9:11]<-0
peaks_reads_ZH11$RT<-0
for(i in 1:nrow(peaks_reads_ZH11))
{
  if(peaks_reads_ZH11[i,5] > 0.9*peaks_reads_ZH11[i,8])
  {peaks_reads_ZH11[i,9]="E"}
  if(peaks_reads_ZH11[i,6] > 0.9*peaks_reads_ZH11[i,8])
  {peaks_reads_ZH11[i,10]="M"}
  if(peaks_reads_ZH11[i,7] > 0.9*peaks_reads_ZH11[i,8])
  {peaks_reads_ZH11[i,11]="L"}
  peaks_reads_ZH11[i,12]<-gsub("0","",paste(peaks_reads_ZH11[i,9],peaks_reads_ZH11[i,10],peaks_reads_ZH11[i,11],sep = ""))
}
peaks_RT_ZH11<-peaks_reads_ZH11[,c(1,2,3,12)]
peaks_RT_ZH11[,5]<-0
for(i in 1:nrow(peaks_RT_ZH11))
{
  if(peaks_RT_ZH11[i,2]==0)
  {peaks_RT_ZH11[i,2]=peaks_RT_ZH11[i,2]+1}
  if(peaks_RT_ZH11[i,4]=="E")
  {peaks_RT_ZH11[i,5]="Name=ES;color=#2250F1;"}
  if(peaks_RT_ZH11[i,4]=="EM")
  {peaks_RT_ZH11[i,5]="Name=ESMS;color=#28C5CC;"}
  if(peaks_RT_ZH11[i,4]=="M")
  {peaks_RT_ZH11[i,5]="Name=MS;color=#1A8A12;"}
  if(peaks_RT_ZH11[i,4]=="ML")
  {peaks_RT_ZH11[i,5]="Name=MSLS;color=#FFFD33;"}
  if(peaks_RT_ZH11[i,4]=="L")
  {peaks_RT_ZH11[i,5]="Name=LS;color=#FB0018;"}
  if(peaks_RT_ZH11[i,4]=="EL")
  {peaks_RT_ZH11[i,5]="Name=ESLS;color=#EA3CF2;"}
  if(peaks_RT_ZH11[i,4]=="EML")
  {peaks_RT_ZH11[i,5]="Name=ESMSLS;color=#FAB427;"}
}
peaks_RT_ZH11[,6]<-"."
peaks_RT_ZH11[,7]<-"."
peaks_RT_ZH11[,8]<-"."
peaks_RT_ZH11[,9]<-"."
peaks_RT_ZH11[,10]<-"peaks"
peaks_RT_ZH11<-peaks_RT_ZH11[,c(1,6,10,2,3,7,8,9,5)]
# write.table(peaks_RT_ZH11,"~/Desktop/ZH11_gff.gff3",row.names = F,quote=F,sep = "\t",col.names = F)