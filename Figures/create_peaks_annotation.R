peaks_reads<-read.table("~/Desktop/R所需文件/idrpeaks_depth_rbq.txt")
peaks_reads[,4]<-peaks_reads[,4]/22706283*1000000
peaks_reads[,5]<-peaks_reads[,5]/6987114*1000000
peaks_reads[,6]<-peaks_reads[,6]/26255468*1000000
peaks_reads[,7]<-peaks_reads[,7]/20195449*1000000
colnames(peaks_reads)<-c("chr","start","end","G1","ES","MS","LS")
#peaks_reads<-peaks_reads[which( !(peaks_reads[,4]>peaks_reads[,5] & peaks_reads[,4]>peaks_reads[,6] & peaks_reads[,4]>peaks_reads[,7])),]
peaks_reads[,8]<-apply(peaks_reads[,5:7],1,max)
peaks_reads[,9:11]<-0
peaks_reads$RT<-0
for(i in 1:nrow(peaks_reads))
{
  if(peaks_reads[i,5] > 0.9*peaks_reads[i,8])
  {peaks_reads[i,9]="E"}
  if(peaks_reads[i,6] > 0.9*peaks_reads[i,8])
  {peaks_reads[i,10]="M"}
  if(peaks_reads[i,7] > 0.9*peaks_reads[i,8])
  {peaks_reads[i,11]="L"}
  peaks_reads[i,12]<-gsub("0","",paste(peaks_reads[i,9],peaks_reads[i,10],peaks_reads[i,11],sep = ""))
}
peaks_RT<-peaks_reads[,c(1,2,3,12)]
peaks_RT[,5]<-0
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
# write.table(peaks_RT,"~/Desktop/rbq_gff.gff3",row.names = F,quote=F,sep = "\t",col.names = F)




#####ZH11
peaks_reads_ZH11<-read.table("~/Desktop/R所需文件/idrpeaks_depth_ZH11.txt")
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