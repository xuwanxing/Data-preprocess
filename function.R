library(dplyr)
library(DESeq2)
library(compositions)
#RNAseq data: variance stabilizing transformation
#reference:https://www.jianshu.com/p/9864be30ebce

data1<-read.table("./raw/pro_A.tsv",header = T)
data1_deal<-varianceStabilizingTransformation(data.matrix(round(data1)))
data1_trans<-t(data1_deal)
write.table(data1_trans,"./normalizaion/expression_pro_A.txt",col.names = T, row.names = T,sep = "\t")
#should revise in excel


#abundance: centered log ratio
#reference:http://events.jianshu.io/p/0763a167f98a
data2<-read.csv("./modified_taxon_A.tsv",header = T,sep = "\t")
data2$COM<-paste(data2$D,data2$P,data2$C,data2$O,data2$F,data2$G,data2$S,sep="_")
data2_deal<-data2[,-c(1,2,3,4,5,6,7)]
rownames(data2_deal)<-data2_deal[,"COM"]
data2_deal<-subset(data2_deal,select = -c(COM))
#write.table(data2_deal,"test.txt",col.names = T, row.names = T, sep = "\t")
colSums(data2_deal)
#filter
precent1 <- as.data.frame(apply(data2_deal, 2, function(x){x/sum(x)}))
#precent2 <-  as.data.frame(t(t(data2_deal)/colSums(data2_deal,na=T))*100)
colSums(precent1)
#colSums(precent2)

cutoff = 0.1
filter <- data.frame(precent1[which(apply(precent1, 1, function(x){
  length(which(x>=0.001))/length(x)}) > cutoff),])
                                                                                               

#clr
#data2_clr<-clr(filter)
filter[filter==0]<-1e-10
data2_clr<-clr(filter)


final<-t(as.data.frame(data2_clr))
final_frame<-as.data.frame(final)

#Standard deviation
#my_sd<-vector()
#for (i in 1:ncol(final_frame)) {
#  my_sd=c(my_sd,sd(final_frame[,i]))
#  
#}
#any(my_sd=0)

write.table(final_frame,"./normalizaion/abundance_taxon_A_test.txt",col.names = T, row.names = T,sep = "\t")
#should revise in excel
