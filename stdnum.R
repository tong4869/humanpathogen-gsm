

len<-read.table("seqstat.txt",sep='\t')
colnames(len)<-c("samp","len","rnum")
library(dplyr)
len1<-mutate(len,seats=(len-50)*rnum/(10^10))
scol<-strsplit(len1$samp,"_")
scol1<-data.frame(do.call(rbind,scol))
len2<-cbind(len1,scol1)
len3<-aggregate(x=len2$seats,by=list(len2$X1),sum)
colnames(len3)<-c("smpname","stdnum")
write.table(len3,"stdnum.txt",sep = '\t')

