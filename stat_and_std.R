#/usr/bin/Rscript std_1229.R --rdir /home/litong/2023dec/workflow/results --thr 2 --speid /home/litong/2023dec/workflow/speid.txt --outdir /home/litong/2023dec/workflow/out --std /home/litong/2023dec/workflow/stdnum.txt

library(getopt)
library(dplyr)
library(tidyr)
library(reshape2)

command=matrix(c( 
  'rdir', 'r', 1,'character',
  'thr', 't', 1, 'integer',
  'speid','i',1,'character',
  'outdir','d',1,'character',
  'std','s','1','character'
),byrow=TRUE, ncol=4)

options=getopt(command)

setwd(options$outdir)



count.freq <- function(outfile,acc){
  stat <- read.table(outfile,header=F,sep='\t')
  stat <- stat[,1:2]
  stat <- as.data.frame(table(stat$V1))
  rownames(stat) <- stat[,1]
  colnames(stat) <- c("v1",acc)
  stat %>% dplyr::select(acc)
}


df_list<-list()
for (i in 1:10){
pat<-paste0("*set", i, ".out")

filelist <- dir(path = options$rdir ,full.names = T,pattern = pat)

sample_gsm_stat <- data.frame()
sample_patho_stat_species <-  data.frame()
sample_patho_stat <-  data.frame()


emptyfile <- c()
for (file in filelist){    
  acc_r <- sub(".*/([^/]+)_set\\d+\\.out", "\\1", file)
  acc <- sub(".*/([^/]+)_\\d+_set\\d+\\.out", "\\1", file)
  if (file.info(file)$size != 0) {
    tmp <- as.data.frame(t(count.freq(file,acc_r)))
    sample_gsm_stat <- bind_rows(sample_gsm_stat,tmp)
  }else{
    emptyfile <- c(emptyfile,acc_r)
  }
}


sample_gsm_stat [is.na(sample_gsm_stat )] <- 0 


accrow<-c()
for (j in 1:nrow(sample_gsm_stat)){
  acc<-strsplit(rownames(sample_gsm_stat),'_')[[j]][1]
  accrow<-c(accrow,acc)
}


stat_species <- cbind(accrow,sample_gsm_stat)
stat_gsm_bind <- aggregate(stat_species[-1],stat_species["accrow"],sum)
rownames(stat_gsm_bind)<-stat_gsm_bind$accrow
stat_gsm_bind <- stat_gsm_bind[-1]

stat_gsm_bind <- as.data.frame(t(stat_gsm_bind))

gsmcol <- rownames(stat_gsm_bind)
gsmcol1<-gsmcol
stat_gsm_bind <- cbind(gsmcol,gsmcol1,stat_gsm_bind)
stat_gsm_bind_spe <- separate(stat_gsm_bind,gsmcol1,c("species",NA),"_")


sample_spe_thr<-data.frame()
thr <- options$thr


for (n in 3:ncol(stat_gsm_bind_spe)) {
  
  acc <- colnames(stat_gsm_bind_spe)[n]
  tmpcol<-data.frame()
  
  for (spe in unique(stat_gsm_bind_spe$species)) {
    
    tmpspe <- stat_gsm_bind_spe %>% dplyr::select("gsmcol","species",acc) %>% filter(species==spe)
    
    if ( sum(tmpspe[,acc]) >= 1 ) {
      
      detected <- filter(tmpspe,.data[[acc]]!=0)
      n_gsm <- length(unique(detected$gsmcol))
      
      if ( n_gsm >= thr){
        tmpline <- c(spe,sum(tmpspe[,acc]))
        tmpcol <- rbind(tmpcol,tmpline)
      }
    }
  }
  
  if (nrow(tmpcol)!=0){
    rownames(tmpcol) <- tmpcol[,1]
    colnames(tmpcol) <- c("delete",acc)
    tmpcol<-dplyr::select(tmpcol,acc)
    tmprow<-as.data.frame(t(tmpcol))
    sample_spe_thr <- bind_rows(sample_spe_thr,tmprow)
  }
}

df <- data.frame(matrix(0, ncol = ncol(sample_spe_thr), nrow = length(emptyfile)))
rownames(df) <- emptyfile
colnames(df) <- colnames(sample_spe_thr)


sample_spe_thr <- rbind(df,sample_spe_thr)


sample_spe_thr [is.na(sample_spe_thr )] <- 0 
sample_spe_thr_t  <- as.data.frame(t(sample_spe_thr ))
species<-rownames(sample_spe_thr_t)
sample_spe_thr_t <- as.data.frame(lapply(sample_spe_thr_t, as.numeric))
rownames(sample_spe_thr_t)<-species


std <- read.table(options$std,sep='\t',header=T)

rownames(std)<-std$smpname

std.table <- matrix()
for(k in 1:ncol(sample_spe_thr_t)){
  acc <- colnames(sample_spe_thr_t)[k]
  stdnum <- std[acc,][1,2]
  tmp <- dplyr::select(sample_spe_thr_t,acc)/stdnum
  std.table <- cbind(std.table,tmp)
}

std.table.1 <- as.data.frame(std.table[,-1, drop = FALSE])

std.table.1$speid<-rownames(std.table.1)

speid<-read.table(options$speid,sep='\t',header = T)

std.table.2<-merge(std.table.1,speid,by="speid")
std.table.2 <- subset(std.table.2, select = -speid)

df_list[[i]] <- std.table.2

}

combine <- bind_rows(df_list, .id = "Iteration")

custom <- function(x) {
  sum(x) / 10
}
result<-aggregate(combine[2:(ncol(combine)-1)],by=list(spename=combine$spename),FUN=custom)

write.table(result,"PathogensStat.txt",sep = '\t',col.names = T,row.names = F)


