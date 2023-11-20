##picrust2数据预处理
rm(list=ls())
####Input data
setwd("D:/D1/2A-SJTU/1E-总结/2023/manuscript/submit材料/mDriver/real_data_new/UC/UC_picrust2")
count <- read.delim("count.txt", row.names = 1, sep = '\t', check.names = FALSE)
metadata <- read.delim("metadata.txt", row.names = 1, sep = '\t', check.names = FALSE)
silva_species <- read.delim("silva_species.txt", row.names = 1, sep = '\t', check.names = FALSE)
count <-  t(count)

count_meta <- merge(metadata,count,by = 'row.names', all = F)
feature <- count_meta[, -c(2:4)]
row.names(feature) <- feature[,1]
feature <- feature[,-1]

col_sums <- colSums(feature)
zero_microbes <- which(col_sums == 0)
filtered_feature <- feature[, -zero_microbes]

#feature_table
### write rownames of data
adjustdata <- function(data) {
  data<-cbind(rownames(data),data)
}
otu_table <- adjustdata(t(filtered_feature))
write.table(otu_table,file ="otu_table.txt",row.names = F,col.names = T, sep = "\t",quote = F)

#otu.fasta
filter_otu <- silva_species[,c(1:2)]
tfiltered_feature <- t(filtered_feature)
merge_filter_otu <- merge(filter_otu,tfiltered_feature,by = 'row.names', all = F)
fasta <- merge_filter_otu[,c(1:2)]
fasta <- paste0(">", fasta$Row.names, "\n", fasta$sequence)
writeLines(fasta, "otu.fasta")


