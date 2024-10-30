setwd("/home/fengchr/staff/2024_4_24_zhk_oneseq/calitas")
data<-read.csv("ttnAGGACAGAGGGTCAGCATGCCA.calitas_searchreference.out.tsv",sep = "\t",header = T,stringsAsFactors = F)
data<-data[data$total_mm_plus_gaps <= 6 & data$guide_gaps <= 2,]
write.table(data,"ttnAGGACAGAGGGTCAGCATGCCA.calitas_searchreference.out.filtered.tsv",sep = "\t",row.names = F)
