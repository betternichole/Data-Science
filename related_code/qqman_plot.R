install.packages("qqman")
library(qqman)
getwd()
pdf(file = '~/Desktop/summer_project/data_process/manhattan_total.pdf',width = 15)
af_rt = read.table("~/Desktop/summer_project/data_process/total_ggplot.txt",header = T)
af_rt
qqq<-0.05/64539
qqq
-log10(qqq)
manhattan(af_rt,chr = "CHR",bp = "BP",p = "X.log10.p.val.",snp = "SNP",col = c("blue4","orange3"),ylim = c(0,10),suggestiveline = F,genomewideline = -log10(qqq))
dev.off()

##dim() dia hf