setwd("~/Desktop")
getwd()
af <- read.table("~/Desktop/summer_project/dataset/af.txt",header=TRUE);
af
colnames(af)[ colnames(af) %in% c('MarkerName', 'chr','pos','Effect','StdErr','P.value') ] <- c('SNP', 'CHR','BP','beta','se','P')
af_new <- af[,c(1,4,5,8,2,3,6,7)]
af_new
snp<-duplicated(af_new$SNP)
af_new[!snp,]
####generate a txt document
af_txt <- af_new[,c(1,2,3,4)]
af_txt
write.table(af_txt,file = 'af_ggplot.txt',sep = '\t',row.names = FALSE,quote=F)
af_vegas <- af_new[,c(1,4)]
af_vegas

a <- read.table("~/Desktop/snp_txt.txt",header=TRUE)
mid <- merge(af_vegas,a,by="SNP")
mid
write.table(mid,file = 'af_vegas.txt',sep = '\t',row.names = FALSE,quote=F)

####calculate total snp
af_new
newdata<-subset(af_new, af_new$P < 4.115233e-09)
newdata
dim(newdata)
####calculate total snp

cad <- read.table("~/Desktop/summer_project/dataset/cad.txt",header=TRUE);
cad
cad <- subset(cad, select=c(markername,chr,bp_hg19,effect_allele,noneffect_allele,beta,se_dgc,p_dgc))
cad
colnames(cad)[ colnames(cad) %in% c('markername', 'chr','bp_hg19','effect_allele','noneffect_allele','allele1','allele2','se_dgc','p_dgc') ] <- c('SNP', 'CHR','BP','Allele1','Allele2','se','P')
cad_new <- cad[,c(1,2,3,8,4,5,6,7)]
cad_new
snp_cad<-duplicated(cad_new$SNP)
cad_new[!snp_cad,]

####calculate total snp
cad_txt <- cad_new[,c(1,2,3,4)]
cad_txt
write.table(cad_txt,file = 'cad_ggplot.txt',sep = '\t',row.names = FALSE,quote=F)
cad_vegas <- cad_new[,c(1,4)]
cad_vegas
mid <- merge(cad_vegas,a,by="SNP")
mid
write.table(mid,file = 'cad_vegas.txt',sep = '\t',row.names = FALSE,quote=F)
cad_new
newdata<-subset(cad_new, cad_new$P < 5.287772e-09)
newdata
dim(newdata)
####calculate total snp

hf <- read.table("~/Desktop/summer_project/dataset/hf",header=TRUE);
hf
hf <- subset(hf, select=-c(EAF,INFO))
hf
colnames(hf)[ colnames(hf) %in% c('EA','NEA','BETA','SE') ] <- c('Allele1','Allele2','beta','se')
hf
hf_new <- hf[,c(1,2,3,8,4,5,6,7)]
hf_new
snp_hf<-duplicated(hf_new$SNP)
hf_new[!snp_hf,]
delete_all <- function(allele){
  result <- nchar(allele)
  new_all <- substr(allele, result, result)
}
all1 <- delete_all(hf_new$Allele1)
all1
hf_new <- data.frame(hf_new,all1)
all2 <- delete_all(hf_new$Allele2)
hf_new <- data.frame(hf_new,all2)
hf_new
hf_new <- subset(hf_new, select=-c(Allele1,Allele2))
colnames(hf_new)[ colnames(hf_new) %in% c('all1','all2') ] <- c('Allele1','Allele2')
hf_new
hf_new <- hf_new[,c(1,2,3,4,7,8,5,6)]
hf_new
####calculate total snp
hf_txt <- hf_new[,c(1,2,3,4)]
hf_txt
write.table(hf_txt,file = 'hf_ggplot.txt',sep = '\t',row.names = FALSE,quote=F)
hf_vegas <- cad_new[,c(1,4)]
hf_vegas
mid <- merge(hf_vegas,a,by="SNP")
mid
write.table(mid,file = 'hf_vegas.txt',sep = '\t',row.names = FALSE,quote=F)
hf_new
newdata<-subset(hf_new, hf_new$P < 7.153542e-09)
newdata
dim(newdata)
####calculate total snp

bmi <- read.table("~/Desktop/summer_project/dataset/bmi.txt",header=TRUE);
bmi
bmi <- subset(bmi, select=-c(Freq_Tested_Allele_in_HRS,N))
bmi
colnames(bmi)[ colnames(bmi) %in% c('POS','Tested_Allele','Other_Allele','BETA','SE') ] <- c('BP','Allele1','Allele2','beta','se')
bmi
bmi_new <- bmi[,c(3,1,2,8,4,5,6,7)]
bmi_new
####calculate total snp
bmi_txt <- bmi_new[,c(1,2,3,4)]
bmi_txt
write.table(bmi_txt,file = 'bmi_ggplot.txt',sep = '\t',row.names = FALSE,quote=F)
bmi_vegas <- bmi_new[,c(1,4)]
bmi_vegas
mid <- merge(bmi_vegas,a,by="SNP")
mid
write.table(mid,file = 'bmi_vegas.txt',sep = '\t',row.names = FALSE,quote=F)
bmi_new
newdata<-subset(bmi_new, bmi_new$P < 2.140165e-08)
newdata
dim(newdata)
####calculate total snp
af_new
cad_new
hf_new
bmi_new
data <- merge(af_new,cad_new,by="SNP",suffixes = c(".af",".cad"))
data2 <- merge(hf_new,bmi_new,by="SNP",suffixes = c(".hf ",".bmi"))
data3 <- merge(data,data2,by="SNP")

cho <- read.table("~/Desktop/summer_project/dataset/cho.txt",header=TRUE);
cho
cho <- subset(cho, select=-c(SNP_hg18,N,Freq.A1.1000G.EUR))
cho
install.packages("tidyverse")
library(tidyverse)
cho_new <- separate(data = cho, col = SNP_hg19, into = c("CHR", "BP"), sep = ":")
cho_new
get_chr <- function(CHR){
  new_str<-substr(CHR,4,5)
}
CHR <- get_chr(cho_new$CHR)
cho_new <- data.frame(cho_new,CHR)
cho_new
cho_new <- subset(cho_new, select=-c(CHR))
cho_new
to_upper <- function(a1){
  A1_upper <- toupper(a1)
}
A1_upper <- to_upper(cho_new$A1)
cho_new <- data.frame(cho_new,A1_upper)
head(cho_new)
A2_upper <- to_upper(cho_new$A2)
cho_new <- data.frame(cho_new,A2_upper)
head(cho_new)
cho_new <- subset(cho_new, select=-c(A1,A2))
cho_new <- cho_new[,c(2,6,1,5,7,8,3,4)]
colnames(cho_new)[ colnames(cho_new) %in% c('rsid','CHR.1','P.value','A1_upper','A2_upper') ] <- c('SNP','CHR','P','Allele1','Allele2')
cho_new
                      
####calculate total snp
cho_txt <- cho_new[,c(1,2,3,4)]
cho_txt     
write.table(cho_txt,file = 'cho_ggplot.txt',sep = '\t',row.names = FALSE,quote=F)

cho_vegas <- cho_new[,c(1,4)]
cho_vegas
mid <- merge(cho_vegas,a,by="SNP")
mid
write.table(mid,file = 'cho_vegas.txt',sep = '\t',row.names = FALSE,quote=F)
cho_new
newdata<-subset(cho_new, cho_new$P < 2.139512e-08)
newdata
dim(newdata)
####calculate total snp

dia <- read.table("~/Desktop/summer_project/dataset/dia.txt",header=TRUE);
dia
#####convert snp
install.packages("remotes")
remotes::install_github("RHReynolds/colochelpR")
library(remotes)
install_github("RHReynolds/colochelpR")
dia_snp <- dia$Chr.Position[1:100000]
dia_snp
m <- convert_loc_to_rs(dia_snp, dbSNP)
write.table(dia_snp,file = 'dia_snp.txt',sep = '\t',row.names = FALSE,quote=F)
#######
dia_new <- separate(data = dia, col = Chr.Position, into = c("CHR", "BP"), sep = ":")
dia_new
dia_new <- subset(dia_new, select=-c(TotalSampleSize))
dia_new
colnames(dia_new)[ colnames(dia_new) %in% c('Effect','StdErr','P.value') ] <- c('beta','se','P')
dia_new
dia_new <- dia_new[,c(1,2,7,3,4,5,6)]
dia_new
####calculate total snp
dia_txt <- dia_new[,c(1,2,3)]
dia_txt
write.table(dia_txt,file = 'dia_ggplot.txt',sep = '\t',row.names = FALSE,quote=F)
dia_new
newdata<-subset(dia_new, dia_new$P < 4.147194e-09)
newdata
dim(newdata)
####calculate total snp
data4 <- merge(cho_new,dia_new,by=c("CHR","BP"),suffixes = c(".cho ",".dia"))
data4
full_SNP <- merge(data3,data4,by = "SNP") 
full_SNP
######need to be annocated
anno_gene <- full_SNP[,c(1,2,3)]
anno_gene
colnames(anno_gene)[ colnames(anno_gene) %in% c('CHR.af', 'BP.af') ] <- c('CHR','BP')
anno_gene
dim(anno_gene)
write.table(anno_gene, "./anno_gene.txt")
######need to be annocated
S_XY_Study <- full_SNP[,c(1,5,6,7,8,14,15,21,22,28,29,35,36,40,41)]
# S_XY_Study <- subset(full_SNP, select=c(SNP,Allele1.af, Allele2.af, beta.af, se.af, beta.cad, se.cad, beta.hf, se.hf, beta.bmi, se.bmi, beta.cho, se.cho, beta.dia, se.dia))
S_XY_Study 
colnames(S_XY_Study)[ colnames(S_XY_Study) %in% c('Allele1.af','Allele2.af') ] <- c('allele_0','allele_1')
S_XY_Study
beta_fun <- function(beta,se,n){
  nor_beta <- beta/ (sqrt(n)*se)
}
nor_af <- beta_fun(S_XY_Study$beta.af,S_XY_Study$se.af,588190)
nor_af
S_XY_Study <- data.frame(S_XY_Study,nor_af)

nor_cad <- beta_fun(S_XY_Study$beta.cad,S_XY_Study$se.cad,187133)
nor_cad
S_XY_Study <- data.frame(S_XY_Study,nor_cad)


nor_hf <- beta_fun(S_XY_Study$beta.hf,S_XY_Study$se.hf,21303)
nor_hf
S_XY_Study <- data.frame(S_XY_Study,nor_hf)


nor_bmi <- beta_fun(S_XY_Study$beta.bmi,S_XY_Study$se.bmi,341698)
nor_bmi
S_XY_Study <- data.frame(S_XY_Study,nor_bmi)


nor_cho <- beta_fun(S_XY_Study$beta.cho,S_XY_Study$se.cho,188577)
nor_cho
S_XY_Study <- data.frame(S_XY_Study,nor_cho)


nor_dia <- beta_fun(S_XY_Study$beta.dia,S_XY_Study$se.dia,159208)
nor_dia
S_XY_Study <- data.frame(S_XY_Study,nor_dia)

S_XY_Study

S_XY_Study_new <- subset(S_XY_Study, select= c(SNP,allele_0,allele_1, nor_af,se.af,nor_cad,se.cad,nor_hf,se.hf.,nor_bmi,se.bmi,nor_cho,se.cho.,nor_dia,se.dia ))

S_XY_Study_new


colnames(S_XY_Study_new)[ colnames(S_XY_Study_new) %in% c('nor_af','se.af','nor_cad','se.cad','nor_hf','se.hf.','nor_bmi','se.bmi','nor_cho','se.cho.', 'nor_dia', 'se.dia') ] <- c('trait1_b','trait1_se','trait2_b','trait2_se','trait3_b','trait3_se','trait4_b','trait4_se','trait5_b','trait5_se','trait6_b','trait6_se')
S_XY_Study_new

## from here #delete the snp which is not be annnocated
a <- read.table("~/Desktop/snp_txt.txt",header=TRUE)
a <- read.table("~/Desktop/final_ld.txt",header=TRUE) ##calculate single snp to mult
a
class(a)
colnames(a)

x1 <- Reduce(intersect,list(a$SNP),S_XY_Study_new$SNP)
x1
rownames(S_XY_Study_new)<-S_XY_Study_new$SNP
S_XY_Study_new[x1,]
S_XY_full_study1<-S_XY_Study_new[x1,-1]
S_XY_full_study1
write.table(S_XY_full_study1,file = 'S_XY_full_study1.txt',sep = '\t',quote=F)
###
library(metaCCA)
data( package = 'metaCCA' )
getwd()
S_xy <- read.table("/Users/yigenannan/Desktop/summer_project/dataset/S_XY_full_study1.txt",head=TRUE); ###
S_xy
trans_factor <- function(allele){
  fac_allele <- as.factor(allele)
}
allele_0_new <- trans_factor(S_xy$allele_0)
S_xy <- data.frame(S_xy,allele_0_new)
head(S_xy)

allele_1_new <- trans_factor(S_xy$allele_1)
S_xy <- data.frame(S_xy,allele_1_new)
head(S_xy)

S_xy <- subset(S_xy, select=-c(allele_0,allele_1))
head(S_xy)

S_xy <- S_xy[,c(13,14,1,2,3,4,5,6,7,8,9,10,11,12)]
S_xy


colnames(S_xy)[ colnames(S_xy) %in% c('allele_0_new','allele_1_new') ] <- c('allele_0','allele_1')
S_xy

head(S_xy)

#######obtain sxx
sxx <- read.table("~/Desktop/S_XX_study.csv",header=TRUE,sep = ",");
sxx <- sxx[,c(3,7,9)]
new_sxx <- as.data.frame(matrix(nrow=327,ncol=327))
new_sxx
snp_txt <- unique(sxx['SNP_A'])
snp_txt
write.table(snp_txt,file = 'snp_txt.txt',sep = '\t',quote=F, row.names = FALSE)

c = unique(unlist(sxx['SNP_A']))

d = unique(unlist(sxx['SNP_B']))

e <- c(c,"rs13056283")


rownames(new_sxx)<-e
colnames(new_sxx)<-e


length(c)
write.table(e,file = 'finall_snp.txt',sep = '\t',quote=F)
c

count___=1
for (i in c){
  
  for (j in d){
    
    if (i==j){
      new_sxx[i,j] <- 1
    }
    else{
      #print(i)
      #print(j)
      temp <- sxx$R2[which(sxx$SNP_A==i & sxx$SNP_B==j)]
      #print(temp)
      
      if (length(temp)==0){
        temp <- sxx$R2[which(sxx$SNP_A==j & sxx$SNP_B==i)]
        #print(temp)
      }
      
      new_sxx[i,j] <- as.numeric(temp)
      new_sxx[j,i] <- as.numeric(temp)
    }
    

  }
  print(count___)
  count___=count___ + 1

}

class(new_sxx)
#######

dim(S_xy)
S_YY_Study_new = estimateSyy(S_XY = S_xy)
S_YY_Study_new
dim(S_YY_Study_new)
result1 <- metaCcaGp(nr_studies=1, S_XY=list(S_xy), std_info=0,
                    S_YY=list(S_YY_Study_new), N=64539)
result1
dim(result1)
rowname <- rownames(result1)
rowname
result1$SNP<-rowname
result1
class(rowname)

#####find the range of r
find_r_result<-subset(result1, result1$'-log10(p-val)' > 6.110852)
find_r_result
r_1 <- find_r_result$r_1
max(r_1)
min(r_1)
#####

#####
manhattan_sxx <- read.table("~/Desktop/manhan_ld.txt",header=TRUE);
manhattan_sxx
manhattan_table <- merge(result1,manhattan_sxx,by = "SNP") 
manhattan_table <- manhattan_table[,c(1,3,5,6)]
manhattan_table

newdata_total<-subset(manhattan_table, manhattan_table$'-log10(p-val)' > 6.110852)
newdata_total
dim(newdata_total)


write.table(manhattan_table,file = 'manhattan_table.txt',sep = '\t',row.names = FALSE,quote=F)

#####
extract <- read.table("~/Desktop/summer_project/data_process/result1.txt",header=TRUE);

extract <- extract[,c(1,2,3)]
dim(extract)
write.table(dia_txt,file = 'dia_ggplot.txt',sep = '\t',row.names = FALSE,quote=F)

######draw the total snp
# s_xx snp注释到基因上 多个snp与多个表型相关分析
snp_id <- rownames(S_xy)
snp_id
class(snp_id)
str <- list()

for (i in snp_id){
  str<-c(str,i)
}
class(str)
str
###
result2 <- metaCcaGp(nr_studies=1, S_XY=list(S_xy), std_info=0, S_YY=list(S_YY_Study_new),N=326, analysis_type=2, SNP_id=c("rs17804371","rs4787006"),S_XX = list(new_sxx))
result2[1:2]
result2[3]
result2[4]
S_XX=list(new_sxx)
