##Sepsis: survivor vs non-survivor
library(Biobase)
library(GEOquery)
library(limma)
library(knitr)
library(ggplot2)
library(maSigPro)
library(COCONUT)
library(hgu133plus2.db)
source('./helpers.r')
setwd("your_path/") 

Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 2)
gsetGSE48080 <- getGEO("GSE48080", destdir = ".",GSEMatrix =TRUE, getGPL=FALSE)
if (length(gsetGSE48080) > 1) idx <- grep("GPL4133", attr(gsetGSE48080, "names")) else idx <- 1
gsetGSE48080<-gsetGSE48080[[idx]]

gsetGSE54514 <- getGEO("GSE54514", destdir = ".",GSEMatrix =TRUE, getGPL=FALSE)
if (length(gsetGSE54514) > 1) idx <- grep("GPL6947", attr(gsetGSE54514, "names")) else idx <- 1
gsetGSE54514<-gsetGSE54514[[idx]]

##Getting GPL data
gplGSE48080<-getGEO("GPL4133",destdir = ".")
gpl_table<-gplGSE48080@dataTable@table

gplGSE54514<-getGEO("GPL6947",destdir = ".")
gpl_table1<-gplGSE54514@dataTable@table


gsetGSE95233 <- getGEO("GSE95233", destdir = "your_path/", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gsetGSE95233) > 1) idx <- grep("GPL570", attr(gsetGSE95233, "names")) else idx <- 1
gsetGSE95233 <- gsetGSE95233[[idx]]

gsetGSE95233 <- probeGeneconv(gsetGSE95233,hgu133plus2SYMBOL,'SYMBOL')

## Making the number of genes equal in both datasets for COCONUT
ex<-exprs(gsetGSE48080)
gpl_gene_symbolGSE48080<-gpl_table[[10]]
gpl_gene_symbolGSE48080<-gpl_gene_symbolGSE48080[as.numeric(rownames(ex))]
rownames(ex)<-gpl_gene_symbolGSE48080
temp<-ex[order(rownames(ex)),]
temp<-temp[12320:45015,]
temp<-temp[!duplicated(rownames(temp)),]
gpl_gene_symbolGSE48080<-rownames(temp)


ex1<-exprs(gsetGSE54514)
gpl_gene_symbolGSE54514<-gpl_table1[[14]]
gpl_gene_symbolGSE54514<-gpl_gene_symbolGSE54514[gpl_table1[[1]] %in% rownames(ex1)]
rownames(ex1)<-gpl_gene_symbolGSE54514
temp1<-ex1[order(rownames(ex1)),]
temp1<-temp1[3758:24840,]
temp1<-temp1[!duplicated(rownames(temp1)),]
gpl_gene_symbolGSE54514<-rownames(temp1)

ex2<-exprs(gsetGSE95233)
temp2<-ex2[order(rownames(ex2)),]
temp2<-temp2[!duplicated(rownames(temp2)),]
gpl_gene_symbolGSE95233<-rownames(temp2)


gpl_gene_symbol_GSE54514GSE48080 <- intersect(gpl_gene_symbolGSE54514,gpl_gene_symbolGSE48080)
gpl_gene_symbol_all <- intersect(gpl_gene_symbol00,gpl_gene_symbolGSE95233)

ex_GSE48080<-temp[gpl_gene_symbol_all,]
ex_GSE54514<-temp1[gpl_gene_symbol_all,]
ex_GSE95233<-temp2[gpl_gene_symbol_all,]

ex_GSE48080 <- 2^ex_GSE48080
ex_GSE54514 <- 2^ex_GSE54514
ex_GSE95233 <- 2^ex_GSE95233

##Separating Healthy controls and patients in both datasets
disease_status<-gsetGSE48080$`outcome:ch1`
temp<-as.numeric(disease_status=="Healthy volunteer")+1
disease_status<-temp + as.numeric(disease_status=="Dead")*(-1)
ctrl<-ex_GSE48080[,disease_status==2]
ex_GSE48080<-ex_GSE48080[,disease_status!=2]

disease_status<-gsetGSE54514$`disease status:ch1`
temp<-as.numeric(disease_status=="healthy")+1
disease_status<-temp + as.numeric(disease_status=="sepsis nonsurvivor")*(-1)
ctrl1<-ex_GSE54514[,disease_status==2]
ex_GSE54514<-ex_GSE54514[,disease_status!=2]

disease_status<-gsetGSE95233$survival.ch1
temp<-as.numeric(disease_status=="NA")+1
disease_status<-temp + as.numeric(disease_status=="Non Survivor")*(-1)
ctrl2<-ex_GSE95233[,disease_status==2]
ex_GSE95233<-ex_GSE95233[,disease_status!=2]
dim(ex_GSE95233)

##Preparing both datasets for COCONUT
GSE48080<-list()
Healthy0.Sepsis1<-c(0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
col_names<-colnames(cbind(log2(ctrl),log2(ex_GSE48080)))
GSE48080$pheno<-matrix(c(Healthy0.Sepsis1,Healthy0.Sepsis1),nrow = 23, ncol=2, dimnames=list(col_names,c("Healthy0.Sepsis1","duplicate")))
GSE48080$genes<-data.frame(cbind(log2(ctrl),log2(ex_GSE48080)))

GSE54514<-list()
col_names<-colnames(cbind(ctrl1,ex_GSE54514))
Healthy0.Sepsis1<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,
                    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                    1,1,1)
GSE54514$pheno<-matrix(c(Healthy0.Sepsis1,Healthy0.Sepsis1),nrow = 163, ncol = 2,dimnames = list(col_names,c("Healthy0.Sepsis1","duplicate")))
GSE54514$genes<-data.frame(cbind(ctrl1,ex_GSE54514))

GSE95233<-list()
Healthy0.Sepsis1<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                    0,0,1,1,
                    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
col_names<-colnames(cbind(ctrl2,ex_GSE95233))
GSE95233$pheno<-matrix(c(Healthy0.Sepsis1,Healthy0.Sepsis1),nrow = 124, ncol=2, dimnames=list(col_names,c("Healthy0.Sepsis1","duplicate")))
GSE95233$genes<-data.frame(cbind(ctrl2,ex_GSE95233))


GSEs.data<-list()
GSEs.data$GSE48080<-GSE48080
GSEs.data$GSE54514<-GSE54514
GSEs.data$GSE95233<-GSE95233

result<-COCONUT(GSEs.data, control.0.col = "Healthy0.Sepsis1")

norm_exGSE48080<-result$COCONUTList$GSE48080$genes
norm_exGSE54514<-result$COCONUTList$GSE54514$genes
norm_exGSE95233<-result$COCONUTList$GSE95233$genes

norm_exGSE48080 <- log2(norm_exGSE48080)
norm_exGSE54514 <- log2(norm_exGSE54514)
norm_exGSE95233 <- log2(norm_exGSE95233)



write.csv(norm_exGSE48080,file = "your_path/coconut_GSE48080.csv")
write.csv(norm_exGSE54514,file = "your_path/coconut_GSE54514.csv")
write.csv(norm_exGSE95233,file = "your_path/coconut_GSE95233.csv")


