library(Biobase)
library(GEOquery)
library(limma)
library(knitr)
library(ggplot2)
library(maSigPro)
library(EnrichmentBrowser)
library(illuminaHumanv3.db) 


#Process the data
gset <- getGEO("GSE54514", destdir = "/Users/mhuan98/Desktop/emory_research/rscripts_test", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL6947", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
gset <- process_GSE54514_labels(gset) 
gset <- probeGeneconv(gset,illuminaHumanv3SYMBOL,'SYMBOL') # go from probe id's to entrex geneID

gset <- gset[,gset$group != "HC"]
gset <- gset[,gset$time != 4]
gset <- gset[,gset$time != 5]


time_interest1 <- c(1,2,3)
for(e in gset$patient_id){
  sub <- gset[,gset$patient_id == e]
  if(  all(time_interest1 %in% sub$time) ==FALSE){
    print(sprintf("%s doesn't have all time point of interest",e) )
    gset<- gset[,gset$patient_id != e]
  }
}

expression_df <- as.data.frame(exprs(gset))

gset_df <- as.data.frame(gset)
gset_df$patient_id <- unlist(gset$patient_id)
gset_df$survivor <- as.integer(gset$group=="S")
gset_df$noSurvivor <- as.integer(gset$group =="NS")

gset_df$time <- gset$time
survival  <- 0
dead <- 0 
for(id in unique(gset_df$patient_id)){
  patient_subset <- gset_df[gset_df$patient_id == id,]
  already =0
  if(all(patient_subset$survivor == 1)){
    already =1
    survival <- survival +1 
  }
  if (all(patient_subset$noSurvivor == 1)){
    if(already ==1){
      print("ERROR")
    }
    dead <- dead + 1
  }
}

#make the table for Masigpro
items_interest <- gset_df[c("patient_id","time","survivor","noSurvivor")]
final <- cbind(expression_df[row.names(items_interest),],items_interest)
tems_interest['other'] = 0  
items_interest[ items_interest$time==1 & items_interest$survivor== 1,'other'] = 1 
items_interest[ items_interest$time ==1 & items_interest$noSurvivor==1,'other'] = 2 
items_interest[items_interest$time ==2 & items_interest$survivor ==1,'other'] = 3 
items_interest[items_interest$time==2 & items_interest$noSurvivor==1,'other'] =4 
items_interest[items_interest$time ==3 & items_interest$survivor ==1,'other'] =5 
items_interest[items_interest$time==3 & items_interest$noSurvivor==1,'other'] =6 
items_interest <- items_interest[c('time','other','survivor','noSurvivor')]


design_matrix <- as.matrix.data.frame(items_interest )

design <- make.design.matrix(design_matrix,degree=1)

#read the coconut normalization file
norm_ex <- read.csv(file = 'your_path/cocomut_GSE54514.csv',row.names = 1)

#steps in Masigpro
fit <- p.vector(norm_ex, design, Q = 0.05, MT.adjust = "fdr")

tstep <- T.fit(
  fit,
  step.method = "backward",
  alfa = 0.05)

sigs <- get.siggenes(
  tstep,
  rsq = 0.3,
  vars = "groups")

genes_masig <- see.genes(
  sigs$sig.genes$noSurvivorvssurvivor,
  show.fit = T,
  dis =design$dis,
  cluster.method="hclust",
  cluster.data = 1,
  k = 4,
  newX11 = F)

write.csv(genes_masig$cut, 'your_path/GSE54514_masig.csv')















