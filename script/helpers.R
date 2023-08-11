library(Biobase)
library(GEOquery)
library(limma)
library(knitr)
library(ggplot2)
library(maSigPro)
library(EnrichmentBrowser)

getNullConvs <- function(items){
  nullid <- c()
  i=1 
  counter =1 
  for( e in items){
    if( is.null(e) ){
      nullid[counter] = i 
      counter = counter +1 
    }
    i  = i + 1
  }
  nullid
}

probeGeneconv <- function(gset,id2probe_ob,ID_NAME){
  #gset <- is the microarray ExperimetnSummary provided by geo 
  #id2probe_obj is the mapping object provided by the illumnia.db library to map probes onto some id or symbol  
  #ID_NAME <- string object speciying the field name of the conversion to be applied. if converting to SYMBOL then ID_NAME=SYMBOL
  probe_names <- rownames(gset)
  #first let's try to remove probes that do not have a map to  ENTREZ IDS 
  mapped_probes <- mappedkeys(id2probe_ob)# Convert to a lis
  probe_gene_map <- as.list(id2probe_ob[mapped_probes])
  ID_obj <- lapply(rownames(gset) , function(x){unlist(probe_gene_map[x][[1]])}) #this gets the numeric entry of the names 
  # probes that do not map to regular gene IDS  are mapped to null. these will be elimated so i can use the probe2gene function. 
  nullid <- getNullConvs(ID_obj)
  null_names <- probe_names[nullid]
  #get the none null names
  ID_obj <-Filter( Negate(is.null) , ID_obj )
  ID_obj <- unlist(ID_obj)  # this is a strange property of R where the outputs of the lappply fucntio nare 
  # nested.  This pulls them out. 
  new_gset <- gset[!(row.names(gset) %in% null_names), ]  # make a new gset obj where only mappable probes are present
  #probe2gene functions requires object be converted to SummarizedExperiment object
  gset_se <- as(new_gset,"SummarizedExperiment")
  rowData(gset_se)$PROBEID <- rownames(gset_se)#rowData(gset_se)$ID
  rowData(gset_se)[ID_NAME] <-  unlist(lapply(rowData(gset_se)$PROBEID,function(x){unlist(probe_gene_map[x][[1]])}))
  print(rowData(gset_se)$PROBEID)
  print(rowData(gset_se)[ID_NAME])
  gset_se_gene <-probe2gene(gset_se,to=ID_NAME)
  #for the sake of simplicity i convert back to the original  datatype 
  mapped_back <- as(gset_se_gene,"ExpressionSet")
  mapped_back
} 
