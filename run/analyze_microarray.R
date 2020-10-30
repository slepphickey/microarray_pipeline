# This is the third script used to analyze microarray data
# time the script
tic <- as.integer(as.POSIXct( Sys.time() ))

#load the libraries
library(oligo)
library(stringr)
library(affycoretools)
library(tidyverse)
library(purrr)
library(qvalue)


args <- commandArgs(TRUE)
# args[1] is the path to the directory containing the downloaded-files dir
# args[2] name for the output folder and output file prefixes
# args[3] is the path to a sample_table.csv file containing rows for the samples to be included in DE analysis
# args[4] is the path to BrainArrayMapping_withGLP.csv
# args[5] is  a comma separated list of the groups each sample belongs to [optional]
# ie for 6 samples "treated, treated, treated, untreated, untreated, untreated"
# if args[5] is not included this information must be added to arg[3] indepenently as a column named "Condition"

# this is a function that used biomart to map mouse genes to human genes one-to-one homologs
# https://www.r-bloggers.com/2016/10/converting-mouse-to-human-gene-names-with-biomart-package/

convertMouseGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("mgi_symbol"), 
                   filters = "mgi_symbol", 
                   values = x , 
                   mart = mouse, 
                   attributesL = c("hgnc_symbol"), 
                   martL = human, 
                   uniqueRows=T)
  
  return(genesV2)
}

prefix = args[2]

dir.name <- paste0(args[1], "/", prefix)
if(!dir.exists(dir.name)) {
  dir.create(dir.name)
}

meta = read.csv(args[3])
BAtoGPL = read.csv(args[4])

if (length(args < 5)) {
  Conditions = unique(meta$Condition)
} else {
  cond = unlist(str_split(args[5], ", "))
  meta$Condition  = cond
}

# https://www.biostars.org/p/283083/#283261
# https://www.biostars.org/p/266507/

# 1) if the GSMs are from the same GSE and platform, normalize together and run through limma
# Test: GSE69078
# 2) if the GSMs are from different GSE and the same platform, normalize together and include GSE as a covariate
# Test: GSE130459, GSE12968
# 3) if GSMs are from different platforms, same species, normalize separately, convert to zscore, merge by shared genes, and include platform as a covariate
# Test GSE131978
# 4) if GSMs are from different species, normalize separately, convert to zscore, map one-to-orthos, merge by shared genes, and include species/GSE as a covariate
# test GSE46301
# 5) GSE with an array type we can't map 
# Test GSE46806
# if batch overlaps with condition throw a warning
# for all make PCA plot with GSE/platform/and condition labels

Conditions = unique(meta$Condition)
GSEs = unique(meta$series_id)
GPLs = unique(meta$gpl)

dataRMA = list()
organism = {}

for (GPL in GPLs){
  if(!GPL %in% BAtoGPL$gpl) {stop(paste("No BrainArray annotation is currently associated with", 
                                        GPL, 
                                        "go to http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/24.0.0/entrezg.asp and see if the chip associated with",
                                        GPL,
                                        "is present under 'Chip'. You can manually add the pdName and dbName name to BrainArrayMapping_withGLP.csv and try again.",
                                        "pdName is the Custom CDF name with '.' instead of '_', and all lowercase.",
                                        "dbName is pdName without '.'"))}
  
  pdName = BAtoGPL$pdName[which(BAtoGPL$gpl == GPL)]
  dbName = BAtoGPL$dbName[which(BAtoGPL$gpl == GPL)]
  organism = c(organism, BAtoGPL$organism[which(BAtoGPL$gpl == GPL)])
  
  install.packages(paste0("http://mbni.org/customcdf/24.0.0/entrezg.download/pd.", 
                          pdName, 
                          "_24.0.0.tar.gz"), 
                   repos = NULL, 
                   type = "source")
  
  install.packages(paste0("http://mbni.org/customcdf/24.0.0/entrezg.download/", 
                          dbName,
                          ".db_24.0.0.tar.gz"), 
                   repos = NULL, 
                   type = "source")
  
  cellFiles = meta$CELfile[which(meta$gpl == GPL)]

  dat <- read.celfiles(cellFiles, pkgname = paste0("pd.", pdName))
  dataRMA[[GPL]] <- rma(dat)
  
  library(paste0(dbName, ".db"), character.only = TRUE)
  dataRMA[[GPL]] <- annotateEset(dataRMA[[GPL]], paste0(dbName, ".db"))
  
}

if(length(GPLs) > 1){
  
  rmaDF = list()
  for(i in 1:length(dataRMA)){
   
    rmaDF[[i]] = exprs(dataRMA[[i]])
    gene_ids = fData(dataRMA[[i]])
    rownames(rmaDF[[i]]) = gene_ids$SYMBOL
    rmaDF[[i]] = rmaDF[[i]][!is.na(rownames(rmaDF[[i]])),]
    
    if(length(unique(organism)) == 2 & 
       "Mus musculus" %in% organism & 
       "Homo sapiens" %in% organism &
       organism[i] == "Mus musculus"){
      
      print("converting mouse gene symbols to human gene symbols -- one-to-one orthologs only")
      mouse_genes = rownames(rmaDF[[i]])
      mouse2human = convertMouseGeneList(mouse_genes)
      # remove mouse genes that map to more than one human gene
      mouse_dup = unique(mouse2human$MGI.symbol[duplicated(mouse2human$MGI.symbol)])
      mouse2human = mouse2human[!mouse2human$MGI.symbol %in% mouse_dup,]
      # remove human genes that map to more than one mouse gene
      hum_dup = unique(mouse2human$HGNC.symbol[duplicated(mouse2human$HGNC.symbol)])
      mouse2human = mouse2human[!mouse2human$HGNC.symbol %in% hum_dup,]
      mouse2human = mouse2human[order(mouse2human$MGI.symbol),]
      # filter for one-to-one matchers
      keep_rows = mouse2human[,1]
      rmaDF[[i]] = rmaDF[[i]][keep_rows,]
      rownames(rmaDF[[i]]) = mouse2human$HGNC.symbol
      rmaDF[[i]] = as.data.frame(rmaDF[[i]])
      rmaDF[[i]]$SYMBOL = rownames(rmaDF[[i]])
      print(paste(nrow(rmaDF[[i]]), "one-to-one human to mouse orthologs kept out of", length(mouse_genes), "original mouse genes"))
      
      } else if (length(unique(organism)) > 3) {
        stop("This script can only combine human and mouse samples right now. Remove other species and process independently")
      } else {
        rmaDF[[i]] = as.data.frame(rmaDF[[i]])
        rmaDF[[i]]$SYMBOL = rownames(rmaDF[[i]])}
    }
  
  finalDF = rmaDF %>%
    purrr::reduce(inner_join, by = "SYMBOL")
  
  rownames(finalDF) = finalDF$SYMBOL
  finalDF$SYMBOL = NULL
  
  } else {
    dataRMA = dataRMA[[1]]
    rmaDF = exprs(dataRMA)
    gene_ids = fData(dataRMA)
    rmaDF = rmaDF[!is.na(rownames(rmaDF)),]
    finalDF = as.data.frame(rmaDF)
}

write.csv(finalDF, file = paste0(dir.name,"/", prefix, "_RMA_df.csv"))
print(paste0(prefix, "_RMA_df.csv saved"))

# DE analysis with limma
# change finalDF back to mat?
library(limma)

meta$CELfile = basename(meta$CELfile)
finalDF = finalDF[,meta$CELfile]

f = factor(meta$Condition, levels = Conditions)

print(paste("data comes from", length(GSEs), "GSE(s)"))
print(paste("data comes from", length(GPLs), "platform(s)"))

if(length(GSEs) > 1 & length(GPLs) > 1) {
  design.mat = model.matrix(~0 + f + series_id + gpl, meta)
  print("removing batch effects caused by GSE and platform")
  
  } else if (length(GSEs) > 1 & length(GPLs) == 1) {
    design.mat = model.matrix(~0 + f + series_id, meta)
    print("removing batch effect caused by GSE")
    
    } else if (length(GSEs) == 1 & length(GPLs) > 1) {

      print("removing batch effect caused by platform")
      design.mat = stats::model.matrix(~0 + f + gpl, meta)
      
     } else {
        design.mat = model.matrix(~0 + f, meta)
        print("no batch effect caused by GSE or platform")
        }

colnames(design.mat)[1:length(Conditions)] = Conditions
print("design.mat done")
data.fit = lmFit(finalDF,design.mat) 
mean_exp = data.fit$coefficients
write.csv(mean_exp, file = paste0(dir.name,"/", prefix, "_batch_corrected_mean_log_expression.csv"))
print(paste0(prefix, "_batch_corrected_mean_log_expression.csv saved"))

compare = {} 
  for(i in 1:(length(Conditions)-1)){
    for(j in 2:length(Conditions)){
      if(i == j) {next}
      to_add = paste0(Conditions[i],"-",Conditions[j])
      compare = c(compare, to_add)
    }
  }

print(paste("comparing these groups:", compare))

contrast.matrix = makeContrasts(contrasts = compare, levels=design.mat)
data.fit.con = contrasts.fit(data.fit,contrast.matrix)
data.fit.eb = eBayes(data.fit.con)

logfc = data.fit.eb$coefficients
colnames(logfc) = paste0("logFC_", colnames(logfc))
pairwise_pval = data.fit.eb$p.value
colnames(pairwise_pval) = paste0("pairwise.t.pval_", colnames(pairwise_pval))

if(length(compare) > 1) {
  anova_pval = as.data.frame(data.fit.eb$F.p.value)
  colnames(anova_pval) = "ANOVA.F.pval"
  results = as.data.frame(cbind(logfc, pairwise_pval, anova_pval))
  results$ANOVA.F.qval = qvalue(results$ANOVA.F.pval)$qvalues
} else {
    results = as.data.frame(cbind(logfc, pairwise_pval))
    results$qval = qvalue(results[,2])$qvalues
}

write.csv(mean_exp, file = paste0(dir.name,"/", prefix, "_de_anaylsis.csv"))
print(paste0(prefix, "_de_anaylsis.csv saved"))

toc <- as.integer(as.POSIXct( Sys.time() ))
print(paste('The time it took in minutes to run the script was',(toc-tic)/60,sep=" "))










