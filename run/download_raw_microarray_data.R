# This is the second script used to analyze microarray data
# authors: Stephanie Hickey (original code), Kewalin Samart (2023 version)

#time the script
tic <- as.integer(as.POSIXct( Sys.time() ))

#load the libraries
library(GEOquery)
library(stringr)
library(tidyverse)

args <- commandArgs(TRUE)
# args[1] is the GSE(s) for the experiments containing the GSMs you want to download and compare in small quotes
# for one GSE: "'GSE69078'"
# for more than one GSE put in double quotes and comma separate: "'GSE69078','GSE69079'"


GSEs = strsplit(args[1], split = ",")[[1]] # list of GSE ids


#meta_TB = read.csv("/Users/kewalinsamart/Documents/GitHub/Drugrep_Tuberculosis/data/microarray/metadata/TB_gene_info.csv")
#GSEs = meta_TB$geo_id


# make the directory to write the files to (needs to be high I/O capable)
dirname_meta <- paste0(getwd(), "/metadata")
if(!dir.exists(dirname_meta)) {
  dir.create(dirname_meta)
}
dirname_raw <- paste0(getwd(), "/raw_data")
if(!dir.exists(dirname_raw)) {
  dir.create(dirname_raw)
}

# get dataset series matrix file
gse = list()
for(i in 1:length(GSEs)){
  gse_i <- getGEO(GSEs[[i]])
  gse[[i]] <- gse_i
}

# obtain list of metadata and save the files

rs_subset <- data.frame()
for(i in 1:length(gse)){
  pd <- pData(gse[[i]][[1]]) # extract phenotype data
  pd$series_id <- GSEs[i] # add series_id column -- in case of multiple GSE ids
  rs <- as.data.frame(pd)
  write.csv(rs,paste0(dirname_meta,"/",GSEs[[i]],"_meta.csv"))
  if(nrow(rs_subset) < 1){
    rs_subset <- rs[,c("series_id","geo_accession","platform_id","supplementary_file")]
  }else{
    rs_subset <- rbind(rs_subset,rs[,c("series_id","geo_accession","platform_id","supplementary_file")])
  }
}

if(nrow(rs_subset) < 1) {stop("no GSMs for the GSEs found in the database")}

rs_subset %>% group_by(series_id) %>% summarise(nGSMs = length(geo_accession))
rs_subset[,"CELfile"] <- NA

# create raw data folder for all the GSE ids
gseids <- unique(rs_subset$series_id)
gsmids <- rs_subset$geo_accession; cat(i, "\t", gsmid, "\n")
for(gseid in gseids){
  dirname <- paste0(dirname_raw, "/",gseid)
  if(!dir.exists(dirname)) {
    dir.create(dirname)
  }
}

for (i in 1:nrow(rs_subset)) {
  gsmid <- rs_subset$geo_accession[i]
  gseid <- rs_subset$series_id[i]
  gplid <- rs_subset$platform_id[i]
  gsm_file <- paste(gseid,"_",gplid,"_", gsmid, sep="")
  print(gsm_file)
  url = paste0("https://www.ncbi.nlm.nih.gov/geo/download/?acc=",gsmid,"&format=file&file=",gsmid,"%2ECEL%2Egz")
  rs_subset[i,"http"] = url
  check_dwld <- tryCatch(download.file(url, quiet=TRUE, method="wget", destfile=paste(dirname_raw, '/', gseid, '/', gsm_file, ".CEL.gz", sep="")),
                        error=function(e) 1, warning=function(w) 2)
  if(check_dwld != 0) {
    print("Unable to download file")
    rs_subset[i,"CELfile"] = "not downloaded"
    next
  }
  rs_subset[i,"CELfile"] = paste0(dirname_raw, '/', gseid, '/', gsm_file, ".CEL.gz")
}

write.csv(rs_subset, file = paste0(dirname_meta,'/marrayTB_sample_table.csv'),row.names = F)

toc <- as.integer(as.POSIXct( Sys.time() ))
print(paste('The time it took in minutes to run the script was',(toc-tic)/60,sep=" "))

