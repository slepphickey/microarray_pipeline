# This is the second script used to analyze microarray data

#time the script
tic <- as.integer(as.POSIXct( Sys.time() ))

#load the libraries
library(GEOquery)
library(GEOmetadb)
library(stringr)
library(tidyverse)

args <- commandArgs(TRUE)
# args[1] is the file path to the database
# args[2] is the GSE(s) for the experiments containing the GSMs you want to download and compare in small quotes
# for one GSE: "'GSE69078'"
# for more than one GSE put in double quotes and comma separate: "'GSE69078','GSE69079'" 
# args[3] path to the dirctory to save the sample table and downloaded CEL files to
# (needs to be high I/O compatiable, don't have trailing forward slash)


con <- dbConnect(SQLite(),args[1])
GSEs = args[2]

# make the directory to write the files to (needs to be high I/O capable)
dirname <- paste0(args[3], "/downloaded-files")
if(!dir.exists(dirname)) {
  dir.create(dirname)
}

rs <- dbGetQuery(con,paste0("select series_id, ",
                           "gsm, ",
                           "gpl, ",
                           "title, ",
                           "submission_date, ",
                           "source_name_ch1, ",
                           "characteristics_ch1, ",
                           "organism_ch1, ",
                           "supplementary_file ",
                           "from gsm where ",
                           "series_id IN (", 
                            GSEs,
                           ")"
                            )
                 )

rs$title <- sapply(rs$title, function(x) gsub("\t", " ", x))
rs$submission_date <- sapply(rs$submission_date, function(x) gsub("\t", " ", x))
rs$source_name_ch1 <- sapply(rs$source_name_ch1, function(x) gsub("\t", " ", x))
rs$characteristics_ch1 <- sapply(rs$characteristics_ch1, function(x) gsub("\t", " ", x))
rs$supplementary_file <- sapply(rs$supplementary_file, function(x) gsub("\t", " ", x))

if(nrow(rs) < 1) {stop("no GSMs for the GSEs found in the database")}

rs %>% group_by(series_id) %>% summarise(nGSMs = length(gsm))
rs[,"CELfile"] <- NA

celfiles = {}
for (i in 1:nrow(rs)) {
  
  gsmid <- rs[i,2]; cat(i, "\t", gsmid, "\n")
  gseid <- rs[i,1]
  
  dirname <- paste0(args[3], "/downloaded-files/",gseid)
  if(!dir.exists(dirname)) {
    dir.create(dirname)
  }
  
  gsm_file <- paste(gseid, "_", gsmid, sep="")
  print(gsm_file)
  
  url <- as.vector(rs[i,9])
  # sometimes there are multiple parts to the URL, select the part that ends in CEL.gz
  # print(url)
  if(grepl(" ", url)) {
    suburl <- unlist(strsplit(url, " "))
    url <- suburl[grep("cel.gz", suburl, ignore.case=TRUE)][1]
  }
  check_dwld <- tryCatch(download.file(url, quiet=TRUE, method="wget", destfile=paste(args[3], '/downloaded-files/', gseid, '/', gsm_file, ".CEL.gz", sep="")),
                         error=function(e) 1, warning=function(w) 2)
  if(check_dwld != 0) { 
    print("Unable to download file")
    rs[i,"CELfile"] = "not downloaded"
    next 
  }
  
  rs[i,"CELfile"] = paste0(args[3], '/downloaded-files/', gseid, '/', gsm_file, ".CEL.gz")
  
}  

GSE_names = unique(rs$series_id)
write.csv(rs, file = paste0(args[3], 
                            '/downloaded-files/', 
                            paste(GSE_names, collapse = "_"), 
                            "_sample_table.csv"), 
          row.names = F)

toc <- as.integer(as.POSIXct( Sys.time() ))
print(paste('The time it took in minutes to run the script was',(toc-tic)/60,sep=" "))




