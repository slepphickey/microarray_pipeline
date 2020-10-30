# This is the first script used to get microarray or RNA-seq data
# this will make a local copy of the database for GEO
# This does not require a qsub and only takes a few minutes to run

args <- commandArgs(TRUE)
# args[1] path to the directory you want the data base in
# needs to be high I/O compatiable, don't have trailing forward slash)

# load the library
library(GEOmetadb)
# get the date in a good format to save database with as database changes over time
today <- Sys.Date()
gooddate <- format(today,format="%Y-%m-%d")
# filename to save the database as
db_filename <- paste(args[1], "/", gooddate,"_GEOmetadb",".sqlite.gz",sep='')
# make a directory that will eventually store the downloaded gms files

if(!file.exists(db_filename)) {
	getSQLiteFile(destfile = db_filename)
} else {
	print('The file already exists')
}




