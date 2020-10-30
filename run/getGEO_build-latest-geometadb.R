# This is the first script used to get microarray or RNA-seq data
# this will make a local copy of the database for GEO
# This does not require a qsub and only takes a few minutes to run

# load the library
library(GEOmetadb)
# get the date in a good format to save database with as database changes over time
today <- Sys.Date()
gooddate <- format(today,format="%Y-%m-%d")
# filename to save the database as
db_filename <- paste(gooddate,"_GEOmetadb",".sqlite.gz",sep='')
# make a directory that will eventually store the downloaded gms files
dirname <- paste0(gooddate,"_downloaded-files")
# check if file exists and if not make a new file
if(!file.exists(db_filename)) {
	getSQLiteFile(destfile = db_filename)
} else {
	print('The file already exists')
}
# make a directory that will store download gsm files
if(!dir.exists(dirname)) {
	dir.create(dirname)
}



