# get the mapping of microarray chips to GPL
# chips as listed in BrainArray v24
# pdName is for installing the pd package from Brainarray
# dbName is for installing the db package from Brainarray
# http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/24.0.0/entrezg.asp
# load the library

library(GEOmetadb)
library(stringr)

BAM = read.delim("BrainArrayMapping.tsv")
chip = BAM$Chip
map = getBiocPlatformMap(con, bioc = chip)
# 45 chips have GPLs mapped to them

BAM$bioc_package = tolower(BAM$Chip)

BAM = merge(map, BAM, by = "bioc_package")
tosave = BAM[c(1,2,3,5,8,9,10)]

write.csv(tosave, file = "BrainArrayMapping_withGLP.csv", row.names = F)

