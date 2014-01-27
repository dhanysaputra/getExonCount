library(Sequgio)
library(RSQLite)
dbfile <- "hg18GFF.sqlite"
mybio <- loadDb(dbfile)
mparam <- MulticoreParam(8)
txdb <- reshapeTxDb(mybio,probelen=50L,with.junctions=T,mcpar=mparam)
write.table(as.data.frame(txdb@unlistData), "txdb.sql", sep="\t")
db <- dbConnect(SQLite(), dbname="hg18.sqlite")
dbWriteTable(conn=db, name="hg18", value="/home/dhany/txdb.sql", row.names=FALSE, header=FALSE, sep="\t")
