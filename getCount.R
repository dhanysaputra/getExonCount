library(Sequgio)
library('RSQLite')
con <- dbConnect(dbDriver("SQLite"), "hg18.sqlite")
line=read.table('reg.chr17', sep="\t", fill=TRUE, as.is=TRUE, comment.char="", quote="")
Recurse<-function(a,n1,n2) { if (n2-n1==1 || n2==n1) {return(n1)} else { if (a<text[as.integer((n2-n1)/2)+as.integer(n1),3]) { return(Recurse(a,n1,as.integer((n2-n1)/2)+n1)) } else {return(Recurse(a,as.integer((n2-n1)/2)+n1,n2))} }}
j=17
text=dbGetQuery(con, paste("SELECT region_id,exon_name,start,cigar,end FROM hg18 where Chromosome='\"chr", "\"' order by start", sep=as.character(j)))
idx=Recurse(as.integer(40344263),0,dim(text)[1])
pos <- lapply(as.integer(line$V4),Recurse,n1=0,n2=dim(text)[1])
checking<-function(position,coordinate) { while(text[position+1,3]==text[position,3]){position=position+1}; if(coordinate+5>text[position+1,3] && coordinate+5<text[position,3] && coordinate+5<text[position,5]) {position=position+1}; return(position)}
newpos <- lapply(pos,checking,coordinate=line$V4)
for(i in 1:22){
	header<-''
	exon <- c()
	exon_txt<-''
	gene_txt<-''
	oldpos=0
	for(i in 1:dim(line)[1]){
		if(header!=strsplit(line[i,1], "\\.")[[1]][2]) {
			header<-strsplit(line[i,1], "\\.")[[1]][2]
			exon_txt<-''
			gene_txt<-''
			oldpos=0
			if ((!grepl("N", text[newpos[[i]][1],4]) && !grepl("N",line[i,6])) || (grepl("N", text[newpos[[i]][1],4]) && grepl("N", line[i,6]) && strsplit(strsplit(text[newpos[[i]][1],4], "M")[[1]][2], "N")[[1]][1]==strsplit(strsplit(line[i,6], "M")[[1]][2], "N")[[1]][1])){
				exon_txt<-text[newpos[[i]][1],2]
				gene_txt<-text[newpos[[i]][1],1]
				oldpos<-line[i,4]
			}
		}
		else {
			if ((!grepl("N", text[newpos[[i]][1],4]) && !grepl("N",line[i,6]) && gene_txt==text[newpos[[i]][1],1]) || (grepl("N", text[newpos[[i]][1],4]) && grepl("N", line[i,6]) && strsplit(strsplit(text[newpos[[i]][1],4], "M")[[1]][2], "N")[[1]][1]==strsplit(strsplit(line[i,6], "M")[[1]][2], "N")[[1]][1] && gene_txt==text[newpos[[i]][1],1])) {
				if(as.integer(oldpos)<as.integer(line[i,4])){
					temp <- paste(noquote(exon_txt), '.', noquote(text[newpos[[i]][1],2]), '__', as.character(gene_txt), sep='')
					if (!temp %in% names(exon)) {exon[temp]=1} else {exon[temp]<-noquote(unname(exon[temp]))+1}
				} else {
					temp <- paste(noquote(text[newpos[[i]][1],2]), '.', noquote(exon_txt), '__', as.character(gene_txt), sep='')
					if (!temp %in% names(exon)) {exon[temp]=1} else {exon[temp]<-noquote(unname(exon[temp]))+1}
				}
			}
		}
	}
}