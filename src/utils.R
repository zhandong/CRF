
### read in a TCGA directory 
read.tcga.rna <- function(fdir)
{
#fdir  = "./TCGA-GBM/Expression-Genes/UNC__AgilentG4502A_07_1/Level_3/"

gnames =c()
cnames=c()
tcgaData = c()

for(i in dir(fdir,full.names=T))
{
	tmp =scan(i,what="string",skip=1)
	tmp = matrix(tmp,ncol=3,byrow=T)
	if(sum(is.na(as.numeric(tmp[,3])))){print(i)}
	if(length(gnames)==0){
		gnames = tmp[,2]
		tcgaData = cbind(tcgaData,as.numeric(tmp[,3]))
		cnames = unique(tmp[,1])
	}else{
		if(sum(gnames!=tmp[,2])){
			print(c("wrong order",i))
		}
		tcgaData = cbind(tcgaData,as.numeric(tmp[,3]))
		cnames = c(cnames,unique(tmp[,1]))
		}
}
rownames(tcgaData)=gnames
colnames(tcgaData)=cnames
return(tcgaData)

}


##### read the cgh 


### read in a TCGA directory 
read.tcga.cgh <- function(fdir)
{
#fdir  = "./TCGA-GBM/Expression-Genes/UNC__AgilentG4502A_07_1/Level_3/"

tcgaData = c()

for(i in dir(fdir,full.names=T))
{
	print(i)
	tmp =scan(i,what="character",skip=1)
	tmp = matrix(tmp,ncol=6,byrow=T)
	tcgaData = rbind(tcgaData,tmp)
}

cghData = data.frame(
	ID = as.character(tcgaData[,1]),
	chrom = as.character(tcgaData[,2]),
	loc.start = as.integer(tcgaData[,3]),
	loc.end = as.integer(tcgaData[,4]),
	num.mar = as.integer(tcgaData[,5]),
	seg.mean = as.numeric(tcgaData[,6]),
	stringsAsFactors=FALSE
)

cghData[is.na(cghData)]=0

library(CNTools)
data('geneInfo')
seg = CNSeg(cghData)
mat= getRS(seg,by="gene",imput=FALSE,XY=FALSE,geneMap=geneInfo,what="median")
return(mat)
}






