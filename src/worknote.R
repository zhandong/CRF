rm(list=ls())
source('CRFfunctions.R')

# # Y = matrix(rnorm(1000),nrow=200)
# X = matrix(rnorm(600),nrow=200)

# tmp=GGM.linearCRF.path.network(Y,X,20)
# est=GGM.linearCRF.path.network.stars(Y,X,20,0.05)

############ load real data. 
source('utils.R')
library(CNTools)
data('geneInfo')
#### preprocess GBM CGH
fdir1 = "../dat/TCGA-GBM/CNV_Array/HMS__HG-CGH-415K_G4124A/Level_3/"
fdir2 = "../dat/TCGA-GBM/CNV_Array/MSKCC__HG-CGH-244A/Level_3/"

mat1 = read.tcga.cgh(fdir1)
mat2 = read.tcga.cgh(fdir2)

##### preprocess GBM array RNA

fdir1 = '../dat/TCGA-GBM/Expression-Genes/BI__HT_HG-U133A/Level_3/'
fdir2 = '../dat/TCGA-GBM/Expression-Genes/UNC__AgilentG4502A_07_1/Level_3/'
fdir3 = '../dat/TCGA-GBM/Expression-Genes/UNC__AgilentG4502A_07_2//Level_3/'

rna.mat1 = read.tcga.rna(fdir1)
rna.mat2 = read.tcga.rna(fdir2)
rna.mat3 = read.tcga.rna(fdir3)

save.image('../dat/tcga-ggm-dat.rdata')

##### match sample ID

load("../dat/tcga-ggm-dat.rdata")
sampleID=c()
sampleID$rna=c(colnames(rna.mat2),colnames(rna.mat3))
sampleID$rna=substring(sampleID$rna,1,15)

tcgaData=c()
tcgaData$rna = cbind(rna.mat2,rna.mat3)
colnames(tcgaData$rna)=sampleID$rna
tcgaData$rna = tcgaData$rna[,!duplicated(sampleID$rna)]


sampleID$dna=c(colnames(rs(mat1)),colnames(rs(mat2)))
sampleID$dna=substring(sampleID$dna,1,15)
tcgaData$dna = cbind(rs(mat1),rs(mat2))
colnames(tcgaData$dna)=sampleID$dna
tcgaData$dna = tcgaData$dna[,!duplicated(sampleID$dna)]

commonID = intersect(colnames(tcgaData$dna),colnames(tcgaData$rna))
tcgaData$dna = tcgaData$dna[,commonID]
tcgaData$rna = tcgaData$rna[,commonID]

###### Grab the EGFR, CDKN2A/B PTEN, CDK4, PDGFRA
covNames = c("EGFR","CDKN2A", "PTEN","CDK4","PDGFRA")
tmp = rs(mat1)
index = which(tmp[,5]%in%covNames,arr.ind=T)

dataX=tcgaData$dna[index,]
dataX=dataX>log2(3/2) | dataX<log2(1/2)

### filtering Y by variance
index2 = rowSds(tcgaData$rna) > quantile(rowSds(tcgaData$rna),0.90,na.rm=T)
dataY = tcgaData$rna[index2,]
dataY = dataY[!is.na(rownames(dataY)),]

pathwayGenes = unique(scan('../dat/c2.cp.v3.1.symbols.gmt.txt',what="character"))
dataY = dataY[rownames(dataY)%in% pathwayGenes,]
dataY = rbind(dataY,tcgaData$rna[covNames,])

save(dataX,dataY,file="../dat/XY-GBM-TCGA.rdata")

####
rm(list=ls())
load('../dat/XY-GBM-TCGA.rdata')
source('CRFfunctions.R')
library(igraph)

GBMNetwork=GGM.linearCRF.path.network.stars(t(dataY),t(dataX),20,0.05,niter=100)
for(i in 1:5){
#	tmp =matrix(   GBMNetwork$path$Ahat[,i,12],nrow=nrow(dataY))
#	rownames(tmp)=rownames(dataY)
#	colnames(tmp)=rownames(dataY)
graphObj=graph.adjacency( matrix(   GBMNetwork$path$Ahat[,i,12],nrow=nrow(dataY)),mode="undirected")
V(graphObj)$name= rownames(dataY)
write.graph(graphObj,file=paste("../dat/GBM_",i,".gml",sep=""),format="gml")
	
}

save.image("../dat/GBM-CRF-results.rdata")
#################


rm(list=ls())
load('../dat/GBM-CRF-results.rdata')
library(igraph)

#GBMNetwork$var<0.001
vd=c()
for(i in 1:5){
	tmp =matrix(GBMNetwork$path$Ahat[,i,12],nrow=nrow(dataY))
	vd=rbind(vd,rowSums(tmp))
}

which(colSums(vd)==max(colSums(vd)),arr.ind=T)

tnames=c("PTEN","CDK4","PDGFR4","EGFR","CDKN2A/B")
for( hubID in order(colSums(vd),decreasing=T)[1:50]){
	
	pdf(file=paste('../pdfs/GBM_GGM_',rownames(dataY)[hubID],".pdf",sep=""),width=10,height=3)
	par(mfrow=c(1,5))

	for(i in 1:5){
	tmp =matrix(GBMNetwork$path$Ahat[,i,12],nrow=nrow(dataY))
	graphObj=graph.adjacency( matrix(   GBMNetwork$path$Ahat[,i,12],nrow=nrow(dataY)),mode="undirected")
	V(graphObj)$name= rownames(dataY)
	subgraphObj= induced.subgraph(graphObj,c(hubID,which(tmp[hubID,]==1,arr.ind=T)))
	plot(subgraphObj,vertex.color="green",vertex.size=50)
	title(main=tnames[i])
	}
		dev.off()

}



