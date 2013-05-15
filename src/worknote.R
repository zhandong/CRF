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
index2 = rowSds(tcgaData$rna) > quantile(rowSds(tcgaData$rna),0.95,na.rm=T)
dataY = tcgaData$rna[index2,]
dataY = dataY[!is.na(rownames(dataY)),]

save(dataX,dataY,file="../dat/XY-GBM-TCGA.rdata")

####
rm(list=ls())
load('../dat/XY-GBM-TCGA.rdata')
source('CRFfunctions.R')

GBMNetwork=GGM.linearCRF.path.network.stars(t(dataY),t(dataX),20,0.05,niter=100)






