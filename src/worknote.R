source('CRFfunctions.R')

Y = matrix(rnorm(1000),nrow=200)
X = matrix(rnorm(600),nrow=200)

tmp=GGM.linearCRF.path.network(Y,X,20)
est=GGM.linearCRF.path.network.stars(Y,X,20,0.05)

############ load real data. 

#### preprocess GBM