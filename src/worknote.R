source('CRFfunctions.R')

Y = matrix(rnorm(100),nrow=20)
X = matrix(rnorm(60),nrow=20)

tmp=GGM.linearCRF.path.network(Y,X,20)
est=GGM.linearCRF.path.network.stars(Y,X,20,0.05)

