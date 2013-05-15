####################
#functions for conditional random fields

require('huge')
require('glmnet')
require('multicore')



########internal function
mykron = function(vec,ind){kronecker(vec[ind],vec[-ind])}


#################
#neighborhood selection for Gaussian linear combo model
#assumes Y is n x p and X is n x q & node is the index of the node of interest
#each of the covariates X should be of the same type
#as the corresponding coeffiicents are regularized by the same amount
#Returns:
#lams = sequence of lambda values
#CoefMat = the full coefficient matrix - (p-1)*q x nlams
#Bhat = coefficients arranged as a 3D array: (p-1) x q x nlams
GGM.linearCRF.path.neighborhood = function(Y,X,node,nlams)
{
  ind = (ncol(X)+1):(ncol(X)+ncol(Y)-1)
  XY = t(apply(cbind(X,Y[,-node]),1,mykron,ind))
  fit = glmnet(XY,Y[,node],family="gaussian",nlambda=nlams)
#  return(list(lams=fit$lambda,Coefhat=fit$beta,Bhat=array(fit$beta,c(ncol(Y)-1,ncol(X),nlams))))  
  output =list(lams=fit$lambda,Coefhat=fit$beta,Bhat=array(fit$beta,c(ncol(Y)-1,ncol(X),nlams)))
  tmp= array(0,dim=c(ncol(Y),ncol(X),nlams))
  if(node ==1){
    tmp[(node+1):ncol(Y),,]= output$Bhat
  } else if(node == ncol(Y)){
  	tmp[1:(node-1),,]=output$Bhat
  } else{
    tmp[1:(node-1),,]=output$Bhat[1:(node-1),,]
    tmp[(node+1):ncol(Y),,]= output$Bhat[node:(ncol(Y)-1),,]
  }
  output$Bhat = tmp;
  return(output)
}



#################
#neighborhood selection for Gaussian linear combo model
#assumes Y is n x p and X is n x q & node is the index of the node of interest
#each of the covariates X should be of the same type
#as the corresponding coeffiicents are regularized by the same amount
#Returns:
#Bhat = coefficients arranged as a 3D array: p^2 x q x nlams
#Ahat = adjcency matrix a 3D array: p^2 x q x nlams

GGM.linearCRF.path.network = function(Y,X,nlams,opt ="and")
{	
	p=ncol(Y)
	q=ncol(X)
	Bhat = array(0,dim=c(p^2,q,nlams))
	
	
	for(i in 1:ncol(Y))
	{
		nodei= GGM.linearCRF.path.neighborhood(Y,X,node=i,nlams=nlams)	
		index = (1:p)+(i-1)*p
		Bhat[index,,] = nodei$Bhat
	}
	
	Ahat = Bhat
	for(i in 1:q){
		for(j in 1:nlams){
			tmp = matrix(Bhat[,i,j],nrow=p,ncol=p,byrow=F)
			if(opt == "and"){
				tmp = sign(abs(tmp)>0.0000001)
				tmp = sign(tmp*t(tmp))
			}else if (opt == "or"){
				tmp = sign(abs(tmp)>0.0000001)
				tmp = sign(tmp+t(tmp))
			}
			Ahat[,i,j]=c(tmp)		
		}
	}
	return(list(Bhat=Bhat,Ahat=Ahat))
}


#################
#star selection for Gaussian linear combo model
#assumes Y is n x p and X is n x q & node is the index of the node of interest
#each of the covariates X should be of the same type
#as the corresponding coeffiicents are regularized by the same amount
#Returns:
#Bhat = coefficients arranged as a 3D array: p^2 x q x nlams

GGM.linearCRF.path.network.stars = function(Y,X,nlams,stars.thresh,niter=100)
{	
	p=ncol(Y)
	q=ncol(X)
	if(nrow(Y)==nrow(X)){
		n = nrow(Y)
	}else{
		print("Number of rows is inconsistent!!!!!")	
	}
	est=c()
	est$merge=list()
	for(i in 1:niter){
		 mes <- paste(c("Conducting Subsampling....in progress:", 
                    floor(100 * i/niter), "%"), collapse = "")
                  cat(mes, "\r")
                  flush.console()

		
		if(n>144){
			subsample.ratio= 10*sqrt(n)/n
		}else{
			subsample.ratio = 0.8
		}
		index = sample(1:nrow(Y),size= floor(n*subsample.ratio),replace=FALSE)
		tmp=GGM.linearCRF.path.network(Y[index,],X[index,],nlams=nlams)
		for(j in 1:nlams){
			if(i ==1){
				est$merge[[j]] =   tmp$Ahat[,,j]/niter 
			}else{
				est$merge[[j]] =   tmp$Ahat[,,j]/niter + est$merge[[j]]
			}
		}
	}
	
	
	for(j in 1:nlams){
		tmp = est$merge[[j]]
		index = 1:p+(0:(p-1))*p
		tmp= tmp[-index,]
		est$var[j] = 4*sum(tmp*(1-tmp))/(p*(p-1)*q)	
	}
	est$opt.index = max(which.max(est$var >=  stars.thresh)[1] - 1, 1)
	est$path = GGM.linearCRF.path.network(Y,X,nlams=nlams)
	
	return(est)
}







#################
#neighborhood selection for Gaussian mean model
#assumes Y is n x p and X is n x q & node is the index of the node of interest
#each of the covariates X should be of the same type
#as the corresponding coeffiicents are regularized by the same amount
#alpha is a scalar that differentially penalizes the coefficients corresponding to X
#penalty term = \lambda*\alpha ||\beta||_1 + \lambda||\theta||_1
#Returns:
#lams = sequence of lambda values
#CoefMat = the full coefficient matrix - (p-1 + q) x nlams
#Bhat = coefficients corresponding to X's - q x nlams
#Thethat = coefficients corresponding to Y's - (p-1) x nlams
GGM.meanCRF.path.neighborhood = function(Y,X,node,nlams,alpha)
{
  pfac = c(rep(1,(ncol(Y)-1)),rep(alpha,ncol(X)))
  fit = glmnet(cbind(Y[,-node],X),Y[,node],family="gaussian",nlambda=nlams,penalty.factor=pfac)
  return(list(lams=fit$lambda,CoefMat=fit$beta,Thethat=fit$beta[1:(p-1),],Bhat=fit$beta[p:(p-1+q),]))    
}


#################
#neighborhood selection for Poisson mean model
#assumes Y is n x p and X is n x q & node is the index of the node of interest
#each of the covariates X should be of the same type
#as the corresponding coeffiicents are regularized by the same amount
#alpha is a scalar that differentially penalizes the coefficients corresponding to X
#penalty term = \lambda*\alpha ||\beta||_1 + \lambda||\theta||_1
#Note that both B and Theta are non-positive
#Returns:
#lams = sequence of lambda values
#CoefMat = the full coefficient matrix - (p-1 + q) x nlams
#Bhat = coefficients corresponding to X's - q x nlams
#Thethat = coefficients corresponding to Y's - (p-1) x nlams
PGM.meanCRF.path.neighborhood = function(x,y,node,alpha,nlams,lmax,startb=0)
{
  X = cbind(y[,-node],x); Y = y[,node];
  pfac = c(1,rep(1,(ncol(y)-1)),rep(alpha,ncol(x)))
  n = nrow(X); p = ncol(X);
  lams = exp(seq(log(lmax),log(.0001),l=nlams));
  if(nlams==1){lams = lmax};
  thr = 1e-6; maxit = 1e6;
  Xt = cbind(t(t(rep(1,n))),X);
  L = max(eigen(t(X)%*%X,only.values=TRUE)$values);
  alphas = 0; Bmat = matrix(0,p,nlams);
  if(sum(startb)==0){Bhat = matrix(runif(p+1),p+1,1); Bhat[1] = 0}else{Bhat=startb}
  for(i in 1:nlams){
    iter = 1; ind = 1;
    while(thr<ind & iter<maxit){
      oldb = Bhat;
      tmp = Bhat - (t(Xt)%*%Y - t(Xt)%*%exp(-Xt%*%Bhat))/L;
      Bhat = matrix(sapply(tmp - pfac*lams[i]/L,max,0),p+1,1); Bhat[1] = tmp[1];
      ind = sum((Bhat - oldb)^2);
      iter = iter + 1;
    }
    alphas[i] = -Bhat[1];
    Bmat[,i] = -Bhat[2:(p+1),drop=FALSE];
  }
  return(list(ahat=alphas,CoefMat=Bmat,lams=lams,Thethat=Bmat[1:(ncol(y)-1),],Bhat=Bmat[ncol(y):(ncol(y)-1+ncol(x)),]))
}



#################
#neighborhood selection for Poisson linear combo model
#assumes Y is n x p and X is n x q & node is the index of the node of interest
#each of the covariates X should be of the same type
#as the corresponding coeffiicents are regularized by the same amount
#Note that both all coefficients as non-positive
#Returns:
#lams = sequence of lambda values
#CoefMat = the full coefficient matrix - (p-1)*q x nlams
#Bhat = coefficients arranged as a 3D array: (p-1) x q x nlams
PGM.linearCRF.path.neighborhood = function(x,y,node,nlams,lmax,startb=0)
{
  ind = (ncol(x)+1):(ncol(x)+ncol(y)-1)
  X = t(apply(cbind(x,y[,-node]),1,mykron,ind))
  Y = y[,node]
  n = nrow(X); p = ncol(X);
  lams = exp(seq(log(lmax),log(.0001),l=nlams));
  if(nlams==1){lams = lmax};
  thr = 1e-6; maxit = 1e6;
  Xt = cbind(t(t(rep(1,n))),X);
  L = max(eigen(t(X)%*%X,only.values=TRUE)$values);
  alphas = 0; Bmat = matrix(0,p,nlams);
  if(sum(startb)==0){Bhat = matrix(runif(p+1),p+1,1); Bhat[1] = 0}else{Bhat=startb}
  for(i in 1:nlams){
    iter = 1; ind = 1;
    while(thr<ind & iter<maxit){
      oldb = Bhat;
      tmp = Bhat - (t(Xt)%*%Y - t(Xt)%*%exp(-Xt%*%Bhat))/L;
      Bhat = matrix(sapply(tmp - lams[i]/L,max,0),p+1,1); Bhat[1] = tmp[1];
      ind = sum((Bhat - oldb)^2);
      iter = iter + 1;
    }
    alphas[i] = -Bhat[1];
    Bmat[,i] = -Bhat[2:(p+1),drop=FALSE];
  }
  return(list(ahat=alphas,CoefMat=Bmat,lams=lams,Bhat=array(Bmat,c(ncol(y)-1,ncol(x),nlams))))
}
