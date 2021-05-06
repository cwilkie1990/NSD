

run.NSDmodel<-function(nIter=100,nBurnIn=10,nChains=2,nThin=10,yData,xData,xPred,coordsData,coordsPred,
                        times.yData=NULL,
                        times.xData=NULL,
                        times.xPred=NULL,
                        times.yPred=NULL,
                        basis.type="bspline",
                        basis.dim=5,
                        period=1,
                        phiAlpha=0.1,phiBeta=0.1,
                        aAlpha=2,bAlpha=1,aBeta=2,bBeta=1,aY=2,bY=1,aC=2,bC=1,aX=2,bX=1,muD=NULL,SigmaD=NULL,
                        sigmaAlphaPrecInit=1,sigmaBetaPrecInit=1,sigmaYPrecInit=1,sigmaCPrecInit=1,alphaInit=NULL,betaInit=NULL,cInit=NULL,sigmaXPrecInit=1,dInit=NULL){
  if(!is.list(times.yData)&inherits(times.yData,c("POSIXt", "POSIXct", "POSIXlt", "Date"))){warning("'times.yData' has been provided in a Date format, but will be converted to a numeric vector");times.yData<-lubridate::decimal_date(times.yData)}
  if(!is.list(times.xData)&inherits(times.xData,c("POSIXt", "POSIXct", "POSIXlt", "Date"))){warning("'times.xData' has been provided in a Date format, but will be converted to a numeric vector");times.xData<-lubridate::decimal_date(times.xData)}
  if(!is.list(times.xPred)&inherits(times.xPred,c("POSIXt", "POSIXct", "POSIXlt", "Date"))){warning("'times.xPred' has been provided in a Date format, but will be converted to a numeric vector");times.xPred<-lubridate::decimal_date(times.xPred)}
  if(!is.list(times.yPred)&inherits(times.yPred,c("POSIXt", "POSIXct", "POSIXlt", "Date"))){warning("'times.yPred' has been provided in a Date format, but will be converted to a numeric vector");times.yPred<-lubridate::decimal_date(times.yPred)}
  if(is.list(times.yData)&inherits(times.yData[[1]],c("POSIXt", "POSIXct", "POSIXlt", "Date"))){warning("'times.yData' has been provided in a list in Date format, but will be converted to a list of numeric vectors");for(i in 1:length(times.yData)){times.yData[[i]]<-lubridate::decimal_date(times.yData[[i]])}}
  if(is.list(times.xData)&inherits(times.xData[[1]],c("POSIXt", "POSIXct", "POSIXlt", "Date"))){warning("'times.xData' has been provided in a list in Date format, but will be converted to a list of numeric vectors");for(i in 1:length(times.xData)){times.xData[[i]]<-lubridate::decimal_date(times.xData[[i]])}}
  if(is.list(times.xPred)&inherits(times.xPred[[1]],c("POSIXt", "POSIXct", "POSIXlt", "Date"))){warning("'times.xPred' has been provided in a list in Date format, but will be converted to a list of numeric vectors");for(i in 1:length(times.xPred)){times.xPred[[i]]<-lubridate::decimal_date(times.xPred[[i]])}}
  if(is.list(times.yPred)&inherits(times.yPred[[1]],c("POSIXt", "POSIXct", "POSIXlt", "Date"))){warning("'times.yPred' has been provided in a list in Date format, but will be converted to a list of numeric vectors");for(i in 1:length(times.yPred)){times.yPred[[i]]<-lubridate::decimal_date(times.yPred[[i]])}}
  if(!is.list(xPred)&!is.matrix(xPred)&!is.data.frame(xPred)){xPred<-matrix(xPred,ncol=1)}
  if(is.null(times.yData)){
    warning("'times.yData' has not been provided. Defaulting to equally-spaced time points.")
    if(is.matrix(yData)|is.data.frame(yData)){
      times.yData<-1:nrow(yData)
    }else{
      if(is.list(yData)){
        times.yData<-vector("list",length(yData))
        for(i in 1:length(yData)){
          times.yData[[i]]<-1:length(yData[[i]])
        }
      }
    }
  }
  if(is.null(times.xData)){
    warning("'times.xData' has not been provided. Defaulting to equally-spaced time points.")
    if(is.matrix(xData)|is.data.frame(xData)){
      times.xData<-1:nrow(xData)
    }else{
      if(is.list(xData)){
        times.xData<-vector("list",length(xData))
        for(i in 1:length(xData)){
          times.xData[[i]]<-1:length(xData[[i]])
        }
      }
    }
  }
  if(is.null(times.xPred)){
    warning("'times.xPred' has not been provided. Defaulting to equally-spaced time points.")
    if(is.matrix(xPred)|is.data.frame(xPred)){
      times.xPred<-1:nrow(xPred)
    }else{
      if(is.list(xPred)){
        times.xPred<-vector("list",length(xPred))
        for(i in 1:length(xPred)){
          times.xPred[[i]]<-1:length(xPred[[i]])
        }
      }
    }
  }
  if(is.null(times.yPred)){
    if(is.numeric(times.yData)){
      warning("'times.yPred' has not been provided. Defaulting to predictions at 'times.yData'")
      times.yPred<-times.yData
    }else{
      if(is.list(times.yData)){
        times.yPred<-times.yData[[1]]
        warning("'times.yPred' has not been provided. Defaulting to predictions at 'times.yData[[1]]'")
      }
    }
  }
  if(is.matrix(yData)|is.data.frame(yData)){
    yData.mat<-yData
    yData<-matrix.to.list(yData)
  }
  if(is.matrix(xData)|is.data.frame(xData)){
    xData.mat<-xData
    xData<-matrix.to.list(xData)
  }
  if(is.matrix(xPred)|is.data.frame(xPred)){
    xPred.mat<-xPred
    xPred<-matrix.to.list(xPred)
  }
  if(!is.null(times.yData)){
    if(is.matrix(times.yData)|is.data.frame(times.yData)){
      times.yData<-matrix.to.list(times.yData)
    }
    if(is.vector(times.yData)){
      if(exists("yData.mat")){
        times.yData<-times.to.list(yData.mat,times.yData)
      }else{
        stop("'yData' has been provided as a list, while 'times.yData' has been provided as a vector")
      }
    }
  }
  if(!is.null(times.xData)){
    if(is.matrix(times.xData)|is.data.frame(times.xData)){
      times.xData<-matrix.to.list(times.xData)
    }
    if(is.vector(times.xData)){
      if(exists("xData.mat")){
        times.xData<-times.to.list(xData.mat,times.xData)
      }else{
        stop("'xData' has been provided as a list, while 'times.xData' has been provided as a vector")
      }
    }
  }
  if(!is.null(times.xPred)){
    if(is.matrix(times.xPred)|is.data.frame(times.xPred)){
      times.xPred<-matrix.to.list(times.xPred)
    }
    if(is.vector(times.xPred)){
      if(exists("xPred.mat")){
        times.xPred<-times.to.list(xPred.mat,times.xPred)
      }else{
        stop("'xPred' has been provided as a list, while 'times.xPred' has been provided as a vector")
      }
    }
  }
  if(is.matrix(times.yPred)|is.data.frame(times.yPred)|is.list(times.yPred)){
    stop("times.yPred must be a numeric vector")
  }
  if(!is.null(times.yData)&!is.null(times.xData)&!is.null(times.xPred)&!is.null(times.yPred)){
    if(basis.type=="bspline"){
      basis<-fda::create.bspline.basis(rangeval = c(min(c(unlist(times.yData),unlist(times.xData),unlist(times.xPred),times.yPred)),
                                               max(c(unlist(times.yData),unlist(times.xData),unlist(times.xPred),times.yPred))),
                                  nbasis = basis.dim)
    }else{
      if(basis.type=="fourier"){
        basis<-fda::create.fourier.basis(rangeval = c(min(c(unlist(times.yData),unlist(times.xData),unlist(times.xPred),times.yPred)),
                                                 max(c(unlist(times.yData),unlist(times.xData),unlist(times.xPred),times.yPred))),
                                    nbasis = basis.dim,
                                    period = period)
      }else{
        stop("'basis.type' must be either 'bspline' or 'fourier'")
      }
    }
    By<-vector("list",length(yData))
    for(i in 1:length(yData)){
      By[[i]]<-fda::eval.basis(times.yData[[i]],basis)
    }
    Bx<-vector("list",length(xData))
    for(i in 1:length(xData)){
      Bx[[i]]<-fda::eval.basis(times.xData[[i]],basis)
    }
    ByPred<-fda::eval.basis(times.yPred,basis)
    BxPred<-vector("list",length(xPred))
    for(i in 1:length(xPred)){
      BxPred[[i]]<-fda::eval.basis(times.xPred[[i]],basis)
    }
  }
  if(is.null(alphaInit)){
    alphaInit<-matrix(0.1,nrow=basis.dim,ncol=length(yData))
  }
  if(is.null(betaInit)){
    betaInit<-matrix(0.1,nrow=basis.dim,ncol=length(yData))
  }
  if(is.null(cInit)){
    cInit<-matrix(0.1,nrow=basis.dim,ncol=length(yData))
  }
  if(is.null(dInit)){
    dInit<-matrix(0.1,nrow=basis.dim,ncol=length(yData))
  }
  if(is.null(muD)){
    muD<-rep(0,basis.dim)
  }
  if(is.null(SigmaD)){
    SigmaD<-100*diag(basis.dim)
  }
  
  coordsPred<-matrix(coordsPred,ncol = 2)
  
  if(is.list(ByPred)){
    ByPred<-matrix(unlist(ByPred),nrow=nrow(ByPred[[1]]),ncol=ncol(ByPred[[1]]))
  }
  m<-ncol(Bx[[1]])
  nData<-nrow(coordsData)
  nPred<-nrow(coordsPred)
  r<-nrow(ByPred)
  output<-vector("list",nChains)
  for(i in 1:nChains){
    output[[i]]<-NSDmodel(nIter,nThin,yData,xData,xPred,coordsData,coordsPred,By,Bx,ByPred,BxPred,phiAlpha,phiBeta,
                           aAlpha,bAlpha,aBeta,bBeta,aY,bY,aC,bC,aX,bX,muD,SigmaD,
                           sigmaAlphaPrecInit,sigmaBetaPrecInit,sigmaYPrecInit,sigmaCPrecInit,alphaInit,betaInit,cInit,sigmaXPrecInit,dInit)
    if(nBurnIn>0){
      output[[i]]<-output[[i]][-(1:(ceiling(nBurnIn/nThin)+1)),]
    }else{
      output[[i]]<-output[[i]][-1,]
    }
    colnames(output[[i]])<-c("sigmaAlphaPrec","sigmaBetaPrec","sigmaYPrec","sigmaCPrec",paste0("alpha[",rep(1:m,nData),",",rep(1:nData,each=m),"]"),paste0("beta[",rep(1:m,nData),",",rep(1:nData,each=m),"]"),paste0("c[",rep(1:m,nData),",",rep(1:nData,each=m),"]"),"sigmaXPrec",paste0("d[",rep(1:m,nData),",",rep(1:nData,each=m),"]"),paste0("alphaPred[",rep(1:m,nPred),",",rep(1:nPred,each=m),"]"),paste0("betaPred[",rep(1:m,nPred),",",rep(1:nPred,each=m),"]"),paste0("dPred[",rep(1:m,nPred),",",rep(1:nPred,each=m),"]"),paste0("cPred[",rep(1:m,nPred),",",rep(1:nPred,each=m),"]"),paste0("yPred[",rep(1:r,nPred),",",rep(1:nPred,each=r),"]"))
    output[[i]]<-coda::mcmc(output[[i]],end = nIter,thin = nThin)
  }
  output<-coda::mcmc.list(output)
  return(output)
}

run.NSDmodelMulti<-function(nIter=100,nBurnIn=10,nChains=2,nThin=10,yData,xData,zData,xPred,zPred,coordsData,coordsPred,
                             times.yData=NULL,
                             times.xData=NULL,
                             times.zData=NULL,
                             times.xPred=NULL,
                             times.zPred=NULL,
                             times.yPred=NULL,
                             basis.type="bspline",
                             basis.dim=5,
                             period=1,
                             phiAlpha=0.1,phiBeta=0.1,phiGamma=0.1,
                             aAlpha=2,bAlpha=1,aBeta=2,bBeta=1,aGamma=2,bGamma=1,aY=2,bY=1,aC=2,bC=1,aX=2,bX=1,aZ=2,bZ=1,muD=NULL,SigmaD=NULL,muE=NULL,SigmaE=NULL,
                             sigmaAlphaPrecInit=1,sigmaBetaPrecInit=1,sigmaGammaPrecInit=1,sigmaYPrecInit=1,sigmaCPrecInit=1,alphaInit=NULL,betaInit=NULL,gammaInit=NULL,cInit=NULL,sigmaXPrecInit=1,sigmaZPrecInit=1,dInit=NULL,eInit=NULL){
  if(!is.list(times.yData)&inherits(times.yData,c("POSIXt", "POSIXct", "POSIXlt", "Date"))){warning("'times.yData' has been provided in a Date format, but will be converted to a numeric vector");times.yData<-lubridate::decimal_date(times.yData)}
  if(!is.list(times.xData)&inherits(times.xData,c("POSIXt", "POSIXct", "POSIXlt", "Date"))){warning("'times.xData' has been provided in a Date format, but will be converted to a numeric vector");times.xData<-lubridate::decimal_date(times.xData)}
  if(!is.list(times.zData)&inherits(times.zData,c("POSIXt", "POSIXct", "POSIXlt", "Date"))){warning("'times.zData' has been provided in a Date format, but will be converted to a numeric vector");times.zData<-lubridate::decimal_date(times.zData)}
  if(!is.list(times.xPred)&inherits(times.xPred,c("POSIXt", "POSIXct", "POSIXlt", "Date"))){warning("'times.xPred' has been provided in a Date format, but will be converted to a numeric vector");times.xPred<-lubridate::decimal_date(times.xPred)}
  if(!is.list(times.zPred)&inherits(times.zPred,c("POSIXt", "POSIXct", "POSIXlt", "Date"))){warning("'times.zPred' has been provided in a Date format, but will be converted to a numeric vector");times.zPred<-lubridate::decimal_date(times.zPred)}
  if(!is.list(times.yPred)&inherits(times.yPred,c("POSIXt", "POSIXct", "POSIXlt", "Date"))){warning("'times.yPred' has been provided in a Date format, but will be converted to a numeric vector");times.yPred<-lubridate::decimal_date(times.yPred)}
  if(is.list(times.yData)&inherits(times.yData[[1]],c("POSIXt", "POSIXct", "POSIXlt", "Date"))){warning("'times.yData' has been provided in a list in Date format, but will be converted to a list of numeric vectors");for(i in 1:length(times.yData)){times.yData[[i]]<-lubridate::decimal_date(times.yData[[i]])}}
  if(is.list(times.xData)&inherits(times.xData[[1]],c("POSIXt", "POSIXct", "POSIXlt", "Date"))){warning("'times.xData' has been provided in a list in Date format, but will be converted to a list of numeric vectors");for(i in 1:length(times.xData)){times.xData[[i]]<-lubridate::decimal_date(times.xData[[i]])}}
  if(is.list(times.zData)&inherits(times.zData[[1]],c("POSIXt", "POSIXct", "POSIXlt", "Date"))){warning("'times.zData' has been provided in a list in Date format, but will be converted to a list of numeric vectors");for(i in 1:length(times.zData)){times.zData[[i]]<-lubridate::decimal_date(times.zData[[i]])}}
  if(is.list(times.xPred)&inherits(times.xPred[[1]],c("POSIXt", "POSIXct", "POSIXlt", "Date"))){warning("'times.xPred' has been provided in a list in Date format, but will be converted to a list of numeric vectors");for(i in 1:length(times.xPred)){times.xPred[[i]]<-lubridate::decimal_date(times.xPred[[i]])}}
  if(is.list(times.zPred)&inherits(times.zPred[[1]],c("POSIXt", "POSIXct", "POSIXlt", "Date"))){warning("'times.zPred' has been provided in a list in Date format, but will be converted to a list of numeric vectors");for(i in 1:length(times.zPred)){times.zPred[[i]]<-lubridate::decimal_date(times.zPred[[i]])}}
  if(is.list(times.yPred)&inherits(times.yPred[[1]],c("POSIXt", "POSIXct", "POSIXlt", "Date"))){warning("'times.yPred' has been provided in a list in Date format, but will be converted to a list of numeric vectors");for(i in 1:length(times.yPred)){times.yPred[[i]]<-lubridate::decimal_date(times.yPred[[i]])}}
  if(!is.list(xPred)&!is.matrix(xPred)&!is.data.frame(xPred)){xPred<-matrix(xPred,ncol=1)}
  if(!is.list(zPred)&!is.matrix(zPred)&!is.data.frame(zPred)){zPred<-matrix(zPred,ncol=1)}
  if(is.null(times.yData)){
    warning("'times.yData' has not been provided. Defaulting to equally-spaced time points.")
    if(is.matrix(yData)|is.data.frame(yData)){
      times.yData<-1:nrow(yData)
    }else{
      if(is.list(yData)){
        times.yData<-vector("list",length(yData))
        for(i in 1:length(yData)){
          times.yData[[i]]<-1:length(yData[[i]])
        }
      }
    }
  }
  if(is.null(times.xData)){
    warning("'times.xData' has not been provided. Defaulting to equally-spaced time points.")
    if(is.matrix(xData)|is.data.frame(xData)){
      times.xData<-1:nrow(xData)
    }else{
      if(is.list(xData)){
        times.xData<-vector("list",length(xData))
        for(i in 1:length(xData)){
          times.xData[[i]]<-1:length(xData[[i]])
        }
      }
    }
  }
  if(is.null(times.zData)){
    warning("'times.zData' has not been provided. Defaulting to equally-spaced time points.")
    if(is.matrix(zData)|is.data.frame(zData)){
      times.zData<-1:nrow(zData)
    }else{
      if(is.list(zData)){
        times.zData<-vector("list",length(zData))
        for(i in 1:length(zData)){
          times.zData[[i]]<-1:length(zData[[i]])
        }
      }
    }
  }
  if(is.null(times.xPred)){
    warning("'times.xPred' has not been provided. Defaulting to equally-spaced time points.")
    if(is.matrix(xPred)|is.data.frame(xPred)){
      times.xPred<-1:nrow(xPred)
    }else{
      if(is.list(xPred)){
        times.xPred<-vector("list",length(xPred))
        for(i in 1:length(xPred)){
          times.xPred[[i]]<-1:length(xPred[[i]])
        }
      }
    }
  }
  if(is.null(times.zPred)){
    warning("'times.zPred' has not been provided. Defaulting to equally-spaced time points.")
    if(is.matrix(zPred)|is.data.frame(zPred)){
      times.zPred<-1:nrow(zPred)
    }else{
      if(is.list(zPred)){
        times.zPred<-vector("list",length(zPred))
        for(i in 1:length(zPred)){
          times.zPred[[i]]<-1:length(zPred[[i]])
        }
      }
    }
  }
  if(is.null(times.yPred)){
    if(is.numeric(times.yData)){
      warning("'times.yPred' has not been provided. Defaulting to predictions at 'times.yData'")
      times.yPred<-times.yData
    }else{
      if(is.list(times.yData)){
        times.yPred<-times.yData[[1]]
        warning("'times.yPred' has not been provided. Defaulting to predictions at 'times.yData[[1]]'")
      }
    }
  }
  if(is.matrix(yData)|is.data.frame(yData)){
    yData.mat<-yData
    yData<-matrix.to.list(yData)
  }
  if(is.matrix(xData)|is.data.frame(xData)){
    xData.mat<-xData
    xData<-matrix.to.list(xData)
  }
  if(is.matrix(zData)|is.data.frame(zData)){
    zData.mat<-zData
    zData<-matrix.to.list(zData)
  }
  if(is.matrix(xPred)|is.data.frame(xPred)){
    xPred.mat<-xPred
    xPred<-matrix.to.list(xPred)
  }
  if(is.matrix(zPred)|is.data.frame(zPred)){
    zPred.mat<-zPred
    zPred<-matrix.to.list(zPred)
  }
  if(!is.null(times.yData)){
    if(is.matrix(times.yData)|is.data.frame(times.yData)){
      times.yData<-matrix.to.list(times.yData)
    }
    if(is.numeric(times.yData)){
      if(exists("yData.mat")){
        times.yData<-times.to.list(yData.mat,times.yData)
      }else{
        stop("'yData' has been provided as a list, while 'times.yData' has been provided as a numeric vector")
      }
    }
  }
  if(!is.null(times.xData)){
    if(is.matrix(times.xData)|is.data.frame(times.xData)){
      times.xData<-matrix.to.list(times.xData)
    }
    if(is.numeric(times.xData)){
      if(exists("xData.mat")){
        times.xData<-times.to.list(xData.mat,times.xData)
      }else{
        stop("'xData' has been provided as a list, while 'times.xData' has been provided as a numeric vector")
      }
    }
  }
  if(!is.null(times.zData)){
    if(is.matrix(times.zData)|is.data.frame(times.zData)){
      times.zData<-matrix.to.list(times.zData)
    }
    if(is.numeric(times.zData)){
      if(exists("zData.mat")){
        times.zData<-times.to.list(zData.mat,times.zData)
      }else{
        stop("'zData' has been provided as a list, while 'times.zData' has been provided as a numeric vector")
      }
    }
  }
  if(!is.null(times.xPred)){
    if(is.matrix(times.xPred)|is.data.frame(times.xPred)){
      times.xPred<-matrix.to.list(times.xPred)
    }
    if(is.numeric(times.xPred)){
      if(exists("xPred.mat")){
        times.xPred<-times.to.list(xPred.mat,times.xPred)
      }else{
        stop("'xPred' has been provided as a list, while 'times.xPred' has been provided as a numeric vector")
      }
    }
  }
  if(!is.null(times.zPred)){
    if(is.matrix(times.zPred)|is.data.frame(times.zPred)){
      times.zPred<-matrix.to.list(times.zPred)
    }
    if(is.numeric(times.zPred)){
      if(exists("zPred.mat")){
        times.zPred<-times.to.list(zPred.mat,times.zPred)
      }else{
        stop("'zPred' has been provided as a list, while 'times.zPred' has been provided as a numeric vector")
      }
    }
  }
  if(is.matrix(times.yPred)|is.data.frame(times.yPred)|is.list(times.yPred)){
    stop("times.yPred must be a numeric vector")
  }
  if(!is.null(times.yData)&!is.null(times.xData)&!is.null(times.zData)&!is.null(times.xPred)&!is.null(times.zPred)&!is.null(times.yPred)){
    if(basis.type=="bspline"){
      basis<-fda::create.bspline.basis(rangeval = c(min(c(unlist(times.yData),unlist(times.xData),unlist(times.zData),unlist(times.xPred),unlist(times.zPred),times.yPred)),
                                               max(c(unlist(times.yData),unlist(times.xData),unlist(times.zData),unlist(times.xPred),unlist(times.zPred),times.yPred))),
                                  nbasis = basis.dim)
    }else{
      if(basis.type=="fourier"){
        basis<-fda::create.fourier.basis(rangeval = c(min(c(unlist(times.yData),unlist(times.xData),unlist(times.zData),unlist(times.xPred),unlist(times.zPred),times.yPred)),
                                                 max(c(unlist(times.yData),unlist(times.xData),unlist(times.zData),unlist(times.xPred),unlist(times.zPred),times.yPred))),
                                    nbasis = basis.dim,
                                    period = period)
      }else{
        stop("'basis.type' must be either 'bspline' or 'fourier'")
      }
    }
    By<-vector("list",length(yData))
    for(i in 1:length(yData)){
      By[[i]]<-fda::eval.basis(times.yData[[i]],basis)
    }
    Bx<-vector("list",length(xData))
    for(i in 1:length(xData)){
      Bx[[i]]<-fda::eval.basis(times.xData[[i]],basis)
    }
    Bz<-vector("list",length(zData))
    for(i in 1:length(zData)){
      Bz[[i]]<-fda::eval.basis(times.zData[[i]],basis)
    }
    ByPred<-fda::eval.basis(times.yPred,basis)
    BxPred<-vector("list",length(xPred))
    BzPred<-vector("list",length(zPred))
    for(i in 1:length(xPred)){
      BxPred[[i]]<-fda::eval.basis(times.xPred[[i]],basis)
    }
    for(i in 1:length(zPred)){
      BzPred[[i]]<-fda::eval.basis(times.zPred[[i]],basis)
    }
  }
  if(is.null(alphaInit)){
    alphaInit<-matrix(0.1,nrow=basis.dim,ncol=length(yData))
  }
  if(is.null(betaInit)){
    betaInit<-matrix(0.1,nrow=basis.dim,ncol=length(yData))
  }
  if(is.null(gammaInit)){
    gammaInit<-matrix(0.1,nrow=basis.dim,ncol=length(yData))
  }
  if(is.null(cInit)){
    cInit<-matrix(0.1,nrow=basis.dim,ncol=length(yData))
  }
  if(is.null(dInit)){
    dInit<-matrix(0.1,nrow=basis.dim,ncol=length(yData))
  }
  if(is.null(eInit)){
    eInit<-matrix(0.1,nrow=basis.dim,ncol=length(yData))
  }
  if(is.null(muD)){
    muD<-rep(0,basis.dim)
  }
  if(is.null(SigmaD)){
    SigmaD<-100*diag(basis.dim)
  }
  if(is.null(muE)){
    muE<-rep(0,basis.dim)
  }
  if(is.null(SigmaE)){
    SigmaE<-100*diag(basis.dim)
  }
  
  coordsPred<-matrix(coordsPred,ncol = 2)
  if(is.list(ByPred)){
    ByPred<-matrix(unlist(ByPred),nrow=nrow(ByPred[[1]]),ncol=ncol(ByPred[[1]]))
  }
  m<-ncol(Bx[[1]])
  nData<-nrow(coordsData)
  nPred<-nrow(coordsPred)
  r<-nrow(ByPred)
  output<-vector("list",nChains)
  for(i in 1:nChains){
    output[[i]]<-NSDmodelMulti(nIter,nThin,yData,xData,zData,xPred,zPred,coordsData,coordsPred,By,Bx,Bz,ByPred,BxPred,BzPred,phiAlpha,phiBeta,phiGamma,
                                aAlpha,bAlpha,aBeta,bBeta,aGamma,bGamma,aY,bY,aC,bC,aX,bX,aZ,bZ,muD,SigmaD,muE,SigmaE,
                                sigmaAlphaPrecInit,sigmaBetaPrecInit,sigmaGammaPrecInit,sigmaYPrecInit,sigmaCPrecInit,alphaInit,betaInit,gammaInit,cInit,sigmaXPrecInit,sigmaZPrecInit,dInit,eInit)
    if(nBurnIn>0){
      output[[i]]<-output[[i]][-(1:(ceiling(nBurnIn/nThin)+1)),]
    }else{
      output[[i]]<-output[[i]][-1,]
    }
    colnames(output[[i]])<-c("sigmaAlphaPrec","sigmaBetaPrec","sigmaGammaPrec","sigmaYPrec","sigmaCPrec",paste0("alpha[",rep(1:m,nData),",",rep(1:nData,each=m),"]"),paste0("beta[",rep(1:m,nData),",",rep(1:nData,each=m),"]"),paste0("gamma[",rep(1:m,nData),",",rep(1:nData,each=m),"]"),paste0("c[",rep(1:m,nData),",",rep(1:nData,each=m),"]"),"sigmaXPrec","sigmaZPrec",paste0("d[",rep(1:m,nData),",",rep(1:nData,each=m),"]"),paste0("e[",rep(1:m,nData),",",rep(1:nData,each=m),"]"),paste0("alphaPred[",rep(1:m,nPred),",",rep(1:nPred,each=m),"]"),paste0("betaPred[",rep(1:m,nPred),",",rep(1:nPred,each=m),"]"),paste0("gammaPred[",rep(1:m,nPred),",",rep(1:nPred,each=m),"]"),paste0("dPred[",rep(1:m,nPred),",",rep(1:nPred,each=m),"]"),paste0("ePred[",rep(1:m,nPred),",",rep(1:nPred,each=m),"]"),paste0("cPred[",rep(1:m,nPred),",",rep(1:nPred,each=m),"]"),paste0("yPred[",rep(1:r,nPred),",",rep(1:nPred,each=r),"]"))
    output[[i]]<-coda::mcmc(output[[i]],end = nIter,thin = nThin)
  }
  output<-coda::mcmc.list(output)
  return(output)
}

