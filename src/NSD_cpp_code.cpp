// A fully Bayesian implementation of nonparametric statistical downscaling (incorporating thinning and 
// using an efficient algorithm for the Cholesky decompositions):

#include <cmath>
#include <RcppArmadillo.h>

//[[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
arma::mat calcDistsArma(arma::mat coordsArma){
  NumericMatrix coords = wrap(coordsArma);
  int nRows = coords.nrow();
  double d;
  NumericMatrix dists(nRows,nRows);
  int i,j;
  
  for(i=0; i<nRows; i++){
    for(j=i+1; j<nRows; j++){
      NumericVector v1 = coords.row(i);
      NumericVector v2 = coords.row(j);
      NumericVector v3 = v1-v2;
      d = sqrt(sum(pow(v3,2)));
      dists(j,i)=d;
      dists(i,j)=d;
    }
  }
  arma::mat distsArma(dists.begin(),dists.nrow(),dists.ncol(),false);
  return distsArma;
}


// [[Rcpp::export]]
arma::mat H0functionArma(double phi, arma::mat dists){
  arma::mat Ha = -phi*dists;
  arma::mat H0 = arma::exp(Ha);
  return H0;
}


// [[Rcpp::export]]
arma::mat choleskyEfficient(arma::mat Amat){
  int n = Amat.n_rows;
  arma::vec z(n);
  for(int l=0; l<n; l++){
    z(l) = arma::min(arma::find(Amat.col(l)));
  }
  int dimA = n*n;
  NumericVector L(dimA);
  NumericVector A = wrap(arma::vectorise(Amat));
  for (int i = 0; i < n; i++){
    for (int j = z(i); j < (i+1); j++){
      double s = 0;
      for (int k = 0; k < j; k++){
        s += L(i * n + k) * L(j * n + k);
      }  
      if(i == j){
        L(i * n + j) = sqrt(A(i * n + i) - s);
      }else{
        L(i * n + j) = (1.0 / L(j * n + j) * (A(i * n + j) - s));
      }
    }
  }
  arma::vec LArma(L.begin(),L.size());
  arma::mat LArmaMat(LArma);
  arma::mat LArmaMat2 = arma::reshape(LArmaMat,n,n);
  return LArmaMat2;
}


// [[Rcpp::export]]
NumericMatrix NSDmodel(int nIter, int nThin, List yData, List xData, List xPred, NumericMatrix coordsData, NumericMatrix coordsPred, List By, List Bx, NumericMatrix ByPred, List BxPred, double phiAlpha, double phiBeta,
                        double aAlpha, double bAlpha, double aBeta, double bBeta, double aY, double bY, double aC, double bC, double aX, double bX, NumericVector muD, NumericMatrix SigmaD,
                        double sigmaAlphaPrecInit, double sigmaBetaPrecInit, double sigmaYPrecInit, double sigmaCPrecInit, NumericMatrix alphaInit, NumericMatrix betaInit, NumericMatrix cInit, double sigmaXPrecInit, NumericMatrix dInit){
  int nToSave = nIter/nThin;
  NumericMatrix Bx1 = Bx[1];
  int m = Bx1.ncol();
  int nData = coordsData.nrow();
  int nPred = coordsPred.nrow();
  // NumericVector q = sapply(By,nrow)
  // int p = Bx.nrow();
  int r = ByPred.nrow();
  
  arma::mat coordsDataArma(coordsData.begin(),coordsData.nrow(),coordsData.ncol(),false);
  arma::mat coordsPredArma(coordsPred.begin(),coordsPred.nrow(),coordsPred.ncol(),false);
  arma::mat coordsAllArma = arma::join_cols(coordsPredArma,coordsDataArma);
  arma::mat distsAllArma = calcDistsArma(coordsAllArma);
  
  arma::mat H0AlphaAll = H0functionArma(phiAlpha,distsAllArma);
  arma::mat H0BetaAll = H0functionArma(phiBeta,distsAllArma);
  arma::mat H0Alpha11 = H0AlphaAll.submat(0,0,(nPred-1),(nPred-1));
  arma::mat H0Alpha12 = H0AlphaAll.submat(0,nPred,(nPred-1),(nData+nPred-1));
  arma::mat H0Alpha21 = H0AlphaAll.submat(nPred,0,(nData+nPred-1),(nPred-1));
  arma::mat H0Alpha22 = H0AlphaAll.submat(nPred,nPred,(nData+nPred-1),(nData+nPred-1));
  arma::mat H0Alpha22Inv = arma::inv_sympd(H0Alpha22);
  arma::mat H0Beta11 = H0BetaAll.submat(0,0,(nPred-1),(nPred-1));
  arma::mat H0Beta12 = H0BetaAll.submat(0,nPred,(nPred-1),(nData+nPred-1));
  arma::mat H0Beta21 = H0BetaAll.submat(nPred,0,(nData+nPred-1),(nPred-1));
  arma::mat H0Beta22 = H0BetaAll.submat(nPred,nPred,(nData+nPred-1),(nData+nPred-1));
  arma::mat H0Beta22Inv = arma::inv_sympd(H0Beta22);
  
  arma::mat SigmaDArma(SigmaD.begin(),SigmaD.nrow(),SigmaD.ncol(),false);
  arma::mat SigmaDArmaInv = arma::inv_sympd(SigmaDArma);
  arma::vec muDArma(muD.begin(),muD.size());
  
  NumericVector q(nData);
  for(int i=0; i<nData; i++){
    NumericMatrix ByI = By[i];
    int qI = ByI.nrow();
    q(i) = qI;
  }
  int sumQ = sum(q);
  
  NumericVector p(nData);
  for(int i=0; i<nData; i++){
    NumericMatrix BxI = Bx[i];
    int pI = BxI.nrow();
    p(i) = pI;
  }
  int sumP = sum(p);
  
  // arma::mat xDataArma(xData.begin(),xData.nrow(),xData.ncol(),false);
  // arma::mat xPredArma(xPred.begin(),xPred.nrow(),xPred.ncol(),false);
  // arma::mat BxArma(Bx.begin(),Bx.nrow(),Bx.ncol(),false);
  arma::mat ByPredArma(ByPred.begin(),ByPred.nrow(),ByPred.ncol(),false);
  
  arma::mat SigmaAlphaPredA = H0Alpha11-H0Alpha12*H0Alpha22Inv*H0Alpha21;
  arma::mat SigmaBetaPredA = H0Beta11-H0Beta12*H0Beta22Inv*H0Beta21;
  arma::mat SigmaAlphaPredB = H0Alpha12*H0Alpha22Inv;
  arma::mat SigmaBetaPredB = H0Beta12*H0Beta22Inv;
  
  arma::mat CholSigmaAlphaPredA = choleskyEfficient(SigmaAlphaPredA);
  arma::mat CholSigmaBetaPredA = choleskyEfficient(SigmaBetaPredA);
  
  int nVar2 = 5+4*m*nData+4*m*nPred+r*nPred;
  
  arma::mat outputArma(nToSave+1,nVar2);
  outputArma(0,0) = sigmaAlphaPrecInit;
  outputArma(0,1) = sigmaBetaPrecInit;
  outputArma(0,2) = sigmaYPrecInit;
  outputArma(0,3) = sigmaCPrecInit;
  arma::mat alphaInitArma(alphaInit.begin(),alphaInit.nrow(),alphaInit.ncol(),false);
  arma::mat betaInitArma(betaInit.begin(),betaInit.nrow(),betaInit.ncol(),false);
  arma::mat cInitArma(cInit.begin(),cInit.nrow(),cInit.ncol(),false);
  arma::mat dInitArma(dInit.begin(),dInit.nrow(),dInit.ncol(),false);
  arma::vec alphaInitArmaVec = arma::vectorise(alphaInitArma);
  arma::vec betaInitArmaVec = arma::vectorise(betaInitArma);
  arma::vec cInitArmaVec = arma::vectorise(cInitArma);
  arma::vec dInitArmaVec = arma::vectorise(dInitArma);
  outputArma(0,arma::span(4,m*nData+3)) = alphaInitArmaVec.t();
  outputArma(0,arma::span(m*nData+4,2*m*nData+3)) = betaInitArmaVec.t();
  outputArma(0,arma::span(2*m*nData+4,3*m*nData+3)) = cInitArmaVec.t();
  outputArma(0,3*m*nData+4) = sigmaXPrecInit;
  outputArma(0,arma::span(3*m*nData+5,4*m*nData+4)) = dInitArmaVec.t();
  
  for(int i=1; i<(nToSave+1); i++){
    arma::mat subchainOutputArma(nThin+1,nVar2);
    subchainOutputArma.row(0) = outputArma.row(i-1);
    for(int j=1; j<(nThin+1); j++){
      // double sigmaAlphaPrecPrev = subchainOutputArma(j-1,0);
      // double sigmaBetaPrecPrev = subchainOutputArma(j-1,1);
      // double sigmaYPrecPrev = subchainOutputArma(j-1,2);
      // double sigmaCPrecPrev = subchainOutputArma(j-1,3);
      arma::mat alphaPrev(m,nData);
      arma::mat betaPrev(m,nData);
      arma::mat cPrev(m,nData);
      arma::mat dPrev(m,nData);
      for(int k=0; k<nData; k++){
        alphaPrev.col(k) = subchainOutputArma(j-1,arma::span(4+k*m,3+(k+1)*m)).t();
        betaPrev.col(k) = subchainOutputArma(j-1,arma::span(4+nData*m+k*m,3+nData*m+(k+1)*m)).t();
        cPrev.col(k) = subchainOutputArma(j-1,arma::span(4+2*nData*m+k*m,3+2*nData*m+(k+1)*m)).t();
        dPrev.col(k) = subchainOutputArma(j-1,arma::span(5+3*nData*m+k*m,4+3*nData*m+(k+1)*m)).t();
      }
      // double sigmaXPrecPrev = subchainOutputArma(j-1,4+3*nData*m);
      
      double aAlphaPrecUpdate = aAlpha+m*nData/2;
      double bAlphaPrecUpdate = bAlpha+0.5*arma::trace(H0Alpha22Inv*alphaPrev.t()*alphaPrev);
      double sigmaAlphaPrecUpdate = R::rgamma(aAlphaPrecUpdate,1/bAlphaPrecUpdate);
      
      double aBetaPrecUpdate = aBeta+m*nData/2;
      double bBetaPrecUpdate = bBeta+0.5*arma::trace(H0Beta22Inv*(betaPrev-arma::ones(m,nData)).t()*(betaPrev-arma::ones(m,nData)));
      double sigmaBetaPrecUpdate = R::rgamma(aBetaPrecUpdate,1/bBetaPrecUpdate);
      
      double aYPrecUpdate = aY+0.5*sumQ;
      NumericVector yMBD(nData);
      for(int k=0; k<nData; k++){
        NumericVector yJ = yData[k];
        arma::vec yJArma(yJ.begin(),yJ.size());
        NumericMatrix ByJ = By[k];
        arma::mat ByJArma(ByJ.begin(),ByJ.nrow(),ByJ.ncol(),false);
        arma::vec cJ = cPrev.col(k);
        double yMBDJ = arma::as_scalar((yJArma-ByJArma*cJ).t()*(yJArma-ByJArma*cJ));
        yMBD(k) = yMBDJ;
      }
      double bYPrecUpdate = bY+0.5*sum(yMBD);
      double sigmaYPrecUpdate = R::rgamma(aYPrecUpdate,1/bYPrecUpdate);
      
      double aCPrecUpdate = aC+0.5*m*nData;
      arma::mat dMABG = cPrev-(alphaPrev+betaPrev%dPrev);
      double bCPrecUpdate = bC+0.5*arma::trace(dMABG.t()*dMABG);
      double sigmaCPrecUpdate = R::rgamma(aCPrecUpdate,1/bCPrecUpdate);
      
      arma::mat SigmaAlphaUpdate = arma::inv_sympd(sigmaAlphaPrecUpdate*H0Alpha22Inv + sigmaCPrecUpdate*arma::eye(nData,nData));
      arma::mat CholSigmaAlphaUpdate = choleskyEfficient(SigmaAlphaUpdate);
      arma::mat AAlphaUpdate(m,nData);
      arma::mat alphaUpdate(m,nData);
      for(int k=0; k<m; k++){
        AAlphaUpdate.row(k) = (sigmaCPrecUpdate*arma::eye(nData,nData)*(cPrev.row(k)-betaPrev.row(k)%dPrev.row(k)).t()).t();
        arma::vec muAlphaUpdateJ = SigmaAlphaUpdate*AAlphaUpdate.row(k).t();
        alphaUpdate.row(k) = muAlphaUpdateJ.t() + arma::randn(1,nData)*CholSigmaAlphaUpdate;
      }
      
      arma::mat betaUpdate(m,nData);
      for(int k=0; k<m; k++){
        arma::vec dPrevJ = dPrev.row(k).t();
        arma::vec alphaUpdateJ = alphaUpdate.row(k).t();
        arma::vec cPrevJ = cPrev.row(k).t();
        arma::mat SigmaBetaUpdateJ = arma::inv_sympd(sigmaBetaPrecUpdate*H0Beta22Inv + sigmaCPrecUpdate*arma::diagmat(dPrevJ%dPrevJ));
        arma::vec ABetaUpdateJ = sigmaBetaPrecUpdate*H0Beta22Inv*arma::ones(nData) + sigmaCPrecUpdate*dPrevJ%(cPrevJ-alphaUpdateJ);
        arma::vec muBetaUpdateJ = SigmaBetaUpdateJ*ABetaUpdateJ;
        betaUpdate.row(k) = muBetaUpdateJ.t() + arma::randn(1,nData)*choleskyEfficient(SigmaBetaUpdateJ);
      }
      
      arma::mat cUpdate(m,nData);
      for(int k=0; k<nData; k++){
        NumericVector yJ = yData[k];
        arma::vec yJArma(yJ.begin(),yJ.size());
        NumericMatrix ByJ = By[k];
        arma::mat ByJArma(ByJ.begin(),ByJ.nrow(),ByJ.ncol(),false);
        arma::vec cJ = cPrev.col(k);
        // int qJ = q(k);
        arma::vec alphaUpdateJ = alphaUpdate.col(k);
        arma::vec betaUpdateJ = betaUpdate.col(k);
        arma::vec dPrevJ = dPrev.col(k);
        arma::mat SigmaCUpdateJ = arma::inv_sympd(sigmaYPrecUpdate*ByJArma.t()*ByJArma + sigmaCPrecUpdate*arma::eye(m,m));
        arma::vec ACUpdateJ = sigmaYPrecUpdate*ByJArma.t()*yJArma + sigmaCPrecUpdate*(alphaUpdateJ+betaUpdateJ%dPrevJ);
        arma::vec muCUpdateJ = SigmaCUpdateJ*ACUpdateJ;
        cUpdate.col(k) = (muCUpdateJ.t() + arma::randn(1,m)*choleskyEfficient(SigmaCUpdateJ)).t();
      }
      
      double aXPrecUpdate = aX+0.5*sumP;
      arma::vec xMCGsq(nData);
      for(int k=0; k<nData; k++){
        NumericVector xJ = xData[k];
        arma::vec xJArma(xJ.begin(),xJ.size());
        NumericMatrix BxJ = Bx[k];
        arma::mat BxJArma(BxJ.begin(),BxJ.nrow(),BxJ.ncol(),false);
        arma::vec dJ = dPrev.col(k);
        arma::vec xMCG = xJArma-BxJArma*dJ;
        xMCGsq(k) = arma::as_scalar(xMCG.t()*xMCG);
      }
      double bXPrecUpdate = bX+0.5*sum(xMCGsq);
      double sigmaXPrecUpdate = R::rgamma(aXPrecUpdate,1/bXPrecUpdate);
      
      arma::mat dUpdate(m,nData);
      for(int k=0; k<nData; k++){
        NumericVector xJ = xData[k];
        arma::vec xJArma(xJ.begin(),xJ.size());
        NumericMatrix BxJ = Bx[k];
        arma::mat BxJArma(BxJ.begin(),BxJ.nrow(),BxJ.ncol(),false);
        arma::vec betaJ = betaUpdate.col(k);
        arma::vec alphaJ = alphaUpdate.col(k);
        arma::vec cJ = cUpdate.col(k);
        arma::mat SigmaDUpdateJ = arma::inv_sympd(SigmaDArmaInv + sigmaXPrecUpdate*BxJArma.t()*BxJArma + sigmaCPrecUpdate*arma::diagmat(betaJ%betaJ));
        arma::vec ADUpdateJ = SigmaDArmaInv*muDArma + sigmaXPrecUpdate*BxJArma.t()*xJArma + sigmaCPrecUpdate*betaJ%(cJ-alphaJ);
        arma::vec muDUpdateJ = SigmaDUpdateJ*ADUpdateJ;
        dUpdate.col(k) = (muDUpdateJ.t() + arma::randn(1,m)*choleskyEfficient(SigmaDUpdateJ)).t();
      }
      
      // double sigmaAlphaPrecUpdate = sigmaAlphaPrecPrev;
      // double sigmaBetaPrecUpdate = sigmaBetaPrecPrev;
      // double sigmaYPrecUpdate = sigmaYPrecPrev;
      // double sigmaCPrecUpdate = sigmaCPrecPrev;
      // arma::mat alphaUpdate = alphaPrev;
      // arma::mat betaUpdate = betaPrev;
      // arma::mat cUpdate = cPrev;
      // double sigmaXPrecUpdate = sigmaXPrecPrev;
      // arma::mat dUpdate = dPrev;
      
      subchainOutputArma(j,0) = sigmaAlphaPrecUpdate;
      subchainOutputArma(j,1) = sigmaBetaPrecUpdate;
      subchainOutputArma(j,2) = sigmaYPrecUpdate;
      subchainOutputArma(j,3) = sigmaCPrecUpdate;
      arma::vec alphaUpdateVec = arma::vectorise(alphaUpdate);
      arma::vec betaUpdateVec = arma::vectorise(betaUpdate);
      arma::vec cUpdateVec = arma::vectorise(cUpdate);
      arma::vec dUpdateVec = arma::vectorise(dUpdate);
      subchainOutputArma(j,arma::span(4,m*nData+3)) = alphaUpdateVec.t();
      subchainOutputArma(j,arma::span(m*nData+4,2*m*nData+3)) = betaUpdateVec.t();
      subchainOutputArma(j,arma::span(2*m*nData+4,3*m*nData+3)) = cUpdateVec.t();
      subchainOutputArma(j,3*m*nData+4) = sigmaXPrecUpdate;
      subchainOutputArma(j,arma::span(3*m*nData+5,4*m*nData+4)) = dUpdateVec.t();
      
      arma::mat alphaPredUpdate(m,nPred);
      arma::mat betaPredUpdate(m,nPred);
      arma::mat dPredUpdate(m,nPred);
      arma::mat cPredUpdate(m,nPred);
      for(int k=0; k<m; k++){
        arma::vec alphaJ = alphaUpdate.row(k).t();
        arma::vec muAlphaPredJ = SigmaAlphaPredB*alphaJ;
        // arma::mat SigmaAlphaPredJ = sigmaAlphaPrecUpdate*SigmaAlphaPredA;
        alphaPredUpdate.row(k) = muAlphaPredJ.t() + arma::randn(1,nPred)*sqrt(1/sigmaAlphaPrecUpdate)*CholSigmaAlphaPredA;
        
        arma::vec betaJ = betaUpdate.row(k).t();
        arma::vec muBetaPredJ = arma::ones(nPred) + SigmaBetaPredB*(betaJ-arma::ones(nData));
        // arma::mat SigmaBetaPredJ = sigmaBetaPrecUpdate*SigmaBetaPredA;
        betaPredUpdate.row(k) = muBetaPredJ.t() + arma::randn(1,nPred)*sqrt(1/sigmaBetaPrecUpdate)*CholSigmaBetaPredA;
      }
      for(int k=0; k<nPred; k++){
        NumericVector xPredJ = xPred[k];
        arma::vec xPredJArma(xPredJ.begin(),xPredJ.size());
        NumericMatrix BxPredJ = BxPred[k];
        arma::mat BxPredJArma(BxPredJ.begin(),BxPredJ.nrow(),BxPredJ.ncol(),false);
        arma::mat SigmaDPredJ = arma::inv_sympd(SigmaDArmaInv+sigmaXPrecUpdate*(BxPredJArma.t()*BxPredJArma));
        arma::vec ADPredJ = SigmaDArmaInv*muDArma+sigmaXPrecUpdate*BxPredJArma.t()*xPredJArma;
        arma::vec muDPredJ = SigmaDPredJ*ADPredJ;
        dPredUpdate.col(k) = (muDPredJ.t()+arma::randn(1,m)*choleskyEfficient(SigmaDPredJ)).t();
      }
      for(int k=0; k<m; k++){
        for(int l=0; l<nPred; l++){
          cPredUpdate(k,l) = R::rnorm(alphaPredUpdate(k,l)+betaPredUpdate(k,l)*dPredUpdate(k,l),sqrt(1/sigmaCPrecUpdate));
        }
      }
      arma::vec alphaPredUpdateVec = arma::vectorise(alphaPredUpdate);
      arma::vec betaPredUpdateVec = arma::vectorise(betaPredUpdate);
      arma::vec dPredUpdateVec = arma::vectorise(dPredUpdate);
      arma::vec cPredUpdateVec = arma::vectorise(cPredUpdate);
      subchainOutputArma(j,arma::span(4*m*nData+5,4*m*nData+4+m*nPred)) = alphaPredUpdateVec.t();
      subchainOutputArma(j,arma::span(4*m*nData+5+m*nPred,4*m*nData+4+2*m*nPred)) = betaPredUpdateVec.t();
      subchainOutputArma(j,arma::span(4*m*nData+5+2*m*nPred,4*m*nData+4+3*m*nPred)) = dPredUpdateVec.t();
      subchainOutputArma(j,arma::span(4*m*nData+5+3*m*nPred,4*m*nData+4+4*m*nPred)) = cPredUpdateVec.t();
      
      arma::mat yPredUpdate(r,nPred);
      for(int k=0; k<nPred; k++){
        arma::vec cJ = cPredUpdate.col(k);
        arma::vec muYPredJ = ByPredArma*cJ;
        // arma::mat SigmaYPredJ = (1/sigmaYPrecUpdate)*arma::eye(r,r);
        arma::mat CholSigmaYPredJ = (1/sqrt(sigmaYPrecUpdate))*arma::eye(r,r);
        yPredUpdate.col(k) = (muYPredJ.t() + arma::randn(1,r)*CholSigmaYPredJ).t();
      }
      arma::vec yPredUpdateVec = arma::vectorise(yPredUpdate);
      subchainOutputArma(j,arma::span(4*m*nData+5+4*m*nPred,4*m*nData+4+4*m*nPred+r*nPred)) = yPredUpdateVec.t();
    }
    outputArma.row(i) = subchainOutputArma.row(nThin);
  }
  NumericMatrix output = wrap(outputArma);
  return output;
}



// Model with 3 data sources:

// [[Rcpp::export]]
NumericMatrix NSDmodelMulti(int nIter, int nThin, List yData, List xData, List zData, List xPred, List zPred, NumericMatrix coordsData, NumericMatrix coordsPred, List By, List Bx, List Bz, NumericMatrix ByPred, List BxPred, List BzPred, double phiAlpha, double phiBeta, double phiGamma,
                            double aAlpha, double bAlpha, double aBeta, double bBeta, double aGamma, double bGamma, double aY, double bY, double aC, double bC, double aX, double bX, double aZ, double bZ, NumericVector muD, NumericMatrix SigmaD, NumericVector muE, NumericMatrix SigmaE,
                            double sigmaAlphaPrecInit, double sigmaBetaPrecInit, double sigmaGammaPrecInit, double sigmaYPrecInit, double sigmaCPrecInit, NumericMatrix alphaInit, NumericMatrix betaInit, NumericMatrix gammaInit, NumericMatrix cInit, double sigmaXPrecInit, double sigmaZPrecInit, NumericMatrix dInit, NumericMatrix eInit){
  int nToSave = nIter/nThin;
  NumericMatrix Bx1 = Bx[1];
  int m = Bx1.ncol();
  int nData = coordsData.nrow();
  int nPred = coordsPred.nrow();
  // NumericVector q = sapply(By,nrow)
  // int p = Bx.nrow();
  // int r = Bz.nrow();
  int qPred = ByPred.nrow();
  
  arma::mat coordsDataArma(coordsData.begin(),coordsData.nrow(),coordsData.ncol(),false);
  arma::mat coordsPredArma(coordsPred.begin(),coordsPred.nrow(),coordsPred.ncol(),false);
  arma::mat coordsAllArma = arma::join_cols(coordsPredArma,coordsDataArma);
  arma::mat distsAllArma = calcDistsArma(coordsAllArma);
  
  arma::mat H0AlphaAll = H0functionArma(phiAlpha,distsAllArma);
  arma::mat H0BetaAll = H0functionArma(phiBeta,distsAllArma);
  arma::mat H0GammaAll = H0functionArma(phiGamma,distsAllArma);
  arma::mat H0Alpha11 = H0AlphaAll.submat(0,0,(nPred-1),(nPred-1));
  arma::mat H0Alpha12 = H0AlphaAll.submat(0,nPred,(nPred-1),(nData+nPred-1));
  arma::mat H0Alpha21 = H0AlphaAll.submat(nPred,0,(nData+nPred-1),(nPred-1));
  arma::mat H0Alpha22 = H0AlphaAll.submat(nPred,nPred,(nData+nPred-1),(nData+nPred-1));
  arma::mat H0Alpha22Inv = arma::inv_sympd(H0Alpha22);
  arma::mat H0Beta11 = H0BetaAll.submat(0,0,(nPred-1),(nPred-1));
  arma::mat H0Beta12 = H0BetaAll.submat(0,nPred,(nPred-1),(nData+nPred-1));
  arma::mat H0Beta21 = H0BetaAll.submat(nPred,0,(nData+nPred-1),(nPred-1));
  arma::mat H0Beta22 = H0BetaAll.submat(nPred,nPred,(nData+nPred-1),(nData+nPred-1));
  arma::mat H0Beta22Inv = arma::inv_sympd(H0Beta22);
  arma::mat H0Gamma11 = H0GammaAll.submat(0,0,(nPred-1),(nPred-1));
  arma::mat H0Gamma12 = H0GammaAll.submat(0,nPred,(nPred-1),(nData+nPred-1));
  arma::mat H0Gamma21 = H0GammaAll.submat(nPred,0,(nData+nPred-1),(nPred-1));
  arma::mat H0Gamma22 = H0GammaAll.submat(nPred,nPred,(nData+nPred-1),(nData+nPred-1));
  arma::mat H0Gamma22Inv = arma::inv_sympd(H0Gamma22);
  
  arma::mat SigmaDArma(SigmaD.begin(),SigmaD.nrow(),SigmaD.ncol(),false);
  arma::mat SigmaDArmaInv = arma::inv_sympd(SigmaDArma);
  arma::vec muDArma(muD.begin(),muD.size());
  
  arma::mat SigmaEArma(SigmaE.begin(),SigmaE.nrow(),SigmaE.ncol(),false);
  arma::mat SigmaEArmaInv = arma::inv_sympd(SigmaEArma);
  arma::vec muEArma(muE.begin(),muE.size());
  
  NumericVector q(nData);
  for(int i=0; i<nData; i++){
    NumericMatrix ByI = By[i];
    int qI = ByI.nrow();
    q(i) = qI;
  }
  int sumQ = sum(q);
  
  NumericVector p(nData);
  for(int i=0; i<nData; i++){
    NumericMatrix BxI = Bx[i];
    int pI = BxI.nrow();
    p(i) = pI;
  }
  int sumP = sum(p);
  
  NumericVector r(nData);
  for(int i=0; i<nData; i++){
    NumericMatrix BzI = Bz[i];
    int rI = BzI.nrow();
    r(i) = rI;
  }
  int sumR = sum(r);
  
  // NumericVector qPred(nPred);
  // for(int i=0; i<nPred; i++){
  //   NumericMatrix ByPredI = ByPred[i];
  //   int qPredI = ByPredI.nrow();
  //   qPred(i) = qPredI;
  // }
  // int sumQPred = sum(qPred);
  
  // NumericVector pPred(nPred);
  // for(int i=0; i<nPred; i++){
  //   NumericMatrix BxPredI = BxPred[i];
  //   int pPredI = BxPredI.nrow();
  //   pPred(i) = pPredI;
  // }
  // int sumPPred = sum(pPred);
  // 
  // NumericVector rPred(nPred);
  // for(int i=0; i<nPred; i++){
  //   NumericMatrix BzPredI = BzPred[i];
  //   int rPredI = BzPredI.nrow();
  //   rPred(i) = rPredI;
  // }
  // int sumRPred = sum(rPred);
  
  // arma::mat xDataArma(xData.begin(),xData.nrow(),xData.ncol(),false);
  // arma::mat zDataArma(zData.begin(),zData.nrow(),zData.ncol(),false);
  // arma::mat xPredArma(xPred.begin(),xPred.nrow(),xPred.ncol(),false);
  // arma::mat zPredArma(zPred.begin(),zPred.nrow(),zPred.ncol(),false);
  // arma::mat BxArma(Bx.begin(),Bx.nrow(),Bx.ncol(),false);
  // arma::mat BzArma(Bz.begin(),Bz.nrow(),Bz.ncol(),false);
  arma::mat ByPredArma(ByPred.begin(),ByPred.nrow(),ByPred.ncol(),false);
  
  arma::mat SigmaAlphaPredA = H0Alpha11-H0Alpha12*H0Alpha22Inv*H0Alpha21;
  arma::mat SigmaBetaPredA = H0Beta11-H0Beta12*H0Beta22Inv*H0Beta21;
  arma::mat SigmaGammaPredA = H0Gamma11-H0Gamma12*H0Gamma22Inv*H0Gamma21;
  arma::mat SigmaAlphaPredB = H0Alpha12*H0Alpha22Inv;
  arma::mat SigmaBetaPredB = H0Beta12*H0Beta22Inv;
  arma::mat SigmaGammaPredB = H0Gamma12*H0Gamma22Inv;
  
  arma::mat CholSigmaAlphaPredA = choleskyEfficient(SigmaAlphaPredA);
  arma::mat CholSigmaBetaPredA = choleskyEfficient(SigmaBetaPredA);
  arma::mat CholSigmaGammaPredA = choleskyEfficient(SigmaGammaPredA);
  
  // int nVar2 = 5+4*m*nData+4*m*nPred+r*nPred;
  int nVar2 = 7+6*m*nData+6*m*nPred+qPred*nPred;
  
  arma::mat outputArma(nToSave+1,nVar2);
  outputArma(0,0) = sigmaAlphaPrecInit;
  outputArma(0,1) = sigmaBetaPrecInit;
  outputArma(0,2) = sigmaGammaPrecInit;
  outputArma(0,3) = sigmaYPrecInit;
  outputArma(0,4) = sigmaCPrecInit;
  arma::mat alphaInitArma(alphaInit.begin(),alphaInit.nrow(),alphaInit.ncol(),false);
  arma::mat betaInitArma(betaInit.begin(),betaInit.nrow(),betaInit.ncol(),false);
  arma::mat gammaInitArma(gammaInit.begin(),gammaInit.nrow(),gammaInit.ncol(),false);
  arma::mat cInitArma(cInit.begin(),cInit.nrow(),cInit.ncol(),false);
  arma::mat dInitArma(dInit.begin(),dInit.nrow(),dInit.ncol(),false);
  arma::mat eInitArma(eInit.begin(),eInit.nrow(),eInit.ncol(),false);
  arma::vec alphaInitArmaVec = arma::vectorise(alphaInitArma);
  arma::vec betaInitArmaVec = arma::vectorise(betaInitArma);
  arma::vec gammaInitArmaVec = arma::vectorise(gammaInitArma);
  arma::vec cInitArmaVec = arma::vectorise(cInitArma);
  arma::vec dInitArmaVec = arma::vectorise(dInitArma);
  arma::vec eInitArmaVec = arma::vectorise(eInitArma);
  outputArma(0,arma::span(5,m*nData+4)) = alphaInitArmaVec.t();
  outputArma(0,arma::span(m*nData+5,2*m*nData+4)) = betaInitArmaVec.t();
  outputArma(0,arma::span(2*m*nData+5,3*m*nData+4)) = gammaInitArmaVec.t();
  outputArma(0,arma::span(3*m*nData+5,4*m*nData+4)) = cInitArmaVec.t();
  outputArma(0,4*m*nData+5) = sigmaXPrecInit;
  outputArma(0,4*m*nData+6) = sigmaZPrecInit;
  outputArma(0,arma::span(4*m*nData+7,5*m*nData+6)) = dInitArmaVec.t();
  outputArma(0,arma::span(5*m*nData+7,6*m*nData+6)) = eInitArmaVec.t();
  
  for(int i=1; i<(nToSave+1); i++){
    arma::mat subchainOutputArma(nThin+1,nVar2);
    subchainOutputArma.row(0) = outputArma.row(i-1);
    for(int j=1; j<(nThin+1); j++){
      // double sigmaAlphaPrecPrev = subchainOutputArma(j-1,0);
      // double sigmaBetaPrecPrev = subchainOutputArma(j-1,1);
      // double sigmaYPrecPrev = subchainOutputArma(j-1,2);
      // double sigmaCPrecPrev = subchainOutputArma(j-1,3);
      arma::mat alphaPrev(m,nData);
      arma::mat betaPrev(m,nData);
      arma::mat gammaPrev(m,nData);
      arma::mat cPrev(m,nData);
      arma::mat dPrev(m,nData);
      arma::mat ePrev(m,nData);
      for(int k=0; k<nData; k++){
        alphaPrev.col(k) = subchainOutputArma(j-1,arma::span(5+k*m,4+(k+1)*m)).t();
        betaPrev.col(k) = subchainOutputArma(j-1,arma::span(5+nData*m+k*m,4+nData*m+(k+1)*m)).t();
        gammaPrev.col(k) = subchainOutputArma(j-1,arma::span(5+2*nData*m+k*m,4+2*nData*m+(k+1)*m)).t();
        cPrev.col(k) = subchainOutputArma(j-1,arma::span(5+3*nData*m+k*m,4+3*nData*m+(k+1)*m)).t();
        dPrev.col(k) = subchainOutputArma(j-1,arma::span(7+4*nData*m+k*m,6+4*nData*m+(k+1)*m)).t();
        ePrev.col(k) = subchainOutputArma(j-1,arma::span(7+5*nData*m+k*m,6+5*nData*m+(k+1)*m)).t();
      }
      // double sigmaXPrecPrev = subchainOutputArma(j-1,4+3*nData*m);
      
      double aAlphaPrecUpdate = aAlpha+m*nData/2;
      double bAlphaPrecUpdate = bAlpha+0.5*arma::trace(H0Alpha22Inv*alphaPrev.t()*alphaPrev);
      double sigmaAlphaPrecUpdate = R::rgamma(aAlphaPrecUpdate,1/bAlphaPrecUpdate);
      
      double aBetaPrecUpdate = aBeta+m*nData/2;
      double bBetaPrecUpdate = bBeta+0.5*arma::trace(H0Beta22Inv*(betaPrev-arma::ones(m,nData)).t()*(betaPrev-arma::ones(m,nData)));
      double sigmaBetaPrecUpdate = R::rgamma(aBetaPrecUpdate,1/bBetaPrecUpdate);
      
      double aGammaPrecUpdate = aGamma+m*nData/2;
      double bGammaPrecUpdate = bGamma+0.5*arma::trace(H0Gamma22Inv*(gammaPrev-arma::ones(m,nData)).t()*(gammaPrev-arma::ones(m,nData)));
      double sigmaGammaPrecUpdate = R::rgamma(aGammaPrecUpdate,1/bGammaPrecUpdate);
      
      double aYPrecUpdate = aY+0.5*sumQ;
      NumericVector yMBD(nData);
      for(int k=0; k<nData; k++){
        NumericVector yJ = yData[k];
        arma::vec yJArma(yJ.begin(),yJ.size());
        NumericMatrix ByJ = By[k];
        arma::mat ByJArma(ByJ.begin(),ByJ.nrow(),ByJ.ncol(),false);
        arma::vec cJ = cPrev.col(k);
        double yMBDJ = arma::as_scalar((yJArma-ByJArma*cJ).t()*(yJArma-ByJArma*cJ));
        yMBD(k) = yMBDJ;
      }
      double bYPrecUpdate = bY+0.5*sum(yMBD);
      double sigmaYPrecUpdate = R::rgamma(aYPrecUpdate,1/bYPrecUpdate);
      
      double aCPrecUpdate = aC+0.5*m*nData;
      arma::mat dMABG = cPrev-(alphaPrev+betaPrev%dPrev+gammaPrev%ePrev);
      double bCPrecUpdate = bC+0.5*arma::trace(dMABG.t()*dMABG);
      double sigmaCPrecUpdate = R::rgamma(aCPrecUpdate,1/bCPrecUpdate);
      
      arma::mat SigmaAlphaUpdate = arma::inv_sympd(sigmaAlphaPrecUpdate*H0Alpha22Inv + sigmaCPrecUpdate*arma::eye(nData,nData));
      arma::mat CholSigmaAlphaUpdate = choleskyEfficient(SigmaAlphaUpdate);
      arma::mat AAlphaUpdate(m,nData);
      arma::mat alphaUpdate(m,nData);
      for(int k=0; k<m; k++){
        AAlphaUpdate.row(k) = (sigmaCPrecUpdate*arma::eye(nData,nData)*(cPrev.row(k)-betaPrev.row(k)%dPrev.row(k)-gammaPrev.row(k)%ePrev.row(k)).t()).t();
        arma::vec muAlphaUpdateJ = SigmaAlphaUpdate*AAlphaUpdate.row(k).t();
        alphaUpdate.row(k) = muAlphaUpdateJ.t() + arma::randn(1,nData)*CholSigmaAlphaUpdate;
      }
      
      arma::mat betaUpdate(m,nData);
      for(int k=0; k<m; k++){
        arma::vec dPrevJ = dPrev.row(k).t();
        arma::vec alphaUpdateJ = alphaUpdate.row(k).t();
        arma::vec gammaPrevJ = gammaPrev.row(k).t();
        arma::vec cPrevJ = cPrev.row(k).t();
        arma::vec ePrevJ = ePrev.row(k).t();
        arma::mat SigmaBetaUpdateJ = arma::inv_sympd(sigmaBetaPrecUpdate*H0Beta22Inv + sigmaCPrecUpdate*arma::diagmat(dPrevJ%dPrevJ));
        arma::vec ABetaUpdateJ = sigmaBetaPrecUpdate*H0Beta22Inv*arma::ones(nData) + sigmaCPrecUpdate*dPrevJ%(cPrevJ-alphaUpdateJ-gammaPrevJ%ePrevJ);
        arma::vec muBetaUpdateJ = SigmaBetaUpdateJ*ABetaUpdateJ;
        betaUpdate.row(k) = muBetaUpdateJ.t() + arma::randn(1,nData)*choleskyEfficient(SigmaBetaUpdateJ);
      }
      
      arma::mat gammaUpdate(m,nData);
      for(int k=0; k<m; k++){
        arma::vec cPrevJ = cPrev.row(k).t();
        arma::vec dPrevJ = dPrev.row(k).t();
        arma::vec ePrevJ = ePrev.row(k).t();
        arma::vec alphaUpdateJ = alphaUpdate.row(k).t();
        arma::vec betaUpdateJ = betaUpdate.row(k).t();
        arma::mat SigmaGammaUpdateJ = arma::inv_sympd(sigmaGammaPrecUpdate*H0Gamma22Inv + sigmaCPrecUpdate*arma::diagmat(ePrevJ%ePrevJ));
        arma::vec AGammaUpdateJ = sigmaGammaPrecUpdate*H0Gamma22Inv*arma::ones(nData) + sigmaCPrecUpdate*ePrevJ%(cPrevJ-alphaUpdateJ-betaUpdateJ%dPrevJ);
        arma::vec muGammaUpdateJ = SigmaGammaUpdateJ*AGammaUpdateJ;
        gammaUpdate.row(k) = muGammaUpdateJ.t() + arma::randn(1,nData)*choleskyEfficient(SigmaGammaUpdateJ);
      }
      
      arma::mat cUpdate(m,nData);
      for(int k=0; k<nData; k++){
        NumericVector yJ = yData[k];
        arma::vec yJArma(yJ.begin(),yJ.size());
        NumericMatrix ByJ = By[k];
        arma::mat ByJArma(ByJ.begin(),ByJ.nrow(),ByJ.ncol(),false);
        arma::vec cJ = cPrev.col(k);
        // int qJ = q(k);
        arma::vec alphaUpdateJ = alphaUpdate.col(k);
        arma::vec betaUpdateJ = betaUpdate.col(k);
        arma::vec gammaUpdateJ = gammaUpdate.col(k);
        arma::vec dPrevJ = dPrev.col(k);
        arma::vec ePrevJ = ePrev.col(k);
        arma::mat SigmaCUpdateJ = arma::inv_sympd(sigmaYPrecUpdate*ByJArma.t()*ByJArma + sigmaCPrecUpdate*arma::eye(m,m));
        arma::vec ACUpdateJ = sigmaYPrecUpdate*ByJArma.t()*yJArma + sigmaCPrecUpdate*(alphaUpdateJ+betaUpdateJ%dPrevJ+gammaUpdateJ%ePrevJ);
        arma::vec muCUpdateJ = SigmaCUpdateJ*ACUpdateJ;
        cUpdate.col(k) = (muCUpdateJ.t() + arma::randn(1,m)*choleskyEfficient(SigmaCUpdateJ)).t();
      }
      
      double aXPrecUpdate = aX+0.5*nData*sumP;
      arma::vec xMCGsq(nData);
      for(int k=0; k<nData; k++){
        NumericVector xJ = xData[k];
        arma::vec xJArma(xJ.begin(),xJ.size());
        // arma::vec xDataJ = xDataArma.col(k);
        arma::vec dJ = dPrev.col(k);
        NumericMatrix BxJ = Bx[k];
        arma::mat BxJArma(BxJ.begin(),BxJ.nrow(),BxJ.ncol(),false);
        arma::vec xMCG = xJArma-BxJArma*dJ;
        xMCGsq(k) = arma::as_scalar(xMCG.t()*xMCG);
      }
      double bXPrecUpdate = bX+0.5*sum(xMCGsq);
      double sigmaXPrecUpdate = R::rgamma(aXPrecUpdate,1/bXPrecUpdate);
      
      double aZPrecUpdate = aZ+0.5*nData*sumR;
      arma::vec zMCHsq(nData);
      for(int k=0; k<nData; k++){
        NumericVector zJ = zData[k];
        arma::vec zJArma(zJ.begin(),zJ.size());
        // arma::vec zDataJ = zDataArma.col(k);
        arma::vec eJ = ePrev.col(k);
        NumericMatrix BzJ = Bz[k];
        arma::mat BzJArma(BzJ.begin(),BzJ.nrow(),BzJ.ncol(),false);
        arma::vec zMCH = zJArma-BzJArma*eJ;
        zMCHsq(k) = arma::as_scalar(zMCH.t()*zMCH);
      }
      double bZPrecUpdate = bZ+0.5*sum(zMCHsq);
      double sigmaZPrecUpdate = R::rgamma(aZPrecUpdate,1/bZPrecUpdate);
      
      arma::mat dUpdate(m,nData);
      for(int k=0; k<nData; k++){
        NumericVector xJ = xData[k];
        arma::vec xJArma(xJ.begin(),xJ.size());
        // arma::vec xDataJ = xDataArma.col(k);
        arma::vec betaJ = betaUpdate.col(k);
        arma::vec alphaJ = alphaUpdate.col(k);
        arma::vec gammaJ = gammaUpdate.col(k);
        arma::vec cJ = cUpdate.col(k);
        arma::vec eJ = ePrev.col(k);
        NumericMatrix BxJ = Bx[k];
        arma::mat BxJArma(BxJ.begin(),BxJ.nrow(),BxJ.ncol(),false);
        arma::mat SigmaDUpdateJ = arma::inv_sympd(SigmaDArmaInv + sigmaXPrecUpdate*BxJArma.t()*BxJArma + sigmaCPrecUpdate*arma::diagmat(betaJ%betaJ));
        arma::vec ADUpdateJ = SigmaDArmaInv*muDArma + sigmaXPrecUpdate*BxJArma.t()*xJArma + sigmaCPrecUpdate*betaJ%(cJ-alphaJ-gammaJ%eJ);
        arma::vec muDUpdateJ = SigmaDUpdateJ*ADUpdateJ;
        dUpdate.col(k) = (muDUpdateJ.t() + arma::randn(1,m)*choleskyEfficient(SigmaDUpdateJ)).t();
      }
      
      arma::mat eUpdate(m,nData);
      for(int k=0; k<nData; k++){
        NumericVector zJ = zData[k];
        arma::vec zJArma(zJ.begin(),zJ.size());
        // arma::vec zDataJ = zDataArma.col(k);
        arma::vec betaJ = betaUpdate.col(k);
        arma::vec alphaJ = alphaUpdate.col(k);
        arma::vec gammaJ = gammaUpdate.col(k);
        arma::vec cJ = cUpdate.col(k);
        arma::vec dJ = dUpdate.col(k);
        NumericMatrix BzJ = Bz[k];
        arma::mat BzJArma(BzJ.begin(),BzJ.nrow(),BzJ.ncol(),false);
        arma::mat SigmaEUpdateJ = arma::inv_sympd(SigmaEArmaInv + sigmaZPrecUpdate*BzJArma.t()*BzJArma + sigmaCPrecUpdate*arma::diagmat(gammaJ%gammaJ));
        arma::vec AEUpdateJ = SigmaEArmaInv*muEArma + sigmaZPrecUpdate*BzJArma.t()*zJArma + sigmaCPrecUpdate*gammaJ%(cJ-alphaJ-betaJ%dJ);
        arma::vec muEUpdateJ = SigmaEUpdateJ*AEUpdateJ;
        eUpdate.col(k) = (muEUpdateJ.t() + arma::randn(1,m)*choleskyEfficient(SigmaEUpdateJ)).t();
      }
      
      // double sigmaAlphaPrecUpdate = sigmaAlphaPrecPrev;
      // double sigmaBetaPrecUpdate = sigmaBetaPrecPrev;
      // double sigmaGammaPrecUpdate = sigmaGammaPrecInit;
      // double sigmaYPrecUpdate = sigmaYPrecPrev;
      // double sigmaCPrecUpdate = sigmaCPrecPrev;
      // double sigmaZPrecUpdate = sigmaZPrecInit;
      // arma::mat alphaUpdate = alphaPrev;
      // arma::mat betaUpdate = betaPrev;
      // arma::mat gammaUpdate = gammaInitArma;
      // arma::mat cUpdate = cPrev;
      // double sigmaXPrecUpdate = sigmaXPrecPrev;
      // arma::mat dUpdate = dPrev;
      // arma::mat eUpdate = eInitArma;
      
      subchainOutputArma(j,0) = sigmaAlphaPrecUpdate;
      subchainOutputArma(j,1) = sigmaBetaPrecUpdate;
      subchainOutputArma(j,2) = sigmaGammaPrecUpdate;
      subchainOutputArma(j,3) = sigmaYPrecUpdate;
      subchainOutputArma(j,4) = sigmaCPrecUpdate;
      arma::vec alphaUpdateVec = arma::vectorise(alphaUpdate);
      arma::vec betaUpdateVec = arma::vectorise(betaUpdate);
      arma::vec gammaUpdateVec = arma::vectorise(gammaUpdate);
      arma::vec cUpdateVec = arma::vectorise(cUpdate);
      arma::vec dUpdateVec = arma::vectorise(dUpdate);
      arma::vec eUpdateVec = arma::vectorise(eUpdate);
      subchainOutputArma(j,arma::span(5,m*nData+4)) = alphaUpdateVec.t();
      subchainOutputArma(j,arma::span(m*nData+5,2*m*nData+4)) = betaUpdateVec.t();
      subchainOutputArma(j,arma::span(2*m*nData+5,3*m*nData+4)) = gammaUpdateVec.t();
      subchainOutputArma(j,arma::span(3*m*nData+5,4*m*nData+4)) = cUpdateVec.t();
      subchainOutputArma(j,4*m*nData+5) = sigmaXPrecUpdate;
      subchainOutputArma(j,4*m*nData+6) = sigmaZPrecUpdate;
      subchainOutputArma(j,arma::span(4*m*nData+7,5*m*nData+6)) = dUpdateVec.t();
      subchainOutputArma(j,arma::span(5*m*nData+7,6*m*nData+6)) = eUpdateVec.t();
      
      arma::mat alphaPredUpdate(m,nPred);
      arma::mat betaPredUpdate(m,nPred);
      arma::mat gammaPredUpdate(m,nPred);
      arma::mat dPredUpdate(m,nPred);
      arma::mat ePredUpdate(m,nPred);
      arma::mat cPredUpdate(m,nPred);
      for(int k=0; k<m; k++){
        arma::vec alphaJ = alphaUpdate.row(k).t();
        arma::vec muAlphaPredJ = SigmaAlphaPredB*alphaJ;
        // arma::mat SigmaAlphaPredJ = sigmaAlphaPrecUpdate*SigmaAlphaPredA;
        alphaPredUpdate.row(k) = muAlphaPredJ.t() + arma::randn(1,nPred)*sqrt(1/sigmaAlphaPrecUpdate)*CholSigmaAlphaPredA;
        
        arma::vec betaJ = betaUpdate.row(k).t();
        arma::vec muBetaPredJ = arma::ones(nPred) + SigmaBetaPredB*(betaJ-arma::ones(nData));
        // arma::mat SigmaBetaPredJ = sigmaBetaPrecUpdate*SigmaBetaPredA;
        betaPredUpdate.row(k) = muBetaPredJ.t() + arma::randn(1,nPred)*sqrt(1/sigmaBetaPrecUpdate)*CholSigmaBetaPredA;
        
        arma::vec gammaJ = gammaUpdate.row(k).t();
        arma::vec muGammaPredJ = arma::ones(nPred) + SigmaBetaPredB*(gammaJ-arma::ones(nData));
        gammaPredUpdate.row(k) = muGammaPredJ.t() + arma::randn(1,nPred)*sqrt(1/sigmaGammaPrecUpdate)*CholSigmaGammaPredA;
      }
      for(int k=0; k<nPred; k++){
        NumericVector xPredJ = xPred[k];
        arma::vec xPredJArma(xPredJ.begin(),xPredJ.size());
        // arma::vec xPredJ = xPredArma.col(k);
        NumericMatrix BxPredJ = BxPred[k];
        arma::mat BxPredJArma(BxPredJ.begin(),BxPredJ.nrow(),BxPredJ.ncol(),false);
        arma::mat SigmaDPredJ = arma::inv_sympd(SigmaDArmaInv+sigmaXPrecUpdate*(BxPredJArma.t()*BxPredJArma));
        arma::vec ADPredJ = SigmaDArmaInv*muDArma+sigmaXPrecUpdate*BxPredJArma.t()*xPredJArma;
        arma::vec muDPredJ = SigmaDPredJ*ADPredJ;
        dPredUpdate.col(k) = (muDPredJ.t()+arma::randn(1,m)*choleskyEfficient(SigmaDPredJ)).t();
        
        NumericVector zPredJ = zPred[k];
        arma::vec zPredJArma(zPredJ.begin(),zPredJ.size());
        // arma::vec zPredJ = zPredArma.col(k);
        NumericMatrix BzPredJ = BzPred[k];
        arma::mat BzPredJArma(BzPredJ.begin(),BzPredJ.nrow(),BzPredJ.ncol(),false);
        arma::mat SigmaEPredJ = arma::inv_sympd(SigmaEArmaInv+sigmaZPrecUpdate*(BzPredJArma.t()*BzPredJArma));
        arma::vec AEPredJ = SigmaEArmaInv*muEArma+sigmaZPrecUpdate*BzPredJArma.t()*zPredJArma;
        arma::vec muEPredJ = SigmaEPredJ*AEPredJ;
        ePredUpdate.col(k) = (muEPredJ.t()+arma::randn(1,m)*choleskyEfficient(SigmaEPredJ)).t();
      }
      for(int k=0; k<m; k++){
        for(int l=0; l<nPred; l++){
          cPredUpdate(k,l) = R::rnorm(alphaPredUpdate(k,l)+betaPredUpdate(k,l)*dPredUpdate(k,l)+gammaPredUpdate(k,l)*ePredUpdate(k,l),sqrt(1/sigmaCPrecUpdate));
        }
      }
      arma::vec alphaPredUpdateVec = arma::vectorise(alphaPredUpdate);
      arma::vec betaPredUpdateVec = arma::vectorise(betaPredUpdate);
      arma::vec gammaPredUpdateVec = arma::vectorise(gammaPredUpdate);
      arma::vec dPredUpdateVec = arma::vectorise(dPredUpdate);
      arma::vec ePredUpdateVec = arma::vectorise(ePredUpdate);
      arma::vec cPredUpdateVec = arma::vectorise(cPredUpdate);
      subchainOutputArma(j,arma::span(6*m*nData+7,6*m*nData+6+m*nPred)) = alphaPredUpdateVec.t();
      subchainOutputArma(j,arma::span(6*m*nData+7+m*nPred,6*m*nData+6+2*m*nPred)) = betaPredUpdateVec.t();
      subchainOutputArma(j,arma::span(6*m*nData+7+2*m*nPred,6*m*nData+6+3*m*nPred)) = gammaPredUpdateVec.t();
      subchainOutputArma(j,arma::span(6*m*nData+7+3*m*nPred,6*m*nData+6+4*m*nPred)) = dPredUpdateVec.t();
      subchainOutputArma(j,arma::span(6*m*nData+7+4*m*nPred,6*m*nData+6+5*m*nPred)) = ePredUpdateVec.t();
      subchainOutputArma(j,arma::span(6*m*nData+7+5*m*nPred,6*m*nData+6+6*m*nPred)) = cPredUpdateVec.t();
      
      arma::mat yPredUpdate(qPred,nPred);
      for(int k=0; k<nPred; k++){
        arma::vec cJ = cPredUpdate.col(k);
        arma::vec muYPredJ = ByPredArma*cJ;
        // arma::mat SigmaYPredJ = (1/sigmaYPrecUpdate)*arma::eye(qPred,qPred);
        arma::mat CholSigmaYPredJ = (1/sqrt(sigmaYPrecUpdate))*arma::eye(qPred,qPred);
        yPredUpdate.col(k) = (muYPredJ.t() + arma::randn(1,qPred)*CholSigmaYPredJ).t();
      }
      arma::vec yPredUpdateVec = arma::vectorise(yPredUpdate);
      subchainOutputArma(j,arma::span(6*m*nData+7+6*m*nPred,6*m*nData+6+6*m*nPred+qPred*nPred)) = yPredUpdateVec.t();
    }
    outputArma.row(i) = subchainOutputArma.row(nThin);
  }
  NumericMatrix output = wrap(outputArma);
  return output;
}







