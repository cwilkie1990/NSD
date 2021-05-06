matrix.to.list<-function(in.mat){
  out.list<-list(ncol(in.mat))
  for(i in 1:ncol(in.mat)){
    out.list[[i]]<-in.mat[!is.na(in.mat[,i]),i]
  }
  return(out.list)
}

times.to.list<-function(in.mat,in.times){
  out.list<-list(ncol(in.mat))
  for(i in 1:ncol(in.mat)){
    out.list[[i]]<-in.times[!is.na(in.mat[,i])]
  }
  return(out.list)
}

summary_NSD<-function(x){
  summary.x<-summary(x)
  preds.rows<-grep("yPred",rownames(summary.x$quantiles))
  n.times.preds<-as.numeric(strsplit(strsplit(rownames(summary.x$quantiles)[nrow(summary.x$quantiles)],split="[",fixed=TRUE)[[1]][2],",",fixed=TRUE)[[1]][1])
  n.locations.preds<-as.numeric(strsplit(strsplit(strsplit(rownames(summary.x$quantiles)[nrow(summary.x$quantiles)],split="[",fixed=TRUE)[[1]][2],",",fixed=TRUE)[[1]][2],"]",fixed=TRUE)[[1]][[1]])
  pred.mat.x<-matrix(summary.x$quantiles[preds.rows,3],nrow=n.times.preds,ncol=n.locations.preds)
  lwrbnd.mat.x<-matrix(summary.x$quantiles[preds.rows,1],nrow=n.times.preds,ncol=n.locations.preds)
  uprbnd.mat.x<-matrix(summary.x$quantiles[preds.rows,5],nrow=n.times.preds,ncol=n.locations.preds)
  out.list<-list(summary.mcmc.list=summary.x,pred.mat=pred.mat.x,lwrbnd.mat=lwrbnd.mat.x,uprbnd.mat=uprbnd.mat.x)
  return(out.list)
}
