# library
library(survival);
library(Matrix);
require(gtools);
library(glasso);
library(matrixcalc);
library(mnormt);
library(Rcpp);
source('R/cOPT/R/opt2d.R');
source('R/cOPT/R/dabrowska.R');
source('R/cOPT/R/linying.R');
source('R/Density.Calc.R');
source('R/opt.MonteCarlo.R');
sourceCpp("R/cOPT/copt.cpp")


Estimate_T_MonteCalro = function(X, Prob, Prob.Dimension = 2, GEN.BOUND=FALSE){
  MC.T = matrix(rep(0, nrow(X)*ncol(X)/2), nrow=nrow(X));
  for(i in 1:length(Prob)){
    Message = sprintf("measure T for %s columns", names(Prob)[i]);
    cat(Message, "\r");
    c.idx = as.numeric(strsplit(names(Prob)[i], '-')[[1]]);
    if(c.idx[1]==c.idx[2]){next;}

    T.Gen = generate_samples(S=Prob[[i]],X_=X[,c(c.idx[1]*2-1,c.idx[1]*2,
                                                 c.idx[2]*2-1,c.idx[2]*2)],
                             GEN.BOUND=GEN.BOUND,N=1000);
    M.T = 0*T.Gen[[1]];
    for( i in 1:length(T.Gen) ) {
      M.T = M.T + T.Gen[[i]]/length(T.Gen);
    }
    MC.T[,c.idx] = MC.T[,c.idx] + M.T/(ncol(X)/2-1);
    flush.console();
  }
  return(MC.T);
}


#############################
# d-Dimensional OPT measure #
#############################
opt_ = function(X, dimension=2, opt.path="copt",
                input.name = "input.surv", output.name = "output.surv", extension = ".txt",
                repeats = FALSE){
  Ncolumn = ncol(X)/2;
  Nrow = nrow(X);
  # nCr
  sequences = combinations(n=Ncolumn, r=dimension, repeats.allowed = repeats);
  Pairs = NULL;
  Prob.Density = list();

  # opt for each pair of columns
  for(i in 1:nrow(sequences)){
    Pairs = c(Pairs, paste(sequences[i,], collapse = "-") );
    # Generate input file for opt
    inputfile = paste(input.name, extension, sep="");
    column.idx = NULL;
    for(j in 1:length(sequences[i,]) ){
      column.idx = c(column.idx, sequences[i,j]*2-1, sequences[i,j]*2);
    }
    write.table(X[,column.idx], file = inputfile, sep = '\t', col.names = FALSE, row.names = FALSE);
    # Set output file name
    outputfile = paste(output.name, extension, sep="");

    # run opt
    Message = sprintf("%d dimension OPT for %s columns", dimension, paste(sequences[i,], collapse = "-"));
    cat(Message, "\r");

    Prob.Density[[length(Prob.Density)+1]] = opt_run(opt.path=opt.path, inputfile, outputfile);
    flush.console();
  }
  names(Prob.Density) = Pairs;

  return(Prob.Density);
}

opt_run = function(opt.path="copt", inputfile, outputfile){

  #system(paste(opt.path, "-d 4 -e 0.05 -o", outputfile, inputfile, sep=" ") );
  cOPT(as.matrix(read.table(inputfile, sep='\t', header = FALSE)))
  if(!file.exists(outputfile)){return(0)}
  Scopt = as.matrix(read.table(outputfile, sep = '\t', header = FALSE) );

  return(Scopt);
}

#########################
# Probability estimates #
#########################
probability_estimation = function(X, method = c("dabrowska", "linying"), repeats = FALSE){
  method = match.arg(method)

  Ncolumn = ncol(X)/2;
  Nrow = nrow(X);
  # nCr
  sequences = combinations(n=Ncolumn, r=2, repeats.allowed = repeats);
  Pairs = NULL;
  Prob.Density = list();

  # Probability estimation for each pair of columns
  for(i in 1:nrow(sequences)){
    Pairs = c(Pairs, paste(sequences[i,], collapse = "-") );
    X1 = X[,c(sequences[i,1]*2-1, sequences[i,1]*2)];
    MIN1 = min(X1[,1]); MAX1 = max(X1[,1]);
    X2 = X[,c(sequences[i,2]*2-1, sequences[i,2]*2)];
    MIN2 = min(X2[,1]); MAX2 = max(X2[,1]);
    A1 = seq(from = MIN1, to = MAX1*1.1, by = (MAX1*1.1-MIN1)/20);
    A2 = seq(from = MIN2, to = MAX2*1.1, by = (MAX2*1.1-MIN2)/20);
    X_ = cbind(X1[,1], X1[,2], X2[,1], X2[,2]);
    Message = sprintf("%s method %s columns", method, paste(sequences[i,], collapse = "-"));
    cat(Message, "\r");
    if(method=="dabrowska"){
      Prob.Density[[length(Prob.Density)+1]] = dabrowska_run(X_, A1, A2);
    }
    else if(method=="linying"){
      Prob.Density[[length(Prob.Density)+1]] = linying_run(X_, A1, A2);
    }

    flush.console();
  }
  names(Prob.Density) = Pairs;

  return(Prob.Density);
}

#######################
# dabrowska estimates #
#######################
dabrowska_run = function(X, A1, A2){
  TM = NULL;
  for( i in A1 ) { for( j in A2 )  TM = rbind(TM,c(i,j)); }
  S = dabrowska(X[,1], X[,3], X[,2], X[,4], TM);
  ss = t(rbind(cbind(matrix(S[,3],nrow=length(A1)),0),0));
  tlist1 = c(A1,Inf);
  tlist2 = c(A2,Inf);
  SS = NULL;
  for( ii in 1:(length(tlist1)-1) ) {
    for( jj in 1:(length(tlist2)-1) ) {
      p = ss[ii,jj] - ss[ii+1,jj] - ss[ii,jj+1] + ss[ii+1,jj+1];
      SS = rbind(SS,c(tlist1[ii],tlist1[ii+1],tlist2[jj],tlist2[jj+1],0,p));
    }
  }
  SS[,5] = SS[,6]/(SS[,2]-SS[,1])/(SS[,4]-SS[,3]);
  Sdb = SS

  return(Sdb);
}

#####################
# linying estimates #
#####################
linying_run = function(X, A1, A2){
  X1 = Surv(X[,1],X[,2]);
  X2 = Surv(X[,3],X[,4]);
  X = cbind(as.data.frame(X1),X2);

  tlist1 = c(A1,Inf);
  tlist2 = c(A2,Inf);
  SS = NULL;
  for( ii in 2:(length(tlist1)) ) {
    for( jj in 2:(length(tlist2)) ) {
      p = linying(X,tlist1[ii-1],tlist2[jj-1]) - linying(X,tlist1[ii-1],tlist2[jj]) -
        linying(X,tlist1[ii],tlist2[jj-1]) + linying(X,tlist1[ii],tlist2[jj]);
      if(is.nan(p)){p=0;}
      if(is.infinite(p)){p=1;}
      if(is.infinite(-p)){p=0;}
      SS = rbind(SS,c(tlist1[ii-1],tlist1[ii],tlist2[jj-1],tlist2[jj],0,p));
    }
  }
  SS[,5] = SS[,6]/(SS[,2]-SS[,1])/(SS[,4]-SS[,3]);
  Sly = SS;

  return(Sly);
}


Measure_covariance = function(l.prob, dim, p.Inf=TRUE){
  COV = matrix(rep(0,dim*dim), nrow=dim);
  MEAN = mean_total(l.prob, dim);
  name.prob = names(l.prob);

  for(i in 1:length(l.prob)){
    c.idx = as.numeric(strsplit(name.prob[i],'-')[[1]]);
    COV[c.idx[1], c.idx[2]] = dist_cov(l.prob[[i]], MEAN[c.idx[1]], MEAN[c.idx[2]], p.Inf=p.Inf, Monte.C=Monte.C)
  }

  RES = COV + t(COV);
  diag(RES) = diag(RES)/2;
  return(RES);
}

rho_glasso = function(Cov, N.sample, l=1000){
  library(gtools);
  dim = ncol(Cov);
  Threshold_list = combinations(dim, r=2);
  Threshold_matrix = matrix(rep(0, dim*dim), nrow=dim);

  rholist = seq(0, max(abs(Cov[upper.tri(Cov)])), length.out = l+1)[-1];
  BIC = rholist;
  for(i in 1:l){
    gl = glasso(Cov, rho = rholist[i])
    wi = gl$wi;
    p_off_d = sum(gl$wi!=0 & col(Cov)<row(Cov));
    #BIC[i] = -2*(gl$loglik)+p_off_d*log(N.sample);
    BIC[i] = N.sample*(-log(determinant(wi)$modulus) + sum(diag(wi*Cov)))+p_off_d*log(N.sample);
    idx = Threshold_matrix==0&wi==0;
    Threshold_matrix[idx] = rholist[i];
    cat(sprintf("(%d/%d)",i,l),'\r');
    flush.console();
  }
  best.rho = rholist[which.min(BIC)];
  Thresholds = cbind(Threshold_list, Threshold_matrix[lower.tri(Threshold_matrix)]);
  #return(list(Thresholds=Thresholds, best.lambda=best.rho));
  return(Thresholds);
}


X.T.cov = function(X){
  dim = ncol(X)/2;
  COV = matrix(rep(0,dim*dim), nrow=dim);

  for(i in 1:dim){
    for(j in 1:dim){
      idx1 = X[,2*i]==1;
      idx2 = X[,2*j]==1;
      idx = idx1&idx2;
      COV[i,j] = cov(X[idx,2*i-1], X[idx,2*j-1]);
    }
  }
  return(COV);
}



naive.cor = function(X, NUM=20){
  nc=ncol(X)/2
  res=NULL
  for(i in 1:(nc-1)){
    for(j in (i+1):nc){
      idx = (X[,i*2]==1)&(X[,j*2]==1)
      if(sum(idx)<NUM) {res = rbind(res, c(i,j,0)) }
      else {res = rbind(res, c(i,j,abs(cor(X[idx,i*2-1],X[idx,j*2-1]))))}
    }
  }
  return(res)
}





