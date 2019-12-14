library(survival);
library(Matrix);
library(gtools);
library(igraph);
library(glasso);
library(matrixcalc);
library(mnormt);
library(Rcpp);
source("R/CNE.sub.R")

CNE = function(survdata, method = c('CNE', 'dabrowska', 'linying', 'corcoef','naive.X'), l=10000, GEN.BOUND=T){
  method = match.arg(method)
  N.sample = nrow(survdata)
  dimension = ncol(survdata)/2
  if(method=='CNE'){
    Den = opt_(survdata, dimension = 2, repeats = TRUE);
    estimated.T = Estimate_T_MonteCalro(survdata, Den, GEN.BOUND = GEN.BOUND);
    estimated.cov = cov(estimated.T);
    estimated.edge = rho_glasso(estimated.cov, N.sample, l=l);
  } else if(method=='dabrowska'){
    Den = probability_estimation(survdata, method = "dabrowska", repeats=TRUE);
    estimated.T = Estimate_T_MonteCalro(survdata, Den, GEN.BOUND = GEN.BOUND);
    estimated.cov = cov(estimated.T);
    estimated.edge = rho_glasso(estimated.cov, N.sample, l=l);
  } else if(method=='linying'){
    Den = probability_estimation(survdata, method = "linying", repeats=TRUE);
    estimated.T = Estimate_T_MonteCalro(survdata, Den, GEN.BOUND = GEN.BOUND);
    estimated.cov = cov(estimated.T);
    estimated.edge = rho_glasso(estimated.cov, N.sample, l=l);
  } else if(method=='corcoef'){
    estimated.edge = naive.cor(survdata, NUM=20)
  } else if(method=='naive.X'){
    X.cov = cov(survdata[,1:dimension*2-1])
    estimated.edge = rho_glasso(X.cov, N.sample, l=l);
  }
  colnames(estimated.edge) = c('Node1', 'Node2', 'rho')
  return(estimated.edge)
}
