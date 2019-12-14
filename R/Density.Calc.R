require(mnormt);
library(gtools);



marginal_dist = function(prob){
  x = sort(unique(c(prob[,1], prob[,2])));
  m.dist = cbind(x[1:(length(x)-1)],x[2:length(x)],rep(0,length(x)-1));
  for(i in 1:nrow(m.dist)){
    idx = which(prob[,1]<=m.dist[i,1]&prob[,2]>=m.dist[i,2]);
    tmp.p = 0;
    if(is.infinite(m.dist[i,2])){
      for(j in idx){
        tmp.p = tmp.p + prob[j,6];
      }
    } else{
      for(j in idx){
        tmp.p = tmp.p + (m.dist[i,2]-m.dist[i,1])*(prob[j,6])/((prob[j,2]-prob[j,1]));
      }     
    }
    m.dist[i,3] = tmp.p;
  }
  return(m.dist);
}

dist_mean = function(prob, p.Inf = TRUE, lambda=1){
  m=0;
  if(sum(is.infinite(prob[,2]))==0){
    m = sum(prob[,3]/2*(prob[,2]+prob[,1]));
  }
  else{
    idx.i = which(!is.infinite(prob[,2]));
    idx.n = which(is.infinite(prob[,2]));
    for(i in idx.i){
      m = m+prob[i,3]/2*(prob[i,2]+prob[i,1]);
    }
    if(p.Inf) {
      for(i in idx.n){
        m = m + prob[i,3]*(prob[i,1]+1/lambda);       
      }
    }
  }
  
  return(m);
}

dist_mean_MC = function(prob, p.Inf=TRUE, B=10000, lambda=1, opt.print=FALSE){
  d = c(0,cumsum(prob[,3]))
  s = runif(B)
  ss = rep(0, B);
  for(i in 1:length(s)){
    for(j in 1:(length(d)-1)){
      if(d[j]<=s[i]&&s[i]<d[j+1]){
        if(is.infinite(prob[j,2])){
          ss[i] = prob[j,1]+rexp(1,lambda);
          break;
        }
        ss[i] = prob[j,1]+(prob[j,2]-prob[j,1])/(d[j+1]-d[j])*(s[i]-d[j]);
        break;
      }
    }
  }
  if(opt.print) {cat("A sample in MC : "); cat(s[1]); cat("\t"); cat(ss[1]); cat("\n");}
  return(mean(ss))
}

xy_mean_dist = function(prob, p.Inf = TRUE, lambda=1){
  x.idx = is.infinite(prob[,2]);
  y.idx = is.infinite(prob[,4]);
  prob_ = prob[!x.idx&!y.idx,]

  Exy_ = prob_[,6]/4*(prob_[,2]+prob_[,1])*(prob_[,4]+prob_[,3]);
  Exy = sum(Exy_);
  if(p.Inf){
    if(sum(x.idx&(!y.idx))>0){
      x.i = which(x.idx&(!y.idx));
      for(i in x.i){
        Exy = Exy+prob[i,6]*(prob[i,1]+1/lambda)*((prob[i,4]+prob[i,3])/2);
      }
    }
    if(sum(y.idx&(!x.idx))>0){
      y.i = which(y.idx&(!x.idx));
      for(i in y.i){
        Exy = Exy+prob[i,6]*(prob[i,3]+1/lambda)*((prob[i,2]+prob[i,1])/2);
      }
    }
    if(sum(x.idx&y.idx)>0){
      xy.i = which(x.idx&y.idx);
      for(i in xy.i){
        Exy = Exy+prob[i,6]*(prob[i,1]+1/lambda)*(prob[i,3]+1/lambda);
      }
    }
  }
  return(Exy);
}

xy_mean_MC = function(prob, p.Inf = TRUE, lambda=1, B=10000, opt.print=FALSE){
  d = c(0,cumsum(prob[,6]));
  s = runif(B)
  sx = rep(0, B);
  sy = rep(0, B);
  for(i in 1:length(s)){
    for(j in 1:(length(d)-1)){
      if(d[j]<=s[i]&&s[i]<d[j+1]){
        if(is.infinite(prob[j,2])&&is.infinite(prob[j,4])){
          sx[i] = prob[j,1]+rexp(1,lambda);
          sy[i] = prob[j,3]+rexp(1,lambda);
          break;
        }
        if(is.infinite(prob[j,2])){
          sx[i] = prob[j,1]+rexp(1,lambda);
          sy[i] = runif(1, prob[j,3], prob[j,4]);
          break;
        }
        if(is.infinite(prob[j,4])){
          sx[i] = runif(1, prob[j,1], prob[j,2]);
          sy[i] = prob[j,3]+rexp(1,lambda);
          break;
        }
        sx[i] = runif(1, prob[j,1], prob[j,2]);
        sy[i] = runif(1, prob[j,3], prob[j,4]);
        break;
      }
    }
  }
  if(opt.print) {
    cat("A sample in MC : "); cat(s[1]); cat("\t(");
    cat(sx[1]); cat(" "); cat(sy[1]); cat(")\n");
  }
  RES = mean(sx*sy);
  return(RES);
}

mean_total = function(l.prob, dim, p.Inf = TRUE, lambda=1){
  MEAN = rep(0,dim);
  name.prob = names(l.prob);
  for(i in 1:length(l.prob)){
    c.idx = as.numeric(strsplit(name.prob[i],'-')[[1]]);
    margin_x = marginal_dist(l.prob[[i]]);
    margin_y = marginal_dist(l.prob[[i]][,c(3,4,1,2,5,6)]);
    MEAN[c.idx[1]] = MEAN[c.idx[1]] + dist_mean(margin_x, p.Inf = p.Inf);
    MEAN[c.idx[2]] = MEAN[c.idx[2]] + dist_mean(margin_y, p.Inf = p.Inf);
  }
  return(MEAN/dim);
}

dist_cov = function(prob, Ex, Ey, p.Inf = TRUE, lambda=1, B=10000, Monte.C = FALSE){
  prob_ = prob[-which(is.infinite(prob[,2])),];
  prob_ = prob_[-which(is.infinite(prob[,4])),];

  if(Monte.C){
    Exy = xy_mean_MC(prob, B=B);
    Ex = dist_mean_MC(margin_x, p.Inf = p.Inf, B=B);
    Ey = dist_mean_MC(margin_y, p.Inf = p.Inf, B=B);
    COV = Exy*(B/(B-1)) - Ex*(B/(B-1))*Ey - Ey*(B/(B-1))*Ex + Ex*Ey;
  } else{
    Exy = xy_mean_dist(prob, p.Inf = p.Inf);
    COV = Exy - Ex*Ey;
  }
  
  return(COV);
}

MY_combinations = function(n, r, v = 1:n, set = TRUE, repeats.allowed=FALSE){

  if(mode(n) != "numeric" || length(n) != 1 
     || n < 1 || (n %% 1) != 0) stop("bad value of n") 
  if(mode(r) != "numeric" || length(r) != 1 
     || r < 1 || (r %% 1) != 0) stop("bad value of r") 
  if(!is.atomic(v) || length(v) < n) 
    stop("v is either non-atomic or too short")
  if( (r > n) & repeats.allowed==FALSE)
    stop("r > n and repeats.allowed=FALSE")
  if(set) {
    
    if (length(v) < n) stop("too few different elements")
  }
  v0 <- vector(mode(v), 0)
  ## Inner workhorse
  if(repeats.allowed)
    sub <- function(n, r, v)
    { 
      if(r == 0) v0 else
        if(r == 1) matrix(v, n, 1) else
          if(n == 1) matrix(v, 1, r) else
            rbind( cbind(v[1], Recall(n, r-1, v)),
                   Recall(n-1, r, v[-1]))
    }
  else
    sub <- function(n, r, v)
    { 
      if(r == 0) v0 else
        if(r == 1) matrix(v, n, 1) else
          if(r == n) matrix(v, 1, n) else
            rbind(cbind(v[1], Recall(n-1, r-1, v[-1])),
                  Recall(n-1, r, v[-1]))
    }
  sub(n, r, v[1:n])
}

   
        
        
        
        

        