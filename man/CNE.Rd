\name{CNE}
\alias{CNE}
\title{CNE}
\usage{
CNE(survdata, method = c('CNE', 'dabrowska', 'linying', 'corcoef','naive.X'), l=10000, cv=FALSE, GEN.BOUND=TRUE)
}
\arguments{
  \item{survdata}{Survival Data including censoring indicator}
  \item{methods}{Network estimation methods. Default 'CNE'}
  \item{l}{Length of penalty parameter. Default 10,000}
  \item{cv}{Cross validation for a penalty parameter. Default FALSE}
  \item{GEN.BOUND}{Monte Carlo simulation option. Default TRUE}
}
\value{
  \item{cv.rho}{Optimal penalty parameter by cv (if cv=T)}
  \item{estimated.edge}{Estimated edges}
}
\description{
Network Estimation for the right-censored times to the occurrence of multiple events data.
}
\examples{
# 3 events
# Network: 1-2 3
library(CNE)
set.seed(12345)
Nsample = 100

T1 = rexp(Nsample,1)
T2 = T1 + rexp(Nsample,1)
T3 = rexp(Nsample,1)
C = matrix(rexp(Nsample*3),ncol=3)

X1 = apply(cbind(T1,C[,1]),1,min)
X2 = apply(cbind(T2,C[,2]),1,min)
X3 = apply(cbind(T3,C[,3]),1,min)

d1 = T1<=C[,1]
d2 = T2<=C[,2]
d3 = T3<=C[,3]

survdata = cbind(X1,d1,X2,d2,X3,d3)
head(survdata)
CNE(survdata)

#      Node1 Node2       rho
# [1,]     1     2 0.28777471
# [2,]     1     3 0.01838880
# [3,]     2     3 0.01565494
}
