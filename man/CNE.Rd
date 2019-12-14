\name{CNE}
\alias{CNE}
\title{CNE}
\usage{
CNE(survdata)
}
\description{
CNE
}
\examples{
# 1-2    3

T1 = rexp(100,1)
T2 = T1 + rexp(100,1)
T3 = rexp(100,1)
C = matrix(rexp(300),ncol=3)

X1 = apply(cbind(T1,C[,1]),1,min)
X2 = apply(cbind(T2,C[,2]),1,min)
X3 = apply(cbind(T3,C[,3]),1,min)

d1 = T1<=C[,1]
d2 = T1<=C[,2]
d3 = T1<=C[,3]

survdata = cbind(X1,d1,X2,d2,X3,d3)

CNE(survdata)
}
