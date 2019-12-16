# 3d Kaplan-Meier Estimator

# library
library(survival);

# functions

pathwise3d.coord2idx = function(k,l,m,K,L,M) {
	return( ((k-1)*(K+1)+l-1)*(L+1)+m );
}

pathwise3d.count.r = function(x,t1,t2,t3) {
	return( sum( x[,1][,1]>t1 & x[,2][,1]>t2 & x[,3][,1]>t3 ) );
}

pathwise3d.count.d = function(x,t1=NULL,t2=NULL,t3=NULL,c1=FALSE,c2=FALSE,c3=FALSE) {
	idx1 = idx2 = idx3 = rep(TRUE,nrow(x));
	if( !is.null(t1) ) { idx1 = x[,1][,1]>t1[1] & x[,1][,1]<=t1[2]; }
	if( !is.null(t2) ) { idx2 = x[,2][,1]>t2[1] & x[,2][,1]<=t2[2]; }
	if( !is.null(t3) ) { idx3 = x[,3][,1]>t3[1] & x[,3][,1]<=t3[2]; }
	idx.range = idx1 & idx2 & idx3;

	idx1 = idx2 = idx3 = rep(TRUE,nrow(x));
	if( c1 ) { idx1 = x[,1][,2]==1; }
	if( c2 ) { idx2 = x[,2][,2]==1; }
	if( c3 ) { idx3 = x[,3][,2]==1; }
	idx.nc = idx1 & idx2 & idx3;

	return( sum( idx.range & idx.nc ) );
}

pathwise3d.count.c = function(x,t1=NULL,t2=NULL,t3=NULL,c1=FALSE,c2=FALSE,c3=FALSE) {

	idx1 = idx2 = idx3 = rep(TRUE,nrow(x));
	if( !is.null(t1) ) { idx1 = x[,1][,1]>t1[1] & x[,1][,1]<=t1[2]; }
	if( !is.null(t2) ) { idx2 = x[,2][,1]>t2[1] & x[,2][,1]<=t2[2]; }
	if( !is.null(t3) ) { idx3 = x[,3][,1]>t3[1] & x[,3][,1]<=t3[2]; }
	idx.range = idx1 & idx2 & idx3;

	idx1 = idx2 = idx3 = rep(FALSE,nrow(x));
	if( c1 ) { idx1 = x[,1][,2]==0; }
	if( c2 ) { idx2 = x[,2][,2]==0; }
	if( c3 ) { idx3 = x[,3][,2]==0; }
	idx.c = idx1 | idx2 | idx3;

	return( sum( idx.range & idx.c ) );
}

# estimating KM on the space
pathwise3d.est = function(d,b1=NULL,b2=NULL,b3=NULL,nb=11,return.info=FALSE) {

	if( is.null(b1) ) { b1 = seq(min(d[,1][,1])*0.9,max(d[,1][,1])*1.1,length=nb); }
	if( is.null(b2) ) { b2 = seq(min(d[,2][,1])*0.9,max(d[,2][,1])*1.1,length=nb); }
	if( is.null(b3) ) { b3 = seq(min(d[,3][,1])*0.9,max(d[,3][,1])*1.1,length=nb); }

	K = length(b1)-1;
	L = length(b2)-1;
	M = length(b3)-1;

	S = matrix(rep(-1,(K+1)*(L+1)*(M+1)*4),ncol=4);

	for( k in 1:(K+1) ) {
		for( l in 1:(L+1) ) {
			for( m in 1:(M+1) ) {
				idx = pathwise3d.coord2idx(k,l,m,K,L,M);
				S[idx,1:3] = c(b1[k],b2[l],b3[m]);
				if( k==1 & l==1 & m==1 ) {
					S[idx,4] = 1;
				} else if( k==1 & l==1 ) {
					r.eff = pathwise3d.count.r(d,b1[k],b2[l],b3[m-1]) - pathwise3d.count.c(d,NULL,NULL,b3[(m-1):m],FALSE,FALSE,TRUE)/2;
					S[idx,4] = S[pathwise3d.coord2idx(k,l,m-1,K,L,M),4] * ( 1 - pathwise3d.count.d(d,NULL,NULL,b3[(m-1):m],FALSE,FALSE,TRUE)/r.eff ); 
					if( !is.finite(S[idx,4]) ) { S[idx,4] = S[pathwise3d.coord2idx(k,l,m-1,K,L,M),4]; }
				} else if( k==1 & m==1 ) {
					r.eff = pathwise3d.count.r(d,b1[k],b2[l-1],b3[m]) - pathwise3d.count.c(d,NULL,b2[(l-1):l],NULL,FALSE,TRUE,FALSE)/2;
					S[idx,4] = S[pathwise3d.coord2idx(k,l-1,m,K,L,M),4] * ( 1 - pathwise3d.count.d(d,NULL,b2[(l-1):l],NULL,FALSE,TRUE,FALSE)/r.eff ); 
					if( !is.finite(S[idx,4]) ) { S[idx,4] = S[pathwise3d.coord2idx(k,l-1,m,K,L,M),4]; }
				} else if( l==1 & m==1 ) {
					r.eff = pathwise3d.count.r(d,b1[k-1],b2[l],b3[m]) - pathwise3d.count.c(d,b1[(k-1):k],NULL,NULL,TRUE,FALSE,FALSE)/2;
					S[idx,4] = S[pathwise3d.coord2idx(k-1,l,m,K,L,M),4] * ( 1 - pathwise3d.count.d(d,b1[(k-1):k],NULL,NULL,TRUE,FALSE,FALSE)/r.eff ); 
					if( !is.finite(S[idx,4]) ) { S[idx,4] = S[pathwise3d.coord2idx(k-1,l,m,K,L,M),4]; }
				} else if( k==1 ) {
					s.hor = S[pathwise3d.coord2idx(k,l,m-1,K,L,M),4];
					if( l < L+1 ) {
						r.eff = pathwise3d.count.r(d,b1[k],b2[l],b3[m-1]) - pathwise3d.count.c(d,NULL,c(b2[l],Inf),b3[(m-1):m],FALSE,FALSE,TRUE)/2;
						s.hor = S[pathwise3d.coord2idx(k,l,m-1,K,L,M),4] * ( 1 - pathwise3d.count.d(d,NULL,c(b2[l],Inf),b3[(m-1):m],FALSE,FALSE,TRUE)/r.eff );
					}
					s.ver = S[pathwise3d.coord2idx(k,l-1,m,K,L,M),4];
					if( m < M+1 ) {
						r.eff = pathwise3d.count.r(d,b1[k],b2[l-1],b3[m]) - pathwise3d.count.c(d,NULL,b2[(l-1):l],c(b3[m],Inf),FALSE,TRUE,FALSE)/2;
						s.ver = S[pathwise3d.coord2idx(k,l-1,m,K,L,M),4] * ( 1 - pathwise3d.count.d(d,NULL,b2[(l-1):l],c(b3[m],Inf),FALSE,TRUE,FALSE)/r.eff );
					}
					if( !is.finite(s.hor) ) { s.hor = S[pathwise3d.coord2idx(k,l,m-1,K,L,M),4]; }
					if( !is.finite(s.ver) ) { s.ver = S[pathwise3d.coord2idx(k,l-1,m,K,L,M),4]; }
					S[idx,4] = mean(c(s.hor,s.ver));
				} else if( l==1 ) {
					s.hor = S[pathwise3d.coord2idx(k,l,m-1,K,L,M),4];
					if( k < K+1 ) {
						r.eff = pathwise3d.count.r(d,b1[k],b2[l],b3[m-1]) - pathwise3d.count.c(d,c(b1[k],Inf),NULL,b3[(m-1):m],FALSE,FALSE,TRUE)/2;
						s.hor = S[pathwise3d.coord2idx(k,l,m-1,K,L,M),4] * ( 1 - pathwise3d.count.d(d,c(b1[k],Inf),NULL,b3[(m-1):m],FALSE,FALSE,TRUE)/r.eff );
					}
					s.ver = S[pathwise3d.coord2idx(k-1,l,m,K,L,M),4];
					if( m < M+1 ) {
						r.eff = pathwise3d.count.r(d,b1[k-1],b2[l],b3[m]) - pathwise3d.count.c(d,b1[(k-1):k],NULL,c(b3[m],Inf),TRUE,FALSE,FALSE)/2;
						s.ver = S[pathwise3d.coord2idx(k-1,l,m,K,L,M),4] * ( 1 - pathwise3d.count.d(d,b1[(k-1):k],NULL,c(b3[m],Inf),TRUE,FALSE,FALSE)/r.eff );
					}
					if( !is.finite(s.hor) ) { s.hor = S[pathwise3d.coord2idx(k,l,m-1,K,L,M),4]; }
					if( !is.finite(s.ver) ) { s.ver = S[pathwise3d.coord2idx(k-1,l,m,K,L,M),4]; }
					S[idx,4] = mean(c(s.hor,s.ver));
				} else if( m==1 ) {
					s.hor = S[pathwise3d.coord2idx(k,l-1,m,K,L,M),4];
					if( k < K+1 ) {
						r.eff = pathwise3d.count.r(d,b1[k],b2[l-1],b3[m]) - pathwise3d.count.c(d,c(b1[k],Inf),b2[(l-1):l],NULL,FALSE,TRUE,FALSE)/2;
						s.hor = S[pathwise3d.coord2idx(k,l-1,m,K,L,M),4] * ( 1 - pathwise3d.count.d(d,c(b1[k],Inf),b2[(l-1):l],NULL,FALSE,TRUE,FALSE)/r.eff );
					}
					s.ver = S[pathwise3d.coord2idx(k-1,l,m,K,L,M),4];
					if( l < L+1 ) {
						r.eff = pathwise3d.count.r(d,b1[k-1],b2[l],b3[m]) - pathwise3d.count.c(d,b1[(k-1):k],c(b2[l],Inf),NULL,TRUE,FALSE,FALSE)/2;
						s.ver = S[pathwise3d.coord2idx(k-1,l,m,K,L,M),4] * ( 1 - pathwise3d.count.d(d,b1[(k-1):k],c(b2[l],Inf),NULL,TRUE,FALSE,FALSE)/r.eff );
					}
					if( !is.finite(s.hor) ) { s.hor = S[pathwise3d.coord2idx(k,l-1,m,K,L,M),4]; }
					if( !is.finite(s.ver) ) { s.ver = S[pathwise3d.coord2idx(k-1,l,m,K,L,M),4]; }
					S[idx,4] = mean(c(s.hor,s.ver));
				} else {
					s.side1 = S[pathwise3d.coord2idx(k-1,l,m,K,L,M),4];
					r.eff = pathwise3d.count.r(d,b1[k-1],b2[l],b3[m]);
					c.eff = pathwise3d.count.c(d,b1[(k-1):k],c(b2[l],Inf),c(b3[m],Inf),TRUE,FALSE,FALSE);
					d.eff = pathwise3d.count.d(d,b1[(k-1):k],c(b2[l],Inf),c(b3[m],Inf),TRUE,FALSE,FALSE);
					s.side1 = S[pathwise3d.coord2idx(k-1,l,m,K,L,M),4] * ( 1 - d.eff/(r.eff-c.eff/2) );	
					if( !is.finite(s.side1) ) { s.side1 = S[pathwise3d.coord2idx(k-1,l,m,K,L,M),4]; }

					s.side2 = S[pathwise3d.coord2idx(k,l-1,m,K,L,M),4];
					r.eff = pathwise3d.count.r(d,b1[k],b2[l-1],b3[m]);
					c.eff = pathwise3d.count.c(d,c(b1[k],Inf),b2[(l-1):l],c(b3[m],Inf),FALSE,TRUE,FALSE);
					d.eff = pathwise3d.count.d(d,c(b1[k],Inf),b2[(l-1):l],c(b3[m],Inf),FALSE,TRUE,FALSE);
					s.side2 = S[pathwise3d.coord2idx(k,l-1,m,K,L,M),4] * ( 1 - d.eff/(r.eff-c.eff/2) );	
					if( !is.finite(s.side2) ) { s.side2 = S[pathwise3d.coord2idx(k,l-1,m,K,L,M),4]; }

					s.side3 = S[pathwise3d.coord2idx(k,l,m-1,K,L,M),4];
					r.eff = pathwise3d.count.r(d,b1[k],b2[l],b3[m-1]);
					c.eff = pathwise3d.count.c(d,c(b1[k],Inf),c(b2[l],Inf),c(b3[(m-1):m]),FALSE,FALSE,TRUE);
					d.eff = pathwise3d.count.d(d,c(b1[k],Inf),c(b2[l],Inf),c(b3[(m-1):m]),FALSE,FALSE,TRUE);
					s.side3 = S[pathwise3d.coord2idx(k,l,m-1,K,L,M),4] * ( 1 - d.eff/(r.eff-c.eff/2) );	
					if( !is.finite(s.side3) ) { s.side3 = S[pathwise3d.coord2idx(k,l,m-1,K,L,M),4]; }

					S[idx,4] = mean(c(s.side1,s.side2,s.side3));
				}
			}
		}
	}
			
	return( S );
}

# linear estimator within partition
pathwise3d.p.est.linear.1d = function(x,y1,y2,x1,x2) {
	return( (y2-y1)/(x2-x1)*x + (x2*y1-x1*y2)/(x2-x1) );
}

pathwise3d.p.est.linear.2d = function(t1,t2,y1,y2,y3,y4,x1,x2,x3,x4) {
	Y = c(y1,y2,y3,y4);
	X = rbind(x1,x2,x3,x4);
	X = cbind(X,X[,1]*X[,2],rep(1,4));
	b = solve(X)%*%Y;
	return( sum( b*c(t1,t2,t1*t2,1) ) );
}

# survival probability, P(T1>=t1,T2>=t2,T3>=t3)
pathwise3d.p = function(S,t1,t2,t3) {

	if( (is.infinite(t1)&t1>0) | (is.infinite(t2)&t2>0) | (is.infinite(t3)&t3>0) ) { return( 0 ); }

	blk1 = unique(S[,1]); blk1 = blk1[is.finite(blk1)];
	blk2 = unique(S[,2]); blk2 = blk2[is.finite(blk2)];
	blk3 = unique(S[,3]); blk3 = blk3[is.finite(blk3)];

	K = length(blk1)-1; 
	L = length(blk2)-1;
	M = length(blk3)-1;

	k = max(c(0,which(blk1-t1 <= 0))); 
	l = max(c(0,which(blk2-t2 <= 0))); 
	m = max(c(0,which(blk3-t3 <= 0))); 
	if( t1 <= 0 ) { k = 0; }
	if( t2 <= 0 ) { l = 0; }
	if( t3 <= 0 ) { m = 0; }

	if( k == 0 ) { k = 1; }
	if( l == 0 ) { l = 1; }
	if( m == 0 ) { m = 1; }

	
	if( k==K+1 & l==L+1 & m==M+1 ) { return( S[pathwise3d.coord2idx(k,l,m,K,L,M),4] ); }
	else if( k==K+1 & l==L+1 & m<=M ) { return( pathwise3d.p.est.linear.1d(t3,S[pathwise3d.coord2idx(k,l,m,K,L,M),4],S[pathwise3d.coord2idx(k,l,m+1,K,L,M),4],blk3[m],blk3[m+1]) ); }
	else if( k==K+1 & l<=L & m==M+1 ) { return( pathwise3d.p.est.linear.1d(t2,S[pathwise3d.coord2idx(k,l,m,K,L,M),4],S[pathwise3d.coord2idx(k,l+1,m,K,L,M),4],blk2[l],blk2[l+1]) ); }
	else if( k<=K & l==L+1 & m==M+1 ) { return( pathwise3d.p.est.linear.1d(t1,S[pathwise3d.coord2idx(k,l,m,K,L,M),4],S[pathwise3d.coord2idx(k+1,l,m,K,L,M),4],blk1[k],blk1[k+1]) ); }
	else if( k==K+1 & l<=L & m<=M ) { return( pathwise3d.p.est.linear.2d(t2,t3, S[pathwise3d.coord2idx(k,l,m,K,L,M),4], S[pathwise3d.coord2idx(k,l+1,m,K,L,M),4], S[pathwise3d.coord2idx(k,l,m+1,K,L,M),4], S[pathwise3d.coord2idx(k,l+1,m+1,K,L,M),4], c(blk2[l],blk3[m]), c(blk2[l+1],blk3[m]), c(blk2[l],blk3[m+1]), c(blk2[l+1],blk3[m+1]) ) ); }
	else if( k<=K & l==L+1 & m<=M ) { return( pathwise3d.p.est.linear.2d(t1,t3, S[pathwise3d.coord2idx(k,l,m,K,L,M),4], S[pathwise3d.coord2idx(k+1,l,m,K,L,M),4], S[pathwise3d.coord2idx(k,l,m+1,K,L,M),4], S[pathwise3d.coord2idx(k+1,l,m+1,K,L,M),4], c(blk1[k],blk3[m]), c(blk1[k+1],blk3[m]), c(blk1[k],blk3[m+1]), c(blk1[k+1],blk3[m+1]) ) ); }
	else if( k<=K & l<=L & m==M+1 ) { return( pathwise3d.p.est.linear.2d(t1,t2, S[pathwise3d.coord2idx(k,l,m,K,L,M),4], S[pathwise3d.coord2idx(k+1,l,m,K,L,M),4], S[pathwise3d.coord2idx(k,l+1,m,K,L,M),4], S[pathwise3d.coord2idx(k+1,l+1,m,K,L,M),4], c(blk1[k],blk2[l]), c(blk1[k+1],blk2[l]), c(blk1[k],blk2[l+1]), c(blk1[k+1],blk2[l+1]) ) ); }

	x1 = blk1[k]; x2 = blk1[k+1];
	y1 = blk2[l]; y2 = blk2[l+1];
	z1 = blk3[m]; z2 = blk3[m+1];

	Y = c(  S[pathwise3d.coord2idx(k,l,m,K,L,M),4],
		S[pathwise3d.coord2idx(k,l,m+1,K,L,M),4],
		S[pathwise3d.coord2idx(k,l+1,m,K,L,M),4],
		S[pathwise3d.coord2idx(k,l+1,m+1,K,L,M),4],
		S[pathwise3d.coord2idx(k+1,l,m,K,L,M),4],
		S[pathwise3d.coord2idx(k+1,l,m+1,K,L,M),4],
		S[pathwise3d.coord2idx(k+1,l+1,m,K,L,M),4],
		S[pathwise3d.coord2idx(k+1,l+1,m+1,K,L,M),4] );
	XX = rep(c(x1,x2),each=4);
	YY = rep(rep(c(y1,y2),each=2),2);
	ZZ = rep(c(z1,z2),4);
	X = cbind(XX,YY,ZZ,XX*YY,YY*ZZ,ZZ*XX,XX*YY*ZZ,rep(1,8));

	b = solve(X)%*%Y;

	return( sum( b*c(t1,t2,t3,t1*t2,t2*t3,t3*t1,t1*t2*t3,1) ) );	
}

# survival probability, P(X in A)
pathwise3d.p.area = function(S,A) {

	p1 = pathwise3d.p(S,A[1],A[3],A[5]);

	p21 = pathwise3d.p(S,A[2],A[3],A[5]);
	p22 = pathwise3d.p(S,A[1],A[4],A[5]);
	p23 = pathwise3d.p(S,A[1],A[3],A[6]);

	p31 = pathwise3d.p(S,A[1],A[4],A[6]);	
	p32 = pathwise3d.p(S,A[2],A[3],A[6]);	
	p33 = pathwise3d.p(S,A[2],A[4],A[5]);	

	p4 = pathwise3d.p(S,A[2],A[4],A[6]);	

	p = p1 - (p21+p22+p23) + (p31+p32+p33) - p4;
	p = round(p*1E10)/1E10;
	return( p );
}

# convert a survival function to density function
pathwise3d.density_grid = function(S,xgrid,ygrid,zgrid) {

	xgrid = c(xgrid,Inf);
	ygrid = c(ygrid,Inf);
	zgrid = c(zgrid,Inf);

	D = rep(0,(length(xgrid)-1)*(length(ygrid)-1)*(length(zgrid)-1));
	DN = rep("",(length(xgrid)-1)*(length(ygrid)-1)*(length(zgrid)-1));
	ii = 1;
	for( i in 2:length(xgrid) ) {
		for( j in 2:length(ygrid) ) {
			for( k in 2:length(zgrid) ) {
				D[ii] = pathwise3d.p.area(S,c(xgrid[i-1],xgrid[i],ygrid[j-1],ygrid[j],zgrid[k-1],zgrid[k]));
				DN[ii] = sprintf("[%.2f,%.2f)x[%.2f,%.2f)x[%.2f,%.2f)", xgrid[i-1],xgrid[i], ygrid[j-1],ygrid[j], zgrid[k-1],zgrid[k]); 
				ii = ii+1;
			}
		}
	}
	names(D) = DN;

	return( D );
}



