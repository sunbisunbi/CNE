# library
library(survival);

##############################################################################
# CenOpt Functions
##############################################################################

# probability in A, P(X in A) when A is finite
opt2d.p.area.sub = function(S,A) {
	P = 0; n = nrow(S);
	a1 = S[,1]; a1[a1<A[1]] = A[1];
	b1 = S[,2]; b1[b1>A[2]] = A[2];
	a2 = S[,3]; a2[a2<A[3]] = A[3];
	b2 = S[,4]; b2[b2>A[4]] = A[4];
	xx = b1-a1; xx[xx<0] = 0;
	yy = b2-a2; yy[yy<0] = 0;
	P = sum(xx*yy*S[,5]);

	return( P );
}

# probability in A, P(X in A)
opt2d.p.area = function(S,A) {
	if( all(is.finite(A)) ) { return( opt2d.p.area.sub(S,A) ); }
	if( is.infinite(A[1]) | is.infinite(A[3]) ) { return( 0 ); }

	# if A is an open space, we recalculate the density by replacing Inf to a large value
	t = S[,1:2]; r1 = range(t[is.finite(t)]);
	t = S[,3:4]; r2 = range(t[is.finite(t)]);
	L = 10*max(c(r1,r2));
	if( is.infinite(A[2]) ) { S[is.infinite(S[,2]),2] = L; A[2] = L; }
	if( is.infinite(A[4]) ) { S[is.infinite(S[,4]),4] = L; A[4] = L; }
	S[,5] = S[,6]/(S[,2]-S[,1])/(S[,4]-S[,3]);
	
	return( opt2d.p.area.sub(S,A) );
}

# probability in A, P(X in A) when A is finite
opt1d.p.area.sub = function(S,A) {
	P = 0; n = nrow(S);
	a1 = S[,1]; a1[a1<A[1]] = A[1];
	b1 = S[,2]; b1[b1>A[2]] = A[2];
	xx = b1-a1; xx[xx<0] = 0;
	P = sum(xx*S[,3]);

	return( P );
}


# probability in A, P(X in A)
opt1d.p.area = function(S,A) {
	if( all(is.finite(A)) ) { return( opt1d.p.area.sub(S,A) ); }
	if( is.infinite(A[1]) ) { return( 0 ); }

	# if A is an open space, we recalculate the density by replacing Inf to a large value
	t = S[,1:2]; r1 = range(t[is.finite(t)]);
	L = 10*max(r1);
	if( is.infinite(A[2]) ) { S[is.infinite(S[,2]),2] = L; A[2] = L; }
	S[,3] = S[,4]/(S[,2]-S[,1]);
	
	return( opt1d.p.area.sub(S,A) );
}

##############################################################################
# MC Functions
##############################################################################

# generate samples
generate_samples = function(S,X_,N=1000,NDIV=11,GEN.BOUND=FALSE) {
  X1 = Surv(X_[,1],X_[,2]); X2 = Surv(X_[,3],X_[,4]);
  X = cbind(as.data.frame(X1),X2);
	if( ncol(S)== 4 ) {
		TALL = generate_samples.1d(S,X,N=N,NDIV=NDIV);
		if( GEN.BOUND ) TALL = generate_samples_boundary.1d(TALL,S);
	} else if( ncol(S) == 6 ) {
		TALL = generate_samples.2d(S,X,N=N,NDIV=NDIV,GEN.Inf=GEN.BOUND);
		#if( GEN.BOUND ) TALL = generate_samples_boundary.2d(TALL,S);
	}

	return( TALL );	
}

generate_samples.2d = function(S,X,N=1000,NDIV=11,GEN.Inf=FALSE,lambda=1) {

	TALL = vector('list',N);
	for( i in 1:N ) {
		TALL[[i]] = matrix(rep(0,nrow(X)*2),ncol=2);
	}

	R = c(range(S[,1]),range(S[,3]));
	dim1 = (R[2]-R[1])/2/(NDIV-1);
	dim2 = (R[4]-R[3])/2/(NDIV-1);

	for( i in 1:nrow(X) ) {
		t1 = X[i,1][,1]; e1 = X[i,1][,2];
		t2 = X[i,2][,1]; e2 = X[i,2][,2];
		if( e1 == 1 & e2 == 1 ) {
			for( k in 1:N ) TALL[[k]][i,] = c(t1,t2);
		} else if( e1 == 0 & e2 == 1 ) {
			bin = c(seq(t1,R[2],length=NDIV),Inf);
			p = rep(0,NDIV);
			for( j in 1:NDIV ) {
				a = c( bin[j], bin[j+1], t2-dim2, t2+dim2 );
				p[j] = opt2d.p.area(S,a);
			}
			if(sum(p)<=0) {p[1]=1;}
			else {p = p/sum(p);}
			cp = c(0,cumsum(p));
			x = t = runif(N); 
			for( j in 1:NDIV ) {
				idx = cp[j]<= x & x < cp[j+1];
				if( sum(idx) == 0 ) next;
				if( !is.infinite(bin[j+1]) ) {
					t[idx] = (bin[j+1]-bin[j])/(cp[j+1]-cp[j])*(x[idx]-cp[j]) + bin[j];
				} else if(GEN.Inf){
				  t[idx] = bin[j] + rexp(1,lambda);
				}	else {
					t[idx] = bin[j];
				}
			}
			for( k in 1:N ) TALL[[k]][i,] = c(t[k],t2);
		} else if( e1 == 1 & e2 == 0 ) {
			bin = c(seq(t2,R[4],length=NDIV),Inf);
			p = rep(0,NDIV);
			for( j in 1:NDIV ) {
				a = c( t1-dim1, t1+dim1, bin[j], bin[j+1] );
				p[j] = opt2d.p.area(S,a);
			}
			if(sum(p)<=0) {p[1]=1;}
			else {p = p/sum(p);}
			cp = c(0,cumsum(p));
			x = t = runif(N); 
			for( j in 1:NDIV ) {
				idx = cp[j]<= x & x < cp[j+1];
				if( sum(idx) == 0 ) next;
				if( !is.infinite(bin[j+1]) ) {
					t[idx] = (bin[j+1]-bin[j])/(cp[j+1]-cp[j])*(x[idx]-cp[j]) + bin[j];
				} else if(GEN.Inf){
				  t[idx] = bin[j] + rexp(1,lambda);
				}	else {
					t[idx] = bin[j];
				}
			}
			for( k in 1:N ) TALL[[k]][i,] = c(t1,t[k]);
		} else if( e1 == 0 & e2 == 0 ) {
			bin1 = c(seq(t1,R[2],length=NDIV),Inf);
			bin2 = c(seq(t2,R[4],length=NDIV),Inf);
			p = matrix(rep(0,NDIV*NDIV),nrow=NDIV);
			for( j in 1:NDIV ) {
				for( jj in 1:NDIV ) {
					a = c( bin1[j],bin1[j+1], bin2[jj],bin2[jj+1]);
					p[j,jj] = opt2d.p.area(S,a);
				}
			}
			if(sum(p)<=0) {p[1,1]=1;}
			else {p = p/sum(p);} 
			p1 = rowSums(p); cp1 = c(0,cumsum(p1)); 
			p2 = colSums(p); cp2 = c(0,cumsum(p2));
			x = t = cbind(runif(N),runif(N));
			for( j1 in 1:NDIV ) {
				for( j2 in 1:NDIV ) {
					idx1 = cp1[j1] <= x[,1] & x[,1] < cp1[j1+1];
					idx2 = cp2[j2] <= x[,2] & x[,2] < cp2[j2+1];
					idx = idx1 & idx2;
					if( sum(idx) == 0 ) next;
					if( !is.infinite(bin1[j1+1]) ) {
						t[idx,1] = (bin1[j1+1]-bin1[j1])/(cp1[j1+1]-cp1[j1])*(x[idx,1]-cp1[j1]) + bin1[j1];
					} else if(GEN.Inf){
					  t[idx,1] = bin1[j1] + rexp(1,lambda);
					} else {
						t[idx,1] = bin1[j1];
					}
					if( !is.infinite(bin2[j2+1]) ) {
						t[idx,2] = (bin2[j2+1]-bin2[j2])/(cp2[j2+1]-cp2[j2])*(x[idx,2]-cp2[j2]) + bin2[j2];
					} else if(GEN.Inf){
					  t[idx,2] = bin2[j2] + rexp(1,lambda);
					} else {
						t[idx,2] = bin2[j2];
					}
				}
			}
			for( k in 1:N ) TALL[[k]][i,] = t[k,];
		}
	}

	return( TALL );
}

generate_samples_boundary.2d = function(TALL,S) {

	R = c(range(S[,1]),range(S[,3]));
	for( kk in 1:length(TALL) ) {
		for( i in 1:nrow(TALL[[kk]]) ) {
			t1 = TALL[[kk]][i,1]; t2 = TALL[[kk]][i,2];	
			if( t1 == R[2] & t2 == R[4] ) {
				TALL[[kk]][i,] = c(NA,NA);
			} else if( t1 == R[2] ) {
				p1 = opt2d.p.area(S,c(0,Inf,t2-0.5,t2+0.5));
				p2 = opt2d.p.area(S,c(t1,Inf,t2-0.5,t2+0.5));
				a = -log(p2/p1)/t1;
				TALL[[kk]][i,1] = rexp(1,a)+t1;
			} else if( t2 == R[4] ) {
				p1 = opt2d.p.area(S,c(t1-0.5,t1+0.5,0,Inf));
				p2 = opt2d.p.area(S,c(t1-0.5,t1+0.5,t2,Inf));
				a = -log(p2/p1)/t2;
				TALL[[kk]][i,2] = rexp(1,a)+t2;
			}
		}
	}

	return( TALL );
}

generate_samples.1d = function(S,X,N=1000,NDIV=11) {

	TALL = vector('list',N);
	for( i in 1:N ) {
		TALL[[i]] = rep(0,nrow(X));
	}

	R = range(S[,1]);
	for( i in 1:nrow(X) ) {
		t1 = X[i][,1]; e1 = X[i][,2];
		if( e1 == 1 ) {
			for( k in 1:N ) TALL[[k]][i] = t1;
		} else if( e1 == 0 ) {
			bin = c(seq(t1,R[2],length=NDIV),Inf);
			p = rep(0,NDIV);
			for( j in 1:NDIV ) {
				a = c( bin[j], bin[j+1] );
				p[j] = opt1d.p.area(S,a);
			}
			p = p/sum(p); cp = c(0,cumsum(p));
			x = t = runif(N); 
			for( j in 1:NDIV ) {
				idx = cp[j]<= x & x < cp[j+1];
				if( sum(idx) == 0 ) next;
				if( !is.infinite(bin[j+1]) ) {
					t[idx] = (bin[j+1]-bin[j])/(cp[j+1]-cp[j])*(x[idx]-cp[j]) + bin[j];
				} else {
					t[idx] = bin[j];
				}
			}
			for( k in 1:N ) TALL[[k]][i] = t[k];
		}
	}

	return( TALL );
}

generate_samples_boundary.1d = function(TALL,S) {
	R = range(S[,1]);
	p = opt1d.p.area(S,c(R[2],Inf));
	a = -log(p)/R[2];
	for( kk in 1:length(TALL) ) {
		idx = TALL[[kk]]==R[2];
		if( sum(idx) == 0 ) next;
		TALL[[kk]][idx] = rexp(sum(idx),a)+R[2];
	}
	return( TALL );
}

