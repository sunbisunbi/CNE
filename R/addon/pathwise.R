# Pathwise Bivariate Survival Estimator
# blocked version of Campbell(1981), "Nonparametric bivariate esitmation with randomly censored data", Biometrika

# estimating a bavariate KM estimator
pathwise.est = function(d,blk1,blk2,method=c('diag','side'),return.info=FALSE) {

	if( missing(method) ) { method='side'; }

	K = length(blk1)-1;
	L = length(blk2)-1;

	D = pathwise.partition_info(d,blk1,blk2);
	S = matrix(rep(-1,(K+1)*(L+1)),nrow=(K+1));
	rownames(S) = blk1;
	colnames(S) = blk2;

	for( k in 1:(K+1) ) {
		for( l in 1:(L+1) ) {
			
			if( k == 1 & l == 1 ) {
				S[k,l] = 1;
			} else if( k == 1 ) {
				r.eff = sum(D$PR[k:K,(l-1):L]) - sum(D$PC2[k:K,l-1])/2; 
				S[k,l] = S[k,l-1] * ( 1 - sum(D$PD2[k:K,l-1])/r.eff );
				if( !is.finite(S[k,l]) ) { S[k,l] = S[k,l-1]; }
			} else if( l == 1 ) {
				r.eff = sum(D$PR[(k-1):K,l:L]) - sum(D$PC1[k-1,l:L])/2; 
				S[k,l] = S[k-1,l] * ( 1 - sum(D$PD1[k-1,l:L])/r.eff );
				if( !is.finite(S[k,l]) ) { S[k,l] = S[k-1,l]; }
			} else {
				r.eff = sum(D$PR[(k-1):K,(l-1):L]) - ( sum(D$PC1[k-1,(l-1):L])+sum(D$PC2[(k-1):K,l-1])-sum(D$PC12[k-1,l-1]) )/2
				s.diag = S[k-1,l-1] * ( 1 - ( sum(D$PD1[k-1,(l-1):L]) + sum(D$PD2[(k-1):K,l-1]) - sum(D$PD12[k-1,l-1]) )/r.eff );

				s.hor = S[k,l-1];
				if( k < K+1 & method=='side' ) {
					r.eff = sum(D$PR[k:K,(l-1):L]) - sum(D$PC2[k:K,l-1])/2; 
					s.hor = S[k,l-1] * ( 1 - sum(D$PD2[k:K,l-1])/r.eff );
				}

				s.ver = S[k-1,l];
				if( l < L+1 & method=='side' ) {
					r.eff = sum(D$PR[(k-1):K,l:L]) - sum(D$PC1[k-1,l:L])/2; 
					s.ver = S[k-1,l] * ( 1 - sum(D$PD1[k-1,l:L])/r.eff );
				}

				if( !is.finite(s.diag) ) { s.diag = S[k-1,l-1]; }
				if( !is.finite(s.hor) ) { s.hor = S[k,l-1]; }
				if( !is.finite(s.ver) ) { s.ver = S[k-1,l]; }

				if( method == 'diag' ) {
					S[k,l] = s.diag;
				} else {
					S[k,l] = mean(c(s.hor,s.ver));
				}
			}
		}
	}

	if( return.info	) {
		return( list(S=S,INFO=D) );
	} else {
		return( S );
	}
}

# get information of each partition
pathwise.partition_info = function(d,blk1,blk2) {
	K = length(blk1)-1; 
	L = length(blk2)-1; 

	F = matrix(rep(0,K*L),nrow=K);
	rownames(F) = rep("",nrow(F));
	colnames(F) = rep("",ncol(F));
	for( k in 1:K ) { rownames(F)[k] = paste("[",blk1[k],",",blk1[k+1],")",sep=""); }
	for( l in 1:L ) { colnames(F)[l] = paste("[",blk2[l],",",blk2[l+1],")",sep=""); }
	
	PR = PD1 = PD2 = PD12 = PC1 = PC2 = PC12 = 0*F;
	for( k in 1:K ) {
		for( l in 1:L ) {
			idx = d[,1][,1]>=blk1[k] & d[,1][,1]<blk1[k+1] & d[,2][,1]>=blk2[l] & d[,2][,1]<blk2[l+1];
			PR[k,l] = sum(idx);
			PD1[k,l] = sum(d[idx,1][,2]==1);
			PD2[k,l] = sum(d[idx,2][,2]==1);
			PD12[k,l] = sum(d[idx,1][,2]==1 & d[idx,2][,2]==1);
			PC1[k,l] = sum(d[idx,1][,2]==0);
			PC2[k,l] = sum(d[idx,2][,2]==0);
			PC12[k,l] = sum(d[idx,1][,2]==0 | d[idx,2][,2]==0);
		}
	}

	return( list(D=D,PR=PR,PD1=PD1,PD2=PD2,PD12=PD12,PC1=PC1,PC2=PC2,PC12=PC12) );
}

# linear estimator within partition
pathwise.p.est.linear.1d = function(x,y1,y2,x1,x2) {
	return( (y2-y1)/(x2-x1)*x + (x2*y1-x1*y2)/(x2-x1) );
}

# survival probability, P(T1>=t1,T2>=t2)
pathwise.p = function(S,t1,t2) {

	if( (is.infinite(t1)&t1>0) | (is.infinite(t2)&t2>0) ) { return( 0 ); }

	blk1 = as.numeric(row.names(S)); blk1 = blk1[is.finite(blk1)];
	blk2 = as.numeric(colnames(S)); blk2 = blk2[is.finite(blk2)];

	k = max(c(0,which(blk1-t1 <= 0)));
	l = max(c(0,which(blk2-t2 <= 0)));

	if( k==0 & l==0 ) { return( 1 ); } 
	else if( k==0 & l<length(blk2) ) { return( pathwise.p.est.linear.1d(t2,S[1,l+1],S[1,l],blk2[l+1],blk2[l]) ); }
	else if( k<length(blk1) & l==0 ) { return( pathwise.p.est.linear.1d(t1,S[k+1,1],S[k,1],blk1[k+1],blk1[k]) ); }
	else if( k==0 & l==length(blk2) ) { return( S[1,l] ); }
	else if( k==length(blk1) & l==0 ) { return( S[k,1] ); }
	else if( k==length(blk1) & l<length(blk2) ) { return( pathwise.p.est.linear.1d(t2,S[k,l+1],S[k,l],blk2[l+1],blk2[l]) ); }
	else if( k<length(blk1) & l==length(blk2) ) { return( pathwise.p.est.linear.1d(t1,S[k+1,l],S[k,l],blk1[k+1],blk1[k]) ); }
	else if( k==length(blk1) & l==length(blk2) ) { return( S[k,l] ); }

	Y = c(S[k,l],S[k+1,l],S[k,l+1],S[k+1,l+1]);
	x1 = blk1[k]; x2 = blk1[k+1];
	y1 = blk2[l]; y2 = blk2[l+1];
	XX = rep(c(x1,x2),2); YY = rep(c(y1,y2),each=2);
	X = cbind(XX,YY,XX*YY,rep(1,4));

	b = solve(X)%*%Y;		
	return( sum( b*c(t1,t2,t1*t2,1) ) );	
}

# survival probability, P(X in A)
pathwise.p.area = function(S,A,PRECISION=1E10) {
	p1 = pathwise.p(S,A[1],A[3]);
	p2 = pathwise.p(S,A[1],A[4]);
	p3 = pathwise.p(S,A[2],A[3]);
	p4 = pathwise.p(S,A[2],A[4]);
	p = p1 - p2 - p3 + p4;
	p = round(p*PRECISION)/PRECISION;
	return( p );
}

# convert a survival function to density function
pathwise.density_grid = function(S,xgrid,ygrid) {

	xgrid = c(xgrid,Inf);
	ygrid = c(ygrid,Inf);

	D = matrix(rep(0,(length(xgrid)-1)*(length(ygrid)-1)),nrow=length(xgrid)-1);
	XN = rep("",length(xgrid)-1);
	YN = rep("",length(ygrid)-1);
	for( i in 2:length(xgrid) ) {
		for( j in 2:length(ygrid) ) {
			D[i-1,j-1] = pathwise.p.area(S,c(xgrid[i-1],xgrid[i],ygrid[j-1],ygrid[j]));
		}
	}
	for( i in 2:length(xgrid) ) { XN[i-1] = sprintf("[%.2f,%.2f)",xgrid[i-1],xgrid[i]); }
	for( i in 2:length(ygrid) ) { YN[i-1] = sprintf("[%.2f,%.2f)",ygrid[i-1],ygrid[i]); }
	row.names(D) = XN;
	colnames(D) = YN;
	
	return( D );
}


