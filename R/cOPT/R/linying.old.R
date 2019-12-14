# Lin-Ying estimator for uni-censored data

linying = function(X,t1,t2) {

	# calculate denominator
	x = apply(cbind(X[,1][,1],X[,2][,1]),1,max);
	e = 1 - X[,1][,2] * X[,2][,2];
	t = max(t1,t2);
	ck = sort(unique(x[e==1]));
	idx = which(ck<t);
	G = 1;
	if( length(idx) > 0 ) {
		for( i in idx ) {
			nk = sum(x>=ck[i]);
			dk = sum(x==ck[i]);
			if( nk > 0 ) { G = G * (1-dk/nk); }
		}
	}

	# calculate numerator
	P = sum(X[,1][,1]>=t1 & X[,2][,1]>=t2)/nrow(X);

	# return
	if( P == 0 | G == 0 ) { return( 0 ); } 

	return( P/G );
}

