##############################################################################
# 4d OPT for censored data
#
# by Junhee Seok
##############################################################################

##############################################################################
# Required Libraries
##############################################################################

library(survival);

##############################################################################
# main OPT functions
##############################################################################

# opt4d
opt4d = function(X,A,max.iter=5,min.err=1E-3,max.depth=5,max.n=1,rho=0.5,alpha=0.5,int.part=FALSE,mixed.mode=FALSE,verb=FALSE) {

	fn.seed = sprintf('%d',round(runif(1)*10000));

	# options
	options = "";
	options = sprintf("%s -s %f-%f,%f-%f,%f-%f,%f-%f",options,A[1],A[2],A[3],A[4],A[5],A[6],A[7],A[8]);
	options = sprintf("%s -i %d -e %f",options,max.iter,min.err);
	options = sprintf("%s -d %d -n %f -r %f -a %f",options,max.depth,max.n,rho,alpha);
	if( int.part ) { options = sprintf("%s -I",options); }
	if( mixed.mode ) { options = sprintf("%s -M",options); }
	out.fn = paste("out",fn.seed,sep='.');
	options = sprintf("%s -o %s",options,out.fn);

	# data file
	d.out = cbind(X[,1][,1],X[,1][,2],X[,2][,1],X[,2][,2],X[,3][,1],X[,3][,2],X[,4][,1],X[,4][,2]);
	data.fn = paste("data",fn.seed,sep='.');
	write.table(d.out,file=data.fn,quote=FALSE,row.names=FALSE,col.names=FALSE,sep='\t');

	# command
	command = sprintf("./copt %s %s",options,data.fn);
	system(command);

	# read output
	S = as.matrix(read.table(out.fn));
	system(sprintf('rm -rf *.%s',fn.seed));

	return( S );	
}

# detail
opt4d.detail = function(X,A,max.iter=5,min.err=1E-3,n.grid=11,max.depth=5,max.n=1,rho=0.5,alpha=0.5,int.part=FALSE,mixed.mode=FALSE,verb=FALSE) {
	S.iter = vector('list',1);
	Dp = Sp = NULL; 
	err = 1; n.iter = 1;
	while( n.iter <= max.iter & err > min.err ) {
		S = opt4d.one(X,A,Sp=Sp,max.depth=max.depth,max.n=max.n,rho=rho,alpha=alpha,int.part=int.part,mixed.mode=mixed.mode,verb=verb);
		D = opt4d.density_grid(S,seq(A[1],A[2],length=n.grid),seq(A[3],A[4],length=n.grid),seq(A[5],A[6],length=n.grid),seq(A[7],A[8],length=n.grid));
		if( !is.null(Dp) ) { err = sqrt(mean((D-Dp)^2)); }
		Sp = S; Dp = D;
		S.iter[[n.iter]] = S;
		cat('iter',n.iter,sprintf('%.5f',err),'\n');
		n.iter = n.iter+1;
	}

	return( list(S=S,S.iter=S.iter) );	
}

# opt4d once
opt4d.one = function(X,A,Sp=NULL,max.depth=5,max.n=1,rho=0.5,alpha=0.5,int.part=FALSE,mixed.mode=FALSE,verb=FALSE) {

	fn.seed = sprintf('%d',round(runif(1)*10000));

	# options
	options = "-O";
	options = sprintf("%s -s %f-%f,%f-%f,%f-%f,%f-%f",options,A[1],A[2],A[3],A[4],A[5],A[6],A[7],A[8]);
	options = sprintf("%s -d %d -n %f -r %f -a %f",options,max.depth,max.n,rho,alpha);
	if( !is.null(Sp) ) {
		fn = paste("prior",fn.seed,sep='.');
		write.table(Sp[,1:10],file=fn,quote=FALSE,row.names=FALSE,col.names=FALSE,sep='\t');
		options = sprintf("%s -p %s",options,fn);
	}
	if( int.part ) { options = sprintf("%s -I",options); }
	if( mixed.mode ) { options = sprintf("%s -M",options); }
	out.fn = paste("out",fn.seed,sep='.');
	options = sprintf("%s -o %s",options,out.fn);

	# data file
	d.out = cbind(X[,1][,1],X[,1][,2],X[,2][,1],X[,2][,2],X[,3][,1],X[,3][,2],X[,4][,1],X[,4][,2]);
	data.fn = paste("data",fn.seed,sep='.');
	write.table(d.out,file=data.fn,quote=FALSE,row.names=FALSE,col.names=FALSE,sep='\t');

	# command
	command = sprintf("./copt %s %s",options,data.fn);
	system(command);

	# read output
	S = as.matrix(read.table(out.fn));

	# remove tmp files
	system(sprintf('rm -rf *.%s',fn.seed));

	return( S );
}

##############################################################################
# probability handlers
##############################################################################

# probability in A, P(X in A) when A is finite
opt4d.p.area.sub = function(S,A) {
        P = 0; n = nrow(S);
        a1 = S[,1]; a1[a1<A[1]] = A[1];
        b1 = S[,2]; b1[b1>A[2]] = A[2];
        a2 = S[,3]; a2[a2<A[3]] = A[3];
        b2 = S[,4]; b2[b2>A[4]] = A[4];
        a3 = S[,5]; a3[a3<A[5]] = A[5];
        b3 = S[,6]; b3[b3>A[6]] = A[6];
        a4 = S[,7]; a4[a4<A[7]] = A[7];
        b4 = S[,8]; b4[b4>A[8]] = A[8];
        xx = b1-a1; xx[xx<0] = 0;
        yy = b2-a2; yy[yy<0] = 0;
        zz = b3-a3; zz[zz<0] = 0;
        ww = b4-a4; ww[ww<0] = 0;
        P = sum(xx*yy*zz*ww*S[,9]);

        return( P );
}

# probability in A, P(X in A)
opt4d.p.area = function(S,A) {
        if( all(is.finite(A)) ) { return( opt4d.p.area.sub(S,A) ); }
        if( is.infinite(A[1]) | is.infinite(A[3]) | is.infinite(A[5]) | is.infinite(A[7]) ) { return( 0 ); }

        # if A is an open space, we recalculate the density by replacing Inf to a large value
        t = S[,1:2]; r1 = range(t[is.finite(t)]);
        t = S[,3:4]; r2 = range(t[is.finite(t)]);
        t = S[,5:6]; r3 = range(t[is.finite(t)]);
        t = S[,7:8]; r4 = range(t[is.finite(t)]);
        L = 10*max(c(r1,r2,r3,r4));
        if( is.infinite(A[2]) ) { S[is.infinite(S[,2]),2] = L; A[2] = L; }
        if( is.infinite(A[4]) ) { S[is.infinite(S[,4]),4] = L; A[4] = L; }
        if( is.infinite(A[6]) ) { S[is.infinite(S[,6]),6] = L; A[6] = L; }
        if( is.infinite(A[8]) ) { S[is.infinite(S[,8]),8] = L; A[8] = L; }
        S[,9] = S[,10]/(S[,2]-S[,1])/(S[,4]-S[,3])/(S[,6]-S[,5])/(S[,8]-S[,7]);

        return( opt4d.p.area.sub(S,A) );
}

# survival probablity P(T1>=t1, T2>=t2)
opt4d.p = function(T,t1,t2,t3) {
        return( opt4d.p.area(T,c(t1,Inf,t2,Inf,t3,Inf)) );
}

##############################################################################
# plotting functions
##############################################################################

##############################################################################
# Misc.
##############################################################################

opt4d.density_grid = function(S,xgrid,ygrid,zgrid,wgrid) {

	xgrid = c(xgrid,Inf);
	ygrid = c(ygrid,Inf);
	zgrid = c(zgrid,Inf);
	wgrid = c(wgrid,Inf);

	D = rep(0,(length(xgrid)-1)*(length(ygrid)-1)*(length(zgrid)-1)*(length(wgrid)-1));
	DN = rep("",(length(xgrid)-1)*(length(ygrid)-1)*(length(zgrid)-1)*(length(wgrid)-1));
	ii = 1;
	for( i in 2:length(xgrid) ) {
		for( j in 2:length(ygrid) ) {
			for( k in 2:length(zgrid) ) {
				for( l in 2:length(wgrid) ) {
					D[ii] = opt4d.p.area(S,c(xgrid[i-1],xgrid[i],ygrid[j-1],ygrid[j],zgrid[k-1],zgrid[k],wgrid[l-1],wgrid[l]));
					DN[ii] = sprintf("[%.2f,%.2f)x[%.2f,%.2f)x[%.2f,%.2f)x[%.2f,%.2f)", xgrid[i-1],xgrid[i], ygrid[j-1],ygrid[j], zgrid[k-1],zgrid[k], wgrid[l-1],wgrid[l]); 
					ii = ii+1;
				}
			}
		}
	}
	names(D) = DN;

	return( D );
}

