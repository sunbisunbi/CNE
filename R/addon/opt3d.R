##############################################################################
# 2d OPT for censored data
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

# opt3d
opt3d = function(X,A,max.iter=5,min.err=1E-3,max.depth=5,max.n=1,rho=0.5,alpha=0.5,int.part=FALSE,mixed.mode=FALSE,verb=FALSE) {

	fn.seed = sprintf('%d',round(runif(1)*10000));

	# options
	options = "";
	options = sprintf("%s -s %f-%f,%f-%f,%f-%f",options,A[1],A[2],A[3],A[4],A[5],A[6]);
	options = sprintf("%s -i %d -e %f",options,max.iter,min.err);
	options = sprintf("%s -d %d -n %f -r %f -a %f",options,max.depth,max.n,rho,alpha);
	if( int.part ) { options = sprintf("%s -I",options); }
	if( mixed.mode ) { options = sprintf("%s -M",options); }
	out.fn = paste("out",fn.seed,sep='.');
	options = sprintf("%s -o %s",options,out.fn);

	# data file
	d.out = cbind(X[,1][,1],X[,1][,2],X[,2][,1],X[,2][,2],X[,3][,1],X[,3][,2]);
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
opt3d.detail = function(X,A,max.iter=5,min.err=1E-3,n.grid=11,max.depth=5,max.n=1,rho=0.5,alpha=0.5,int.part=FALSE,mixed.mode=FALSE,verb=FALSE) {
	S.iter = vector('list',1);
	Dp = Sp = NULL; 
	err = 1; n.iter = 1;
	while( n.iter <= max.iter & err > min.err ) {
		S = opt3d.one(X,A,Sp=Sp,max.depth=max.depth,max.n=max.n,rho=rho,alpha=alpha,int.part=int.part,mixed.mode=mixed.mode,verb=verb);
		D = opt3d.density_grid(S,seq(A[1],A[2],length=n.grid),seq(A[3],A[4],length=n.grid),seq(A[5],A[6],length=n.grid));
		if( !is.null(Dp) ) { err = sqrt(mean((D-Dp)^2)); }
		Sp = S; Dp = D;
		S.iter[[n.iter]] = S;
		cat('iter',n.iter,sprintf('%.5f',err),'\n');
		n.iter = n.iter+1;
	}

	return( list(S=S,S.iter=S.iter) );	
}

# opt3d once
opt3d.one = function(X,A,Sp=NULL,max.depth=5,max.n=1,rho=0.5,alpha=0.5,int.part=FALSE,mixed.mode=FALSE,verb=FALSE) {

	fn.seed = sprintf('%d',round(runif(1)*10000));

	# options
	options = "-O";
	options = sprintf("%s -s %f-%f,%f-%f,%f-%f",options,A[1],A[2],A[3],A[4],A[5],A[6]);
	options = sprintf("%s -d %d -n %f -r %f -a %f",options,max.depth,max.n,rho,alpha);
	if( !is.null(Sp) ) {
		fn = paste("prior",fn.seed,sep='.');
		write.table(Sp[,1:8],file=fn,quote=FALSE,row.names=FALSE,col.names=FALSE,sep='\t');
		options = sprintf("%s -p %s",options,fn);
	}
	if( int.part ) { options = sprintf("%s -I",options); }
	if( mixed.mode ) { options = sprintf("%s -M",options); }
	out.fn = paste("out",fn.seed,sep='.');
	options = sprintf("%s -o %s",options,out.fn);

	# data file
	d.out = cbind(X[,1][,1],X[,1][,2],X[,2][,1],X[,2][,2],X[,3][,1],X[,3][,2]);
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
opt3d.p.area.sub = function(S,A) {
        P = 0; n = nrow(S);
        a1 = S[,1]; a1[a1<A[1]] = A[1];
        b1 = S[,2]; b1[b1>A[2]] = A[2];
        a2 = S[,3]; a2[a2<A[3]] = A[3];
        b2 = S[,4]; b2[b2>A[4]] = A[4];
        a3 = S[,5]; a3[a3<A[5]] = A[5];
        b3 = S[,6]; b3[b3>A[6]] = A[6];
        xx = b1-a1; xx[xx<0] = 0;
        yy = b2-a2; yy[yy<0] = 0;
        zz = b3-a3; zz[zz<0] = 0;
        P = sum(xx*yy*zz*S[,7]);

        return( P );
}

# probability in A, P(X in A)
opt3d.p.area = function(S,A) {
        if( all(is.finite(A)) ) { return( opt3d.p.area.sub(S,A) ); }
        if( is.infinite(A[1]) | is.infinite(A[3]) | is.infinite(A[5]) ) { return( 0 ); }

        # if A is an open space, we recalculate the density by replacing Inf to a large value
        t = S[,1:2]; r1 = range(t[is.finite(t)]);
        t = S[,3:4]; r2 = range(t[is.finite(t)]);
        t = S[,5:6]; r3 = range(t[is.finite(t)]);
        L = 10*max(c(r1,r2,r3));
        if( is.infinite(A[2]) ) { S[is.infinite(S[,2]),2] = L; A[2] = L; }
        if( is.infinite(A[4]) ) { S[is.infinite(S[,4]),4] = L; A[4] = L; }
        if( is.infinite(A[6]) ) { S[is.infinite(S[,6]),6] = L; A[6] = L; }
        S[,7] = S[,8]/(S[,2]-S[,1])/(S[,4]-S[,3])/(S[,6]-S[,5]);

        return( opt3d.p.area.sub(S,A) );
}

# survival probablity P(T1>=t1, T2>=t2)
opt3d.p = function(T,t1,t2,t3) {
        return( opt3d.p.area(T,c(t1,Inf,t2,Inf,t3,Inf)) );
}


# obtain a marginal distribution
opt3d.marginal = function(S,IND) {

        SM = NULL;
        if( IND == 1 ) {
                t = unique(sort(as.vector(S[,1:2])));
                for( i in 2:length(t) ) {
                        p = opt3d.p.area(S,c(t[i-1],t[i],0,Inf,0,Inf));
                        SM = rbind(SM,c(t[i-1],t[i],p/(t[i]-t[i-1]),p));
                }
        } else if( IND == 2 ) {
                t = unique(sort(as.vector(S[,3:4])));
                for( i in 2:length(t) ) {
                        p = opt3d.p.area(S,c(0,Inf,t[i-1],t[i],0,Inf));
                        SM = rbind(SM,c(t[i-1],t[i],p/(t[i]-t[i-1]),p));
                }
        } else if( IND == 3 ) {
                t = unique(sort(as.vector(S[,5:6])));
                for( i in 2:length(t) ) {
                        p = opt3d.p.area(S,c(0,Inf,0,Inf,t[i-1],t[i]));
                        SM = rbind(SM,c(t[i-1],t[i],p/(t[i]-t[i-1]),p));
                }
        } else if( IND == 12 ) {
                t1 = unique(sort(as.vector(S[,1:2])));
                t2 = unique(sort(as.vector(S[,3:4])));
                for( i in 2:length(t1) ) {
                        for( j in 2:length(t2) ) {
                                p = opt3d.p.area(S,c(t1[i-1],t1[i],t2[j-1],t2[j],0,Inf));
                                SM = rbind(SM,c(t1[i-1],t1[i],t2[j-1],t2[j],p/(t1[i]-t1[i-1])/(t2[j]-t2[j-1]),p));
                        }
                }
       } else if( IND == 23 ) {
                t1 = unique(sort(as.vector(S[,3:4])));
                t2 = unique(sort(as.vector(S[,5:6])));
                for( i in 2:length(t1) ) {
                        for( j in 2:length(t2) ) {
                                p = opt3d.p.area(S,c(0,Inf,t1[i-1],t1[i],t2[j-1],t2[j]));
                                SM = rbind(SM,c(t1[i-1],t1[i],t2[j-1],t2[j],p/(t1[i]-t1[i-1])/(t2[j]-t2[j-1]),p));
                        }
                }
        } else if( IND == 13 ) {
                t1 = unique(sort(as.vector(S[,1:2])));
                t2 = unique(sort(as.vector(S[,5:6])));
                for( i in 2:length(t1) ) {
                        for( j in 2:length(t2) ) {
                                p = opt3d.p.area(S,c(t1[i-1],t1[i],0,Inf,t2[j-1],t2[j]));
                                SM = rbind(SM,c(t1[i-1],t1[i],t2[j-1],t2[j],p/(t1[i]-t1[i-1])/(t2[j]-t2[j-1]),p));
                        }
                }
        }

        return( SM );
}


##############################################################################
# plotting functions
##############################################################################

##############################################################################
# Misc.
##############################################################################

opt3d.density_grid = function(S,xgrid,ygrid,zgrid) {

	xgrid = c(xgrid,Inf);
	ygrid = c(ygrid,Inf);
	zgrid = c(zgrid,Inf);

	D = rep(0,(length(xgrid)-1)*(length(ygrid)-1)*(length(zgrid)-1));
	DN = rep("",(length(xgrid)-1)*(length(ygrid)-1)*(length(zgrid)-1));
	ii = 1;
	for( i in 2:length(xgrid) ) {
		for( j in 2:length(ygrid) ) {
			for( k in 2:length(zgrid) ) {
				D[ii] = opt3d.p.area(S,c(xgrid[i-1],xgrid[i],ygrid[j-1],ygrid[j],zgrid[k-1],zgrid[k]));
				DN[ii] = sprintf("[%.2f,%.2f)x[%.2f,%.2f)x[%.2f,%.2f)", xgrid[i-1],xgrid[i], ygrid[j-1],ygrid[j], zgrid[k-1],zgrid[k]); 
				ii = ii+1;
			}
		}
	}
	names(D) = DN;

	return( D );
}

