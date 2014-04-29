#valeurs extremes
#premier exemple : loi uniforme

N=100000;
bn=100;
theta=8;
X=runif(N,0,theta);

T=cummax(X);
S=cummin(X);

#sous-echantillonage : methode de Monte Carlo
#tirage sans remise de B echantillons
B=1000;
Tstar=rep(0,B);
for (k in 1:B ){
	i=sample(1:N,bn,replace=FALSE);
	Tstar[k]=max(X[i]);
	}

plot(Tstar[order(Tstar)],(1/B)*(1:B),col='blue',type='s',xlim=c(0,10))
title(main='FDR empirique du max obtenue par Bootstrap')