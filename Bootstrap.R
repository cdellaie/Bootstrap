#Valeurs extremes

#package qui propose une fonction simulant des pareto
#rpareto(N,location,shape)
library(VGAM)

# Methode sous echantillonage par MC pour obtenir une estimation de la loi du maximum
# Tirage sans remise de B echantillons de x=data, bn est la taille des sous-echantillons

BootMC=function(x,B,bn){
N=length(x);
Tstar=rep(0,B);
for (k in 1:B ){
	i=sample(1:N,bn,replace=FALSE);
	Tstar[k]=max(x[i]);
	}
return(Tstar);
}

#Premier exemple : loi uniforme

B=100;
N=100000;
theta=8;
X=runif(N,0,theta);

T=cummax(X);
S=cummin(X);

Tstar=BootMC(X,1000,100);
plot(Tstar[order(Tstar)],(1/B)*(1:B),col='blue',type='s',xlim=c(0,10))
title(main='FDR empirique du max d'une uniforme obtenue par Bootstrap')


#loi normale
B=1000;
N=100000;
m=0;
sigma=1;
Z=rnorm(N,m,sigma);

Zmax=BootMC(Z,B,100);
plot(Zmax[order(Zmax)],(1/B)*(1:B),col='blue',type='s')
title(main='FDR empirique du max d'une normale obtenue par Bootstrap')

#loi de Pareto : queue epaisse
B=1000;
N=100000;
gamma=1;
alpha=1;
P=rpareto(N,gamma,alpha);

Pmax=BootMC(P,B,100);
plot(Pmax[order(Pmax)],(1/B)*(1:B),col='blue',type='s')
title(main='FDR empirique du max d une pareto obtenue par Bootstrap')





