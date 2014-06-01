		####################
		# Valeurs extremes #
		####################

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

##################################
# Premier exemple : loi uniforme #
##################################

B=100;
N=100000;
theta=8;
X=runif(N,0,theta);

T=cummax(X);
S=cummin(X);

Tstar=BootMC(X,1000,100);
plot(Tstar[order(Tstar)],(1/B)*(1:B),col='blue',type='s',xlim=c(0,10))
title(main='FDR empirique du max d une uniforme obtenue par Bootstrap')

#################
#  loi normale ##
#################

B=1000;
N=100000;
m=0;
sigma=1;
Z=rnorm(N,m,sigma);

Zmax=BootMC(Z,B,100);
plot(Zmax[order(Zmax)],(1/B)*(1:B),col='blue',type='s')
title(main='FDR empirique du max d une normale obtenue par Bootstrap')

#################################
# loi de Pareto : queue epaisse #
#################################

B=1000;
N=100000;
gamma=1;
alpha=1;
P=rpareto(N,gamma,alpha);

Pmax=BootMC(P,B,100);
plot(Pmax[order(Pmax)],(1/B)*(1:B),col='blue',type='s')
title(main='FDR empirique du max d une pareto obtenue par Bootstrap')

#########################
# 16/05/2014            #
# Pdf de Guillaume      #
# Test sur une normale  #
#########################

e=exp(-1);
B=1000;
N=100000;
bn=100;
m=0;
sigma=1;
t=seq(0,10,by= 0.1);

Z=rnorm(N,m,sigma);
a=quantile(Z,1-(e*N)^(-1))-quantile(Z,1-(1/N));
b=quantile(Z,1-(1/N));

Zmax=(BootMC(Z,B,bn)-b)/a;
Ztild=Zmax[order(Zmax)];
plot(Ztild,(1/B)*(1:B),col='blue',type='s')
lines(t,exp(-exp(-t)),col=rgb(0.6,0.3,0.1),type='l');
title(main='FDR empirique du max d une normale obtenue par Bootstrap')

###########################################################
# 01/06/2014                                              #
# estimation du gamma par regression, articles de Bertail #
###########################################################

Subsample=function(x,bn){
N=length(x);
q=N-bn+1;
Tstar=rep(0,q);
for (k in 1:q ){
	i=k:(k+bn-1);
	Tstar[k]=max(x[i]);
	}
return(Tstar);
}

N=10;
I=16;
bn= N**seq(1,4,by=0.2);
y=rep(0,I);
t=0.5;

############
# Uniforme #
############

theta=10;
U=runif(N**5,0,theta);

for (k in 1:I ){
	Z=Subsample(U,bn[k]);
	y[k]=log(abs(quantile(Z,t)));
	}

#regression
plot(y ~ log(bn)+0,pch=4,col="blue"); 
lin<-lm(y ~ log(bn)+0);
abline(lin, col="red");
title("Estimation par régression");

C=sum( log(bn)- mean(log(bn)) ) ;
gamma=((sum(y-mean(y)))*C)/C**2;    #0.4

##############
# Gaussienne #
##############

X=rnorm(N**5,0,1);
for (k in 1:I ){
	Z=Subsample(X,bn[k]);
	y[k]=log(abs(quantile(Z,t)));
	}

#regression
plot(y ~ log(bn),pch=4,col="blue"); 
lin<-lm(y ~ log(bn));
abline(lin, col="red");
title("Estimation par régression");

# à la main (on trouve un coefficient différent !)
C=sum( log(bn)- mean(log(bn)) ) ;
gamma=((sum(y-mean(y)))*C)/C**2;    #0.3

###
## Section : The extreme value statistics revisited
## Estimation par difference de quantiles
###
t1=0.5;
t2=0.25;
for (k in 1:I ){
	Z=Subsample(X,bn[k]);
	y[k]=log(abs(quantile(Z,t1)-quantile(Z,t2)));
	}

plot(y ~ log(bn),pch=4,col="blue"); 
lin<-lm(y ~ log(bn));
abline(lin, col="red");
title("Estimation par régression");
# intercept : -0.8158 gamma: -0.1108 gamma a la main : -0.3




