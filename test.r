N=100000; # taille de l'échantillon

X = rexp(N); # echantillon

# Bootstrap Efron
B=999; maxB=c(); Y=seq(-5,6,0.01);
for (i in 1:B) maxB[i] = max(sample(X,size=N,replace=T));
plot(Y,ecdf(maxB-log(m))(Y),type='l',ylim=c(0,1)); lines(Y,exp(-exp(-Y)),type='l');

# Moon
B=999; m=1000; maxB=c(); Y=seq(-5,6,0.01);
for (i in 1:B) maxB[i] = max(sample(X,size=m,replace=T));
plot(Y,ecdf(maxB-log(N))(Y),type='l',ylim=c(0,1)); lines(Y,exp(-exp(-Y)),type='l');

# subSampling
B=999; m=1000; maxB=c(); Y=seq(-5,6,0.01);
for (i in 1:B) maxB[i] = max(sample(X,size=m,replace=F));
plot(Y,ecdf(maxB-log(m))(Y),type='l',ylim=c(0,1)); lines(Y,exp(-exp(-Y)),type='l');

# estimation des coefficients de normalisation
m=1000; bm = quantile(X,1-1/m,type=1); am = quantile(X,1-1/(exp(1)*m),type=1) - bm; 

B=999; maxB=c(); Y=seq(-5,6,0.01);
for (i in 1:B) maxB[i] = max(sample(X,size=m,replace=T));
plot(Y,ecdf((maxB-bm)/am)(Y),type='l',ylim=c(0,1)); lines(Y,exp(-exp(-Y)),type='l');

B=999; maxB=c(); Y=seq(-5,6,0.01);
for (i in 1:B) maxB[i] = max(sample(X,size=m,replace=F));
plot(Y,ecdf((maxB-bm)/am)(Y),type='l',ylim=c(0,1)); lines(Y,exp(-exp(-Y)),type='l');



N=100000; # taille de l'échantillon
m=1000;

X = rlogis(N); # echantillon

# Bootstrap Efron
B=999; maxB=c(); Y=seq(-5,6,0.01);
for (i in 1:B) maxB[i] = max(sample(X,size=N,replace=T));
bn = quantile(X,1-1/N,type=1); an = quantile(X,1-1/(exp(1)*N),type=1) - bn; 
plot(Y,ecdf((maxB-bn)/an)(Y),type='l',ylim=c(0,1)); lines(Y,exp(-exp(-Y)),type='l');

# estimation des coefficients de normalisation
m=1000; bm = quantile(X,1-1/m,type=1); am = quantile(X,1-1/(exp(1)*m),type=1) - bm; 

B=999; m=1000; maxB=c(); Y=seq(-5,6,0.01);
for (i in 1:B) maxB[i] = max(sample(X,size=m,replace=T));
plot(Y,ecdf((maxB-bm)/am)(Y),type='l',ylim=c(0,1)); lines(Y,exp(-exp(-Y)),type='l');

B=999; m=1000; maxB=c(); Y=seq(-5,6,0.01);
for (i in 1:B) maxB[i] = max(sample(X,size=m,replace=F));
plot(Y,ecdf((maxB-bm)/am)(Y),type='l',ylim=c(0,1)); lines(Y,exp(-exp(-Y)),type='l');





N=100000; # taille de l'échantillon
m=1000;

X = rnorm(N); # echantillon

# Bootstrap Efron
B=999; maxB=c(); Y=seq(-5,6,0.01);
for (i in 1:B) maxB[i] = max(sample(X,size=N,replace=T));
bn = quantile(X,1-1/N,type=1); an = quantile(X,1-1/(exp(1)*N),type=1) - bn; 
plot(Y,ecdf((maxB-bn)/an)(Y),type='l',ylim=c(0,1)); lines(Y,exp(-exp(-Y)),type='l');

# estimation des coefficients de normalisation
m=1000; bm = quantile(X,1-1/m,type=1); am = quantile(X,1-1/(exp(1)*m),type=1) - bm; 

B=999; m=1000; maxB=c(); Y=seq(-5,6,0.01);
for (i in 1:B) maxB[i] = max(sample(X,size=m,replace=T));
plot(Y,ecdf((maxB-bm)/am)(Y),type='l',ylim=c(0,1)); lines(Y,exp(-exp(-Y)),type='l');

B=999; m=1000; maxB=c(); Y=seq(-5,6,0.01);
for (i in 1:B) maxB[i] = max(sample(X,size=m,replace=F));
plot(Y,ecdf((maxB-bm)/am)(Y),type='l',ylim=c(0,1)); lines(Y,exp(-exp(-Y)),type='l');





# Calcul distance de Kolmogorov

X = rexp(N); # echantillon

B=9999; m=20; maxB=c(); Y=seq(-5,6,0.01);
for (i in 1:B) maxB[i] = max(sample(X,size=m,replace=F));
plot(Y,ecdf(maxB-log(m))(Y),type='l',ylim=c(0,1)); lines(Y,exp(-exp(-Y)),type='l');

Y=sort(maxB-log(m)); max(abs((1:length(Y))/length(Y)-exp(-exp(-Y))));
