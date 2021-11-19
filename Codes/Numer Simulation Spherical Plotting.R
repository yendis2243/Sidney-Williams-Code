setwd("C:/Users/Sidney Williams/Documents/New folder/Spherical Stefan/")
library(latex2exp)

LWD=3

zeta=read.csv("s(time)m0.01.txt")
tau=read.csv("timem0.01.txt")
L=length(zeta[,1])
z=c()
t=c()
tstar=0
for(i in 1:L){
	if((zeta[i,1] != 0)&(zeta[i,1]>0)){
		if(length(z)==0){
			tstar=tau[i,1]
		}
		t=append(t,tau[i,1]-tstar)
		z=append(z,zeta[i,1])	
	}
	if((tau[i,1]-tstar)>=10){
		break
	}
}

zeta=read.csv("s(time)m0.1.txt")
tau=read.csv("timem0.1.txt")
L=length(zeta[,1])
z2=c()
t2=c()
tstar=0
for(i in 1:L){
	if((zeta[i,1] != 0)&(zeta[i,1]>0)){
		if(length(z2)==0){
			tstar=tau[i,1]
		}
		t2=append(t2,tau[i,1]-tstar)
		z2=append(z2,zeta[i,1])	
	}
	if((tau[i,1]-tstar)>=10){
		break
	}
}

zeta=read.csv("s(time)m1.txt")
tau=read.csv("timem1.txt")
L=length(zeta[,1])
z3=c()
t3=c()
tstar=0
for(i in 1:L){
	if((zeta[i,1] != 0)&(zeta[i,1]>0)){
		if(length(z3)==0){
			tstar=tau[i,1]
		}
		t3=append(t3,tau[i,1]-tstar)
		z3=append(z3,zeta[i,1])	
	}
	if((tau[i,1]-tstar)>=10){
		break
	}
}

zeta=read.csv("s(time)m10.txt")
tau=read.csv("timem10.txt")
L=length(zeta[,1])
z4=c()
t4=c()
tstar=0
for(i in 1:L){
	if((zeta[i,1] != 0)&(zeta[i,1]>0)){
		if(length(z4)==0){
			tstar=tau[i,1]
		}
		t4=append(t4,tau[i,1]-tstar)
		z4=append(z4,zeta[i,1])	
	}
	if((tau[i,1]-tstar)>=10){
		break
	}
}

t3=append(t3,10)
z3=append(z3,sqrt(1-6/7))
t4=append(t4,10)
z4=append(z4,sqrt(1-6/7))

par(mar=c(5,5,4,2)+.1)
plot(t,z,col="red",type='l',ylim=c(0,0.4),xlim=c(0,0.6),lwd=LWD,xlab=TeX('$\\tau$-nondimensional time'),ylab=TeX('$\\zeta$-nondimensional interface'),cex.lab=1.6,cex.axis=1.3)
points(t2,z2,col="red",type='l',lwd=LWD)
points(t3,z3,col="red",type='l',lwd=LWD)
points(t4,z4,col="red",type='l',lwd=LWD)


setwd("C:/Users/Sidney Williams/Pictures/Spherical/Spherical Data")

Sts=c(0.01,0.1,1,10)
titlem=c("Spherical_Data_Melting_Quasi_St=","Spherical_Data_Melting_St=")
titles=c("Spherical_Data_Solidification_Quasi_St=","Spherical_Data_Solidification_St=")
file_title=function(title,St){
	x=c(title,St,".txt")
	myfile=toString(x)
return(myfile)
}

dat1m=read.csv(file_title(titlem[1],Sts[1]))
dat2m=read.csv(file_title(titlem[1],Sts[2]))
dat3m=read.csv(file_title(titlem[1],Sts[3]))
dat4m=read.csv(file_title(titlem[1],Sts[4]))
dat5m=read.csv(file_title(titlem[2],Sts[1]))
dat6m=read.csv(file_title(titlem[2],Sts[2]))
dat7m=read.csv(file_title(titlem[2],Sts[3]))
dat8m=read.csv(file_title(titlem[2],Sts[4]))

dat1s=read.csv(file_title(titles[1],Sts[1]))
dat2s=read.csv(file_title(titles[1],Sts[2]))
dat3s=read.csv(file_title(titles[1],Sts[3]))
dat4s=read.csv(file_title(titles[1],Sts[4]))
dat5s=read.csv(file_title(titles[2],Sts[1]))
dat6s=read.csv(file_title(titles[2],Sts[2]))
dat7s=read.csv(file_title(titles[2],Sts[3]))
dat8s=read.csv(file_title(titles[2],Sts[4]))

melt_0.01_quasi_tau=dat1m[1:10001,2]
melt_0.01_quasi_eta=dat1m[1:10001,3]
melt_0.1_quasi_tau=dat2m[1:10001,2]
melt_0.1_quasi_eta=dat2m[1:10001,3]
melt_1_quasi_tau=dat3m[1:10001,2]
melt_1_quasi_eta=dat3m[1:10001,3]
melt_10_quasi_tau=dat4m[1:10001,2]
melt_10_quasi_eta=dat4m[1:10001,3]
melt_0.01_tau=dat5m[1:10001,2]
melt_0.01_eta=dat5m[1:10001,3]
melt_0.1_tau=dat6m[1:10001,2]
melt_0.1_eta=dat6m[1:10001,3]
melt_1_tau=dat7m[1:10001,2]
melt_1_eta=dat7m[1:10001,3]
melt_10_tau=dat8m[1:10001,2]
melt_10_eta=dat8m[1:10001,3]

solid_0.01_quasi_tau=dat1s[1:10001,2]
solid_0.01_quasi_eta=dat1s[1:10001,3]
solid_0.1_quasi_tau=dat2s[1:10001,2]
solid_0.1_quasi_eta=dat2s[1:10001,3]
solid_1_quasi_tau=dat3s[1:10001,2]
solid_1_quasi_eta=dat3s[1:10001,3]
solid_10_quasi_tau=dat4s[1:10001,2]
solid_10_quasi_eta=dat4s[1:10001,3]
solid_0.01_tau=dat5s[1:10001,2]
solid_0.01_eta=dat5s[1:10001,3]
solid_0.1_tau=dat6s[1:10001,2]
solid_0.1_eta=dat6s[1:10001,3]
solid_1_tau=dat7s[1:10001,2]
solid_1_eta=dat7s[1:10001,3]
solid_10_tau=dat8s[1:10001,2]
solid_10_eta=dat8s[1:10001,3]

points(melt_0.01_quasi_tau,melt_0.01_quasi_eta,type="l",lwd=LWD,lty=3,col="purple")
points(melt_0.1_quasi_tau,melt_0.1_quasi_eta,type="l",lwd=LWD,lty=3,col="purple")
points(melt_1_quasi_tau,melt_1_quasi_eta,type="l",lwd=LWD,lty=3,col="purple")
points(melt_10_quasi_tau,melt_10_quasi_eta,type="l",lwd=LWD,lty=3,col="purple")
points(melt_0.01_tau,melt_0.01_eta,type="l",lwd=LWD,lty=2)
points(melt_0.1_tau,melt_0.1_eta,type="l",lwd=LWD,lty=2)
points(melt_1_tau,melt_1_eta,type="l",lwd=LWD,lty=2)
points(melt_10_tau,melt_10_eta,type="l",lwd=LWD,lty=2)


text(x=0.5,y=0.06,labels=TeX('$St=0.01$'),cex=1.2)
text(x=0.4,y=0.14,labels=TeX('$St=0.1$'),cex=1.2)
text(x=0.3,y=0.31,labels=TeX('$St=1$'),cex=1.2)
text(x=0.12,y=0.4,labels=TeX('$St=10$'),cex=1.2)

legend(0.35, 0.25, legend=c("Numerical Simulation", "Series Solution", "Quasi-Static Solution"), col=c("red", "black", "purple"), lty=1:3, cex=1,lwd=LWD)