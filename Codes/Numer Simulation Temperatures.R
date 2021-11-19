setwd("C:/Users/Sidney Williams/Documents/New folder/Spherical Stefan/")
library(latex2exp)
library(pracma)

LWD=3
case=1

temp=read.csv("tempm1.txt")
temp2=read.csv("tempm2.txt")
temp3=read.csv("tempm3.txt")
r=seq(1/1400,1,by=1/1400)
temps=c()
temps2=c()
temps3=c()

L=length(temp[,1])

for(i in 1:L){
	temps=append(temps,temp[i,1])
	temps2=append(temps2,temp2[i,1])
	temps3=append(temps3,temp3[i,1])
}

par(mar=c(5,5,4,2)+.1)
plot(r,temps,type='l',,xlab=TeX('$\\eta$-nondimensional radius'),ylab=TeX('$\\theta_{liq},\\theta_{sol}$-nondimensional temperature'),lwd=LWD,cex.lab=1.6,cex.axis=1.3,ylim=c(0,1.3))
points(r,temps2,type='l',lwd=LWD)
points(r,temps3,type='l',lwd=LWD)

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

m=0
L=length(t3)
for(i in 1:L){
	if(t3[i]==0.26233994E+00-tstar){
		m=i
	}

}



Q=7
St=1
nmax=10
xi0=0.01

j0=function(z,y){
	if (y==0){
		a=1
	}
	else{
		a=sin(z*y)/(z*y)
	}
return(a)
}

j1=function(z,y){
	a=sin(z*y)/(z*y)^2-cos(z*y)/(z*y)
return(a)
}

y0=function(z,y){
	a=-cos(z*y)/(z*y)
}

y1=function(z,y){
	a=-cos(z*y)/(z*y)^2-sin(z*y)/(z*y)
}

fn=function(y,z){
	a=j0(z,y)-j0(y,1)/y0(y,1)*y0(z,y)
return(a)
}

dfn=function(y,z){
	a=y*(j0(y,1)/y0(y,1)*y1(z,y)-j1(z,y))
return(a)
}

LIQss=function(x,y){
	a=-1/6*Q*x^2+1/6*Q*y^2+1
return(a)
}

C1=function(y){
	a=(1/6*Q*y^2-1/6*Q+1)/(1-1/y)
return(a)
}

C2=function(y){
	a=C1(y)+1/6*Q
return(a)
}

SOLss=function(x,y){
	a=-1/6*Q*x^2-C1(y)/x+C2(y)
return(a)
}

if(case==1){
solid_profile=function(zeta,delta_tau){
L=zeta
	eta=linspace(zeta,1,50)
		
	theta1=zeros(1,50)
    	theta_ss=zeros(1,50)
    	for(i in 1:50){
        #if (eta[i]<zeta_old){
            theta1[i]=1
        #} else{
            #theta1[i]=interp1(eta_old,theta,eta[i])
        #}
	theta_ss[i]=SOLss(eta[i],zeta)
    }
    
    URTS=zeros(1,nmax)
    B=zeros(1,nmax)
    intagand=zeros(1,50)
    intagand2=zeros(1,50)
    for(n in 1:nmax){
	  URTS[n]=n*pi/(zeta-1)
        for(i in 1:50){
            intagand[i]=((xi0^2)-((eta[i]^2))+1-theta_ss[i])*fn(URTS[n],eta[i])*(eta[i])^2
        }
	  for(r in 1:50){
		intagand2[r]=(fn(URTS[n],eta[r]))^2*(eta[r])^2
	  }
        #inttrap is a numerical integration, B_n (equ. 29)
        B[n]=trapz(eta,intagand)/trapz(eta,intagand2)
    }
dtdn_part=zeros(1,nmax)
dtdn_tr=zeros(1,50)
dtdn=zeros(1,50)
dtdn_ss=zeros(1,50)
for(i in 1:50){
    for(n in 1:nmax){
        dtdn_part[n]=B[n]*fn(URTS[n],eta[i])*(exp(-(URTS[n])^2*delta_tau))
    }
dtdn_tr[i]=sum(dtdn_part)
dtdn_ss[i]=-1/6*Q*(eta[i])^2-(C1(zeta))*1/(eta[i])+C2(zeta)
dtdn[i]=dtdn_tr[i]+dtdn_ss[i]
}
return(dtdn)
}

liquid_profile=function(zeta,delta_tau){
L=zeta
eta=linspace(0,zeta, 50)

theta_ss=zeros(1,50)
for (i in 1:50){
	theta_ss[i]=LIQss(eta[i],zeta)
}

lambda=zeros(1,nmax)
A=zeros(1,nmax)
intagrand=zeros(1,50)
intagrand2=zeros(1,50)
for(n in 1:nmax){
	lambda[n]=pi*n/L
	for(i in 1:50){
		intagrand[i]=(theta_liquid[i]-theta_ss[i])*j0(lambda[n],eta[i])*(eta[i])^2
	}
	for(r in 1:50){
		intagrand2[r]=(j0(lambda[n],eta[r]))^2*(eta[r])^2

	}
	A[n]=trapz(eta,intagrand)/trapz(eta,intagrand2)
}
dtdn_part=zeros(1,nmax)
dtdn_tr=zeros(1,50)
dtdn=zeros(1,50)
dtdn_ss=zeros(1,50)
for(i in 1:50){
    for(n in 1:nmax){
        dtdn_part[n]=A[n]*j0(lambda[n],eta[i])*exp(-(lambda[n])^2*delta_tau)
}
dtdn_tr[i]=sum(dtdn_part)
dtdn_ss[i]=-1/6*Q*(eta[i])^2+1/6*Q*zeta^2+1
dtdn[i]=dtdn_tr[i]+dtdn_ss[i]
}
return(dtdn)
}

eta_solid=linspace(xi0,1, 50)
eta_liquid=linspace(0,xi0, 50)

#preallocates space to save time
theta_solid=zeros(1,50)
theta_liquid=zeros(1,50)

for(i in 1:50){
	theta_solid[i]=SOLss(eta_solid[i],xi0)

	theta_liquid[i]=LIQss(eta_liquid[i],xi0)
}

liquid_profile2=liquid_profile(z3[149],t3[149])
solid_profile2=solid_profile(z3[149],t3[149])
liquid_profile3=liquid_profile(z3[260],t3[260])
solid_profile3=solid_profile(z3[260],t3[260])
liquid_profile4=liquid_profile(z3[511],t3[511])
solid_profile4=solid_profile(z3[511],t3[511])

eta_temp_liq=linspace(0,z3[149],50)
eta_temp_sol=linspace(z3[149],1,50)
eta_temp_liq1=linspace(0,z3[260],50)
eta_temp_sol1=linspace(z3[260],1,50)
eta_temp_liq2=linspace(0,z3[511],50)
eta_temp_sol2=linspace(z3[511],1,50)

points(eta_temp_liq,liquid_profile2,type="l",lty=2,col="green",lwd=LWD)
points(eta_temp_sol,solid_profile2,type="l", lty=2,col="red",lwd=LWD)
points(eta_temp_liq1,liquid_profile3,type="l",lty=2,col="green",lwd=LWD)
points(eta_temp_sol1,solid_profile3,type="l", lty=2,col="red",lwd=LWD)
points(eta_temp_liq2,liquid_profile4,type="l",lty=2,col="green",lwd=LWD)
points(eta_temp_sol2,solid_profile4,type="l", lty=2,col="red",lwd=LWD)

#arrows(x0=0.7,y0=0.8,x1=0.9,y1=1.1)
text(x=0.01,y=0.94,labels=TeX('$\\tau$'),cex=1.3)
text(x=0.08,y=0.94,labels="=0.1021",cex=1)
text(x=0.01,y=1.08,labels=TeX('$\\tau$'),cex=1.3)
text(x=0.08,y=1.08,labels="=0.1905",cex=1)
text(x=0.01,y=1.2,labels=TeX('$\\tau$'),cex=1.3)
text(x=0.08,y=1.2,labels="=2.719",cex=1)

legend(0.1, 0.3, legend=c("Liquid Region", "Solid Region","Numerical Solution"),col=c("green", "red","black"), lty=c(2,2,1), cex=1,lwd=LWD,pt.cex=1)
}

if(case==2){
solid_profile=function(zeta,delta_tau){
L=zeta
	tg=50
	eta=linspace(zeta,1,tg)
		
	theta1=zeros(1,tg)
    	theta_ss=zeros(1,tg)
    	for(i in 1:tg){
        #if (eta[i]<zeta_old){
            theta1[i]=1
        #} else{
            #theta1[i]=interp1(eta_old,theta,eta[i])
        #}
	theta_ss[i]=SOLss(eta[i],zeta)
    }
    
    URTS=zeros(1,nmax)
    B=zeros(1,nmax)
    intagand=zeros(1,tg)
    intagand2=zeros(1,tg)
    for(n in 1:nmax){
	  URTS[n]=n*pi/(zeta-1)
        for(i in 1:tg){
            intagand[i]=(theta1[i]-theta_ss[i])*fn(URTS[n],eta[i])*(eta[i])^2
        }
	  for(r in 1:tg){
		intagand2[r]=(fn(URTS[n],eta[r]))^2*(eta[r])^2
	  }
        #inttrap is a numerical integration, B_n (equ. 29)
        B[n]=trapz(eta,intagand)/trapz(eta,intagand2)
    }
dtdn_part=zeros(1,nmax)
dtdn_tr=zeros(1,tg)
dtdn=zeros(1,tg)
dtdn_ss=zeros(1,tg)
for(i in 1:tg){
    for(n in 1:nmax){
        dtdn_part[n]=B[n]*fn(URTS[n],eta[i])*(exp(-(URTS[n])^2*delta_tau))
    }
dtdn_tr[i]=sum(dtdn_part)
dtdn_ss[i]=-1/6*Q*(eta[i])^2-(C1(zeta))*1/(eta[i])+C2(zeta)
dtdn[i]=dtdn_tr[i]+dtdn_ss[i]
}
return(dtdn)
}

liquid_profile=function(zeta,delta_tau){
L=zeta
tg=50

eta=linspace(0,zeta, tg)
theta_ss=zeros(1,tg)
for (i in 1:tg){
	theta_ss[i]=LIQss(eta[i],zeta)
}

lambda=zeros(1,nmax)
A=zeros(1,nmax)
intagrand=zeros(1,tg)
intagrand2=zeros(1,tg)
for(n in 1:nmax){
	lambda[n]=pi*n/L
	for(i in 1:tg){
		intagrand[i]=(theta_liquid[i]-theta_ss[i])*j0(lambda[n],eta[i])*(eta[i])^2
	}
	for(r in 1:tg){
		intagrand2[r]=(j0(lambda[n],eta[r]))^2*(eta[r])^2

	}
	A[n]=trapz(eta,intagrand)/trapz(eta,intagrand2)
}
dtdn_part=zeros(1,nmax)
dtdn_tr=zeros(1,tg)
dtdn=zeros(1,tg)
dtdn_ss=zeros(1,tg)
for(i in 1:tg){
    for(n in 1:nmax){
        dtdn_part[n]=A[n]*j0(lambda[n],eta[i])*exp(-(lambda[n])^2*delta_tau)
}
dtdn_tr[i]=sum(dtdn_part)
dtdn_ss[i]=-1/6*Q*(eta[i])^2+1/6*Q*zeta^2+1
dtdn[i]=dtdn_tr[i]+dtdn_ss[i]
}
return(dtdn)
}

eta_solid=linspace(xi0,1, 50)
eta_liquid=linspace(0,xi0, 50)

#preallocates space to save time
theta_solid=zeros(1,50)
theta_liquid=zeros(1,50)

for(i in 1:50){
	theta_solid[i]=SOLss(eta_solid[i],xi0)

	theta_liquid[i]=LIQss(eta_liquid[i],xi0)
}

liquid_profile2=liquid_profile(z3[547],t3[547])
solid_profile2=solid_profile(z3[547],t3[547])
liquid_profile3=liquid_profile(z3[694],t3[694])
solid_profile3=solid_profile(z3[694],t3[694])
liquid_profile4=liquid_profile(z3[860],t3[860])
solid_profile4=solid_profile(z3[860],t3[860])

eta_temp_liq=linspace(0,z3[547],50)
eta_temp_sol=linspace(z3[547],1,50)
eta_temp_liq1=linspace(0,z3[694],50)
eta_temp_sol1=linspace(z3[694],1,50)
eta_temp_liq2=linspace(0,z3[860],50)
eta_temp_sol2=linspace(z3[860],1,50)

points(eta_temp_liq,liquid_profile2,type="l",lty=2,col="green",lwd=LWD)
points(eta_temp_sol,solid_profile2,type="l", lty=2,col="red",lwd=LWD)
points(eta_temp_liq1,liquid_profile3,type="l",lty=2,col="green",lwd=LWD)
points(eta_temp_sol1,solid_profile3,type="l", lty=2,col="red",lwd=LWD)
points(eta_temp_liq2,liquid_profile4,type="l",lty=2,col="green",lwd=LWD)
points(eta_temp_sol2,solid_profile4,type="l", lty=2,col="red",lwd=LWD)

#arrows(x0=1,y0=1.3,x1=0.75,y1=1)
text(x=0.5,y=1.2,labels=TeX('$\\tau$'),cex=1.3)
text(x=0.57,y=1.2,labels="=0.2623",cex=1)
text(x=0.11,y=1.33,labels=TeX('$\\tau$'),cex=1.3)
text(x=0.18,y=1.33,labels="=0.4950",cex=1)
text(x=0.3,y=0.9,labels=TeX('$\\tau$'),cex=1.3)
text(x=0.37,y=0.9,labels="=3.0002",cex=1)
legend(0.1, 0.7, legend=c("Liquid Region", "Solid Region","Numerical Solution"),col=c("green", "red","black"), lty=c(2,2,1), cex=1,lwd=LWD,pt.cex=1)
}
######################
#Recorded Times
#m1=0.19054427E+00, 260
#m2=0.10205394E+00, 149
#m3=0.27190586E+01, 511
#s1=0.26233994E+00, 547
#s2=0.49503619E+00, 694
#s3=0.30002895E+01, 860