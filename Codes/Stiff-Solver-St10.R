library(pracma)
library(latex2exp)
library(deSolve)
Q=7
St=1
nmax=10
xi0=0.01
case=2

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

solid_temp=function(delta_tau,zeta,zeta_old,theta){

	L=zeta
	eta=linspace(zeta,1,50)
	eta_old=linspace(zeta_old,1,50)
	
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
    for(n in 1:nmax){
        dtdn_part[n]=B[n]*dfn(URTS[n],zeta)*(exp(-(URTS[n])^2*delta_tau))
    }

    dtdn_tr=sum(dtdn_part)
    dtdn_ss=-1/3*Q*zeta+(C1(zeta))*1/zeta^2
    dtdn=dtdn_tr+dtdn_ss
return(dtdn)
}

liquid_temp=function(delta_tau,zeta,theta){

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
		intagrand[i]=(1-theta_ss[i])*j0(lambda[n],eta[i])*(eta[i])^2
	}
	for(r in 1:50){
		intagrand2[r]=(j0(lambda[n],eta[r]))^2*(eta[r])^2

	}
	A[n]=trapz(eta,intagrand)/trapz(eta,intagrand2)
}
dtdn_part=zeros(1,nmax)
    for(n in 1:nmax){
        dtdn_part[n]=-A[n]*lambda[n]*j1(lambda[n],zeta)*exp(-(lambda[n])^2*delta_tau)
}
dtdn_tr=sum(dtdn_part)
dtdn_ss=-1/3*Q*zeta
dtdn=dtdn_tr+dtdn_ss
return(dtdn)
}

solid_profile=function(zeta,delta_tau,theta){
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

liquid_profile=function(zeta,delta_tau,theta){
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
		intagrand[i]=(1-theta_ss[i])*j0(lambda[n],eta[i])*(eta[i])^2
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
	theta_solid[i]=1
	theta_liquid[i]=(xi0^2) - ((eta_solid[i]^2)) +1
}

maxnum=50000
zeta=zeros(1,maxnum)
tau=zeros(1,maxnum)
delta_tau_vec=zeros(1,maxnum)
#saves the initial position and time into vectors
zeta[1]=xi0
tau[1]=0
delta_tau=1/100000
delta_tau_vec[1]=delta_tau
zeta_old=zeta[1]

interface=function(tau,zeta,St){
with(as.list(c(zeta, St)), {
a=solid_temp(tau,zeta,zeta_old,theta_solid)
b=liquid_temp(tau,zeta,theta_liquid)
zetadot=St*(a-b)
return(list(zetadot))
})
}
if(case==1){
T0=0
TF=10
dT=0.001
Ts=seq(from=T0, to=TF, by=dT)
L=length(Ts)
XIs=rep(0,L)
XI=xi0
for(W in seq_along(Ts)){
XIs[W]=XI
dC=dT*St*(1/XI^2*(1/6*Q*XI^2-1/6*Q+1)/(1-1/XI))
XI=XI+dC
}

times=seq(0,10,by=0.01)
solution=ode(0.01,times,interface,St,method="lsoda")
plot(solution,main="",xlab=TeX('$\\tau$'),ylab=TeX('$\\zeta$'))
points(Ts,XIs,type="l",col="red")
solutiont=numeric(length(solution[,1]))
solutionz=numeric(length(solution[,1]))
start=1
for(i in 1:length(solution[,1])){
	solutiont[i]=solution[i,1]
	solutionz[i]=solution[i,2]
	if(i<6){
		start=start+4
	}
	else{
		start=start+49
	}

	print(i)
	print(start)
}
spherical.data=data.frame(solutiont,solutionz)
spherical.data2=data.frame(Ts,XIs)
setwd("C:/Users/Sidney Williams/Pictures/Spherical/Spherical Data")
x=c("Spherical_Data_Melting_St=",St,".txt")
y=c("Spherical_Data_Melting_Quasi_St=",St,".txt")
myfile=toString(x)
myfile2=toString(y)
write.csv(spherical.data,myfile)
write.csv(spherical.data2,myfile2)
}

if(case==2){
times=seq(0,10,by=0.001)
solution=ode(0.01,times,interface,St,method="lsoda")

N=length(solution[,2])-1
eta_temp_liq=linspace(0,solution[1,2],50)
eta_temp_sol=linspace(solution[1,2],1,50)
eta_temp_liq2=linspace(0,solution[50,2],50)
eta_temp_sol2=linspace(solution[50,2],1,50)
eta_temp_liq3=linspace(0,solution[N,2],50)
eta_temp_sol3=linspace(solution[N,2],1,50)
eta_temp_liq4=linspace(0,solution[200,2],50)
eta_temp_sol4=linspace(solution[200,2],1,50)
eta_temp_liq5=linspace(0,solution[500,2],50)
eta_temp_sol5=linspace(solution[500,2],1,50)
eta_temp_liq6=linspace(0,solution[1000,2],50)
eta_temp_sol6=linspace(solution[1000,2],1,50)
eta_temp_liq7=linspace(0,solution[2000,2],50)
eta_temp_sol7=linspace(solution[2000,2],1,50)


liquid_profile1=liquid_profile(solution[1,2],0,theta_liquid)

#Time 1 Temp
liquid_profile2=liquid_profile(solution[1,2],solution[1,1],theta_liquid)
solid_profile2=solid_profile(solution[1,2],solution[1,1],theta_solid)

#Time 2 Temp
liquid_profile3=liquid_profile(solution[50,2],solution[50,1],theta_liquid)
solid_profile3=solid_profile(solution[50,2],solution[50,1],theta_solid)

#Time 3 Temp
liquid_profile5=liquid_profile(solution[200,2],solution[200,1],theta_liquid)
solid_profile5=solid_profile(solution[200,2],solution[200,1],theta_solid)

#Time 4 Temp
liquid_profile6=liquid_profile(solution[500,2],solution[500,1],theta_liquid)
solid_profile6=solid_profile(solution[500,2],solution[500,1],theta_solid)

#Time 5 Temp
liquid_profile7=liquid_profile(solution[1000,2],solution[1000,1],theta_liquid)
solid_profile7=solid_profile(solution[1000,2],solution[1000,1],theta_solid)

#Time 6 Temp
liquid_profile8=liquid_profile(solution[2000,2],solution[2000,1],theta_liquid)
solid_profile8=solid_profile(solution[2000,2],solution[2000,1],theta_solid)

#Time 7 Temp
liquid_profile4=liquid_profile(solution[N,2],solution[N,1],theta_liquid)
solid_profile4=solid_profile(solution[N,2],solution[N,1],theta_solid)

if(St==0.1){
eta_temp_liq=linspace(0,solution[5,2],50)
eta_temp_sol=linspace(solution[5,2],1,50)
eta_temp_liq2=linspace(0,solution[20,2],50)
eta_temp_sol2=linspace(solution[20,2],1,50)
eta_temp_liq3=linspace(0,solution[N,2],50)
eta_temp_sol3=linspace(solution[N,2],1,50)
eta_temp_liq4=linspace(0,solution[50,2],50)
eta_temp_sol4=linspace(solution[50,2],1,50)
eta_temp_liq5=linspace(0,solution[100,2],50)
eta_temp_sol5=linspace(solution[100,2],1,50)
eta_temp_liq6=linspace(0,solution[150,2],50)
eta_temp_sol6=linspace(solution[150,2],1,50)
eta_temp_liq7=linspace(0,solution[300,2],50)
eta_temp_sol7=linspace(solution[300,2],1,50)

#Time 1 Temp
liquid_profile2=liquid_profile(solution[5,2],solution[5,1],theta_liquid)
solid_profile2=solid_profile(solution[5,2],solution[5,1],theta_solid)

#Time 2 Temp
liquid_profile3=liquid_profile(solution[20,2],solution[20,1],theta_liquid)
solid_profile3=solid_profile(solution[20,2],solution[20,1],theta_solid)

#Time 3 Temp
liquid_profile5=liquid_profile(solution[50,2],solution[50,1],theta_liquid)
solid_profile5=solid_profile(solution[50,2],solution[50,1],theta_solid)

#Time 4 Temp
liquid_profile6=liquid_profile(solution[100,2],solution[100,1],theta_liquid)
solid_profile6=solid_profile(solution[100,2],solution[100,1],theta_solid)

#Time 5 Temp
liquid_profile7=liquid_profile(solution[150,2],solution[150,1],theta_liquid)
solid_profile7=solid_profile(solution[150,2],solution[150,1],theta_solid)

#Time 6 Temp
liquid_profile8=liquid_profile(solution[300,2],solution[300,1],theta_liquid)
solid_profile8=solid_profile(solution[300,2],solution[300,1],theta_solid)

#Time 7 Temp
liquid_profile4=liquid_profile(solution[N,2],solution[N,1],theta_liquid)
solid_profile4=solid_profile(solution[N,2],solution[N,1],theta_solid)
}

LWD=3

par(mar=c(5,5,4,2)+.1)
plot(eta_temp_liq4,liquid_profile5,type="l",lty=2,col="blue",ylim=c(0,1.3),xlim=c(0,1),xlab=TeX('$\\eta$-nondimensional radius'),ylab=TeX('$\\theta_{liq},\\theta_{sol}$-nondimensional temperature'),lwd=LWD,cex.lab=1.6,cex.axis=1.3)
#points(eta_temp_sol,solid_profile2,type="l",lty=2,col="blue")
points(eta_temp_liq2,liquid_profile3,type="l",lty=2,col="blue",lwd=LWD)
points(eta_temp_sol2,solid_profile3,type="l", col="red",lwd=LWD)
#points(eta_temp_liq3,liquid_profile4,type="l",col="red",lwd=LWD)
#points(eta_temp_sol3,solid_profile4,type="l",lty=2, col="blue",lwd=LWD)
#points(eta_temp_liq4,liquid_profile5,type="l",col="red")
points(eta_temp_sol4,solid_profile5,type="l",col="red",lwd=LWD)
#points(eta_temp_liq5,liquid_profile6,type="l",col="red")
#points(eta_temp_sol5,solid_profile6,type="l",lty=2,col="blue")
#points(eta_temp_liq6,liquid_profile7,type="l",col="red")
#points(eta_temp_sol6,solid_profile7,type="l",lty=2,col="blue")
points(eta_temp_liq7,liquid_profile8,type="l",lty=2,col="blue",lwd=LWD)
points(eta_temp_sol7,solid_profile8,type="l",col="red",lwd=LWD)
#grid(nx=NULL,ny=NULL,col="lightgray")
arrows(x0=0.7,y0=0.8,x1=0.9,y1=1.1)
text(x=0.56,y=1,labels=TeX('$\\tau$'),cex=1.3)
text(x=0.68,y=1,labels="increases",cex=1.3)
legend(0.1, 0.3, legend=c("Liquid Region", "Solid Region"),col=c("blue", "red"), lty=2:1, cex=1,lwd=LWD,pt.cex=1)
}
