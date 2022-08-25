

library(dynr)
library(data.table)
source("dynrMi.R")

missing = 3 # 2 = missing at random (MAR); 
            # 3 = not missing at random (NMAR) - missing at the item level;
            # 4 = NMAR - missing at the latent level

mi = 4 # 1 = LD; 2 = MI-MV; 3 = PMI-MV; 4 = MI-FS
T = 100 # number of time points
n = 30 # number of subjects
r1 = 10 # seed

# dynamic model:
#eta1(t)=a1*eta1(t-1)+b1*eta2(t-1)+c1*x1+d1*x2+u1(t)
#eta2(t)=b2*eta1(t-1)+a2*eta2(t-1)+c2*x1+d2*x2+u2(t)
# measurement model:
#y1(t)=int1+eta1(t)+e1(t)
#y2(t)=int2+lamda21*eta1(t)+e2(t)
#y3(t)=int3+lambda31*eta1(t)+e3(t)
#y4(t)=int4+eta2(t)+e4(t)
#y5(t)=int5+lambda52*eta2(t)+e5(t)
#y6(t)=int6+lambda62*eta2(t)+e6(t)


set.seed(r1)
nt = T
n=n 
ne=2 # number of states
nx=2 # number of fixed regressors
nti=nt+100 #for dropping the first 100 observations


# simualte auxiliary variables: z1,z2
z1=runif(n*nti,-3,3) 
z2=runif(n*nti,-3,3) 
z1l<-matrix(z1,nrow=n,byrow=TRUE)
z2l<-matrix(z2,nrow=n,byrow=TRUE)

# simulate covariates: x1,x2

# binary: x1
theta=rep(NA,(nti*n))
#theta[1]=rbinom(1,1,.5)
theta[1]=rnorm(1,0,1)
x1=rep(0,(nti*n))
for (i in 2:(nti*n)){
  theta[i]=0.8*theta[i-1]+0.4*z1[i]
  pi[i]= 1 / (1 + exp(-theta[i]))
  x1[i]= rbinom(1,1,pi[i])
}

# continuous: x2 (mean=0)
x2=rep(0,(nti*n))
for (i in 1:(nti*n)){
  x2[i]=0.6*z2[i]+rnorm(1,0,0.5)  #avoid loggedevents warnings in mice
}

x1l=matrix(x1,nrow=n,byrow=TRUE)
x2l=matrix(x2,nrow=n,byrow=TRUE)


# simulate latent factors: eta1 and eta2
a1=0.5 ; b1=-0.2; 
b2=-0.3 ; a2=.5;

c1=.3 ; c2=-.3;
d1=.5 ; d2=-.4;

int1=3
int2=3
int3=3
int4=3
int5=3
int6=3

lambda21=2
lambda31=1
lambda52=2
lambda62=1

#cholesky of Q
Q=matrix(c(2, 0.5, #.3
           0.5, 6),ne,ne,byrow=T)
Qs=Q
if(sum(diag(Q))>0) Qs=chol(Q)

sd_e1=1
sd_e2=1
sd_e3=1
sd_e4=1
sd_e5=1
sd_e6=1

eta1l=matrix(NA,nrow=n,ncol=nti)
eta2l=matrix(NA,nrow=n,ncol=nti)
y1l=matrix(NA,nrow=n,ncol=nti)
y2l=matrix(NA,nrow=n,ncol=nti)
y3l=matrix(NA,nrow=n,ncol=nti)
y4l=matrix(NA,nrow=n,ncol=nti)
y5l=matrix(NA,nrow=n,ncol=nti)
y6l=matrix(NA,nrow=n,ncol=nti)

for (i in 1:n){
  
  eta1l[i,1]=rnorm(1,0,1)
  eta2l[i,1]=rnorm(1,0,1)
  y1l[i,1]=rnorm(1,int1+eta1l[i,1],sd_e1)
  y2l[i,1]=rnorm(1,int2+lambda21*eta1l[i,1],sd_e2)
  y3l[i,1]=rnorm(1,int3+lambda31*eta1l[i,1],sd_e3)
  y4l[i,1]=rnorm(1,int4+eta2l[i,1],sd_e4)
  y5l[i,1]=rnorm(1,int5+lambda52*eta2l[i,1],sd_e5)
  y6l[i,1]=rnorm(1,int6+lambda62*eta2l[i,1],sd_e6)
  
  for (t in 2:nti){
    utemp=t(rnorm(ne)%*%Qs)
    temp=matrix(c(a1,b1,b2,a2),nrow=ne,byrow=T) %*% 
      matrix(c(eta1l[i,t-1],eta2l[i,t-1]),nrow=ne)+
      matrix(c(c1,d1,c2,d2),nrow=nx,byrow=T) %*%
      matrix(c(x1l[i,t],x2l[i,t]),nrow=nx)+utemp
    eta1l[i,t]=temp[1,1]
    eta2l[i,t]=temp[2,1]    
    y1l[i,t]=rnorm(1,int1+eta1l[i,t],sd_e1)
    y2l[i,t]=rnorm(1,int2+lambda21*eta1l[i,t],sd_e2)
    y3l[i,t]=rnorm(1,int3+lambda31*eta1l[i,t],sd_e3)
    y4l[i,t]=rnorm(1,int4+eta2l[i,t],sd_e4)
    y5l[i,t]=rnorm(1,int5+lambda52*eta2l[i,t],sd_e5)
    y6l[i,t]=rnorm(1,int6+lambda62*eta2l[i,t],sd_e6)
  }
}


# check stationary
#plot(1:nti, eta1l[1,], type='l', ylim = c(min(eta1l),max(eta1l)),main = "eta1")
#for (i in 2:n){lines(1:nti, eta1l[i,], type='l')}
#plot(1:nti, eta2l[1,], type='l', ylim = c(min(eta2l),max(eta2l)),main = "eta2")
#for (i in 2:n){lines(1:nti, eta2l[i,], type='l')}


#drop first 100 time points
z1full=matrix(z1l[,c(101:nti)],nrow=n)
z2full=matrix(z2l[,c(101:nti)],nrow=n)
x1full=matrix(x1l[,c(101:nti)],nrow=n)
x2full=matrix(x2l[,c(101:nti)],nrow=n)
eta1full=matrix(eta1l[,c(101:nti)],nrow=n)
eta2full=matrix(eta2l[,c(101:nti)],nrow=n)
y1full=matrix(y1l[,c(101:nti)],nrow=n)
y2full=matrix(y2l[,c(101:nti)],nrow=n)
y3full=matrix(y3l[,c(101:nti)],nrow=n)
y4full=matrix(y4l[,c(101:nti)],nrow=n)
y5full=matrix(y5l[,c(101:nti)],nrow=n)
y6full=matrix(y6l[,c(101:nti)],nrow=n)

#full data
if(missing==0){
  pfadata=matrix(NA,ncol=14,nrow=0)
  for (i in 1:n){
    temp=cbind(rep(i,nt),seq(1:nt),x1full[i,],x2full[i,],eta1full[i,],eta2full[i,],
               z1full[i,],z2full[i,],y1full[i,],y2full[i,],y3full[i,],y4full[i,],y5full[i,],y6full[i,])
    pfadata=rbind(pfadata,temp)
  }
  colnames(pfadata)=c("ID","Time","x1","x2","eta1","eta2","z1","z2","y1","y2","y3","y4","y5","y6")
  pfadata=as.data.frame(pfadata)
  pfadata[,3]=as.factor(pfadata[,3])
  ddata=pfadata
}



#generate missing data

missing_gen<-function(phi0,nmareta1=0,nmareta2=0,
                      nmarx1=0,marx11=0.6,marx12=0.6,
                      nmarx2=0,marx21=0.6,marx22=0.6,
                      nmary1=0,mary11=0.6,mary12=0.6,
                      nmary2=0,mary21=0.6,mary22=0.6,
                      nmary3=0,mary31=0.6,mary32=0.6,
                      nmary4=0,mary41=0.6,mary42=0.6,
                      nmary5=0,mary51=0.6,mary52=0.6,
                      nmary6=0,mary61=0.6,mary62=0.6){
logitx1=matrix(rep(NA,nt*n),nrow=n)
prx1=matrix(rep(NA,nt*n),nrow=n)
rx1=matrix(rep(NA,nt*n),nrow=n)

logitx2=matrix(rep(NA,nt*n),nrow=n)
prx2=matrix(rep(NA,nt*n),nrow=n)
rx2=matrix(rep(NA,nt*n),nrow=n)

logity1=matrix(rep(NA,nt*n),nrow=n)
pry1=matrix(rep(NA,nt*n),nrow=n)
ry1=matrix(rep(NA,nt*n),nrow=n)

logity2=matrix(rep(NA,nt*n),nrow=n)
pry2=matrix(rep(NA,nt*n),nrow=n)
ry2=matrix(rep(NA,nt*n),nrow=n)

logity3=matrix(rep(NA,nt*n),nrow=n)
pry3=matrix(rep(NA,nt*n),nrow=n)
ry3=matrix(rep(NA,nt*n),nrow=n)

logity4=matrix(rep(NA,nt*n),nrow=n)
pry4=matrix(rep(NA,nt*n),nrow=n)
ry4=matrix(rep(NA,nt*n),nrow=n)

logity5=matrix(rep(NA,nt*n),nrow=n)
pry5=matrix(rep(NA,nt*n),nrow=n)
ry5=matrix(rep(NA,nt*n),nrow=n)

logity6=matrix(rep(NA,nt*n),nrow=n)
pry6=matrix(rep(NA,nt*n),nrow=n)
ry6=matrix(rep(NA,nt*n),nrow=n)


for (i in 1:n){
  for (t in 2:nt){ #avoid missingness in the first time point
    logitx1[i,t]=phi0+marx11*z1full[i,t]+marx12*z2full[i,t]+nmarx1*x1full[i,t]
    prx1[i,t]=exp(logitx1[i,t])/(1+exp(logitx1[i,t]))
    rx1[i,t]=rbinom(1,1,prx1[i,t])
    
    logitx2[i,t]=phi0+marx21*z1full[i,t]+marx22*z2full[i,t]+nmarx2*x2full[i,t]
    prx2[i,t]=exp(logitx2[i,t])/(1+exp(logitx2[i,t]))
    rx2[i,t]=rbinom(1,1,prx2[i,t])
    
    logity1[i,t]=phi0+mary11*z1full[i,t]+mary12*z2full[i,t]+
      nmary1*y1full[i,t]+nmareta1*eta1full[i,t]
    pry1[i,t]=exp(logity1[i,t])/(1+exp(logity1[i,t]))
    ry1[i,t]=rbinom(1,1,pry1[i,t])
    
    logity2[i,t]=phi0+mary21*z1full[i,t]+mary22*z2full[i,t]+
      nmary2*y2full[i,t]+nmareta1*eta1full[i,t]
    pry2[i,t]=exp(logity2[i,t])/(1+exp(logity2[i,t]))
    ry2[i,t]=rbinom(1,1,pry2[i,t])
    
    logity3[i,t]=phi0+mary31*z1full[i,t]+mary32*z2full[i,t]+
      nmary3*y3full[i,t]+nmareta1*eta1full[i,t]
    pry3[i,t]=exp(logity3[i,t])/(1+exp(logity3[i,t]))
    ry3[i,t]=rbinom(1,1,pry3[i,t])
    
    logity4[i,t]=phi0+mary41*z1full[i,t]+mary42*z2full[i,t]+
      nmary4*y4full[i,t]+nmareta2*eta2full[i,t]
    pry4[i,t]=exp(logity4[i,t])/(1+exp(logity4[i,t]))
    ry4[i,t]=rbinom(1,1,pry4[i,t])
    
    logity5[i,t]=phi0+mary51*z1full[i,t]+mary52*z2full[i,t]+
      nmary5*y5full[i,t]+nmareta2*eta2full[i,t]
    pry5[i,t]=exp(logity5[i,t])/(1+exp(logity5[i,t]))
    ry5[i,t]=rbinom(1,1,pry5[i,t])
    
    logity6[i,t]=phi0+mary61*z1full[i,t]+mary62*z2full[i,t]+
      nmary6*y6full[i,t]+nmareta2*eta2full[i,t]
    pry6[i,t]=exp(logity6[i,t])/(1+exp(logity6[i,t]))
    ry6[i,t]=rbinom(1,1,pry6[i,t])
    
    
  }#end of t loop
  rx1[i,1]<-0
  rx2[i,1]<-0
  ry1[i,1]<-0
  ry2[i,1]<-0
  ry3[i,1]<-0
  ry4[i,1]<-0
  ry5[i,1]<-0
  ry6[i,1]<-0
}#end of i loop
return(list(rx1=rx1,rx2=rx2,ry1=ry1,ry2=ry2,ry3=ry3,ry4=ry4,ry5=ry5,ry6=ry6))
}


##MAR missing data
if (missing ==2){
  
  mar<-missing_gen(phi0=-1.1)
  
  print(sum(mar[[1]])/(n*nt))
  print(sum(mar[[2]])/(n*nt))
  print(sum(mar[[3]])/(n*nt))
  print(sum(mar[[4]])/(n*nt))
  print(sum(mar[[5]])/(n*nt))
  print(sum(mar[[6]])/(n*nt))
  print(sum(mar[[7]])/(n*nt))
  print(sum(mar[[8]])/(n*nt))
  
  marmissx1<-x1full
  marmissx1[mar[[1]]==1]<-NA
  marmissx2<-x2full
  marmissx2[mar[[2]]==1]<-NA
  marmissy1<-y1full
  marmissy1[mar[[3]]==1]<-NA
  marmissy2<-y2full
  marmissy2[mar[[4]]==1]<-NA
  marmissy3<-y3full
  marmissy3[mar[[5]]==1]<-NA
  marmissy4<-y4full
  marmissy4[mar[[6]]==1]<-NA
  marmissy5<-y5full
  marmissy5[mar[[7]]==1]<-NA
  marmissy6<-y6full
  marmissy6[mar[[8]]==1]<-NA
  
  #create dataset for dynr
  ddata=matrix(NA,ncol=14,nrow=0)
  
  #MAR data
  for (i in 1:n){
    temp=cbind(rep(i,nt),seq(1:nt),marmissx1[i,],marmissx2[i,],eta1full[i,],eta2full[i,],
               z1full[i,],z2full[i,],marmissy1[i,],marmissy2[i,],marmissy3[i,],
               marmissy4[i,],marmissy5[i,],marmissy6[i,])
    ddata=rbind(ddata,temp)
  }
  
  colnames(ddata)=c("ID","Time","x1","x2","eta1","eta2","z1","z2","y1","y2","y3","y4","y5","y6")
  ddata=as.data.frame(ddata)
  
  ddata[,3]=as.factor(ddata[,3]) 
  #save(ddata, file="mardata.Rdata")
}


##NMAR missing data
if(missing ==3){
  nmar<-missing_gen(phi0=-0.7,
                      nmarx1=-0.8)
 
  print(sum(nmar[[1]])/(n*nt))
  nmarmissx1<-x1full
  nmarmissx1[nmar[[1]]==1]<-NA
  
  nmar<-missing_gen(phi0=-1,
                    nmarx2=-0.8)
 
  print(sum(nmar[[2]])/(n*nt))
  nmarmissx2<-x2full
  nmarmissx2[nmar[[2]]==1]<-NA

  nmar<-missing_gen(phi0=1,
                      nmary1=-0.6,
                      nmary2=-0.6,
                      nmary3=-0.6)
  
  print(sum(nmar[[3]])/(n*nt))
  print(sum(nmar[[4]])/(n*nt))
  print(sum(nmar[[5]])/(n*nt))
  

  nmarmissy1<-y1full
  nmarmissy1[nmar[[3]]==1]<-NA
  nmarmissy2<-y2full
  nmarmissy2[nmar[[4]]==1]<-NA
  nmarmissy3<-y3full
  nmarmissy3[nmar[[5]]==1]<-NA

  nmar<-missing_gen(phi0=-3,
                      nmary4=0.6,
                      nmary5=0.6,
                      nmary6=0.6)

  print(sum(nmar[[6]])/(n*nt))
  print(sum(nmar[[7]])/(n*nt))
  print(sum(nmar[[8]])/(n*nt))
  
  nmarmissy4<-y4full
  nmarmissy4[nmar[[6]]==1]<-NA
  nmarmissy5<-y5full
  nmarmissy5[nmar[[7]]==1]<-NA
  nmarmissy6<-y6full
  nmarmissy6[nmar[[8]]==1]<-NA

  #create dataset for dynr
  ddata=matrix(NA,ncol=14,nrow=0)
  
  #NMAR data
  for (i in 1:n){
    temp=cbind(rep(i,nt),seq(1:nt),nmarmissx1[i,],nmarmissx2[i,],eta1full[i,],eta2full[i,],
               z1full[i,],z2full[i,],nmarmissy1[i,],nmarmissy2[i,],nmarmissy3[i,],
               nmarmissy4[i,],nmarmissy5[i,],nmarmissy6[i,])
    ddata=rbind(ddata,temp)
  }
  
  colnames(ddata)=c("ID","Time","x1","x2","eta1","eta2","z1","z2","y1","y2","y3","y4","y5","y6")
  ddata=as.data.frame(ddata)
  
  ddata[,3]=as.factor(ddata[,3]) 
  #save(ddata, file = "nmardata.Rdata")
}

if(missing ==4){
  
  nmar<-missing_gen(phi0=-0.7,
                      nmarx1=-0.8)
  print(sum(nmar[[1]])/(n*nt))
  nmarmissx1<-x1full
  nmarmissx1[nmar[[1]]==1]<-NA
  
  
  nmar<-missing_gen(phi0=-1,
                    nmarx2=-0.8)
  print(sum(nmar[[2]])/(n*nt))
  nmarmissx2<-x2full
  nmarmissx2[nmar[[2]]==1]<-NA
  
 
  nmar<-missing_gen(phi0=-1,nmareta1=-0.6)
  print(sum(nmar[[3]])/(n*nt))
  print(sum(nmar[[4]])/(n*nt))
  print(sum(nmar[[5]])/(n*nt))
  

  nmarmissy1<-y1full
  nmarmissy1[nmar[[3]]==1]<-NA
  nmarmissy2<-y2full
  nmarmissy2[nmar[[4]]==1]<-NA
  nmarmissy3<-y3full
  nmarmissy3[nmar[[5]]==1]<-NA


  nmar<-missing_gen(phi0=-1,nmareta2=0.6)
  print(sum(nmar[[6]])/(n*nt))
  print(sum(nmar[[7]])/(n*nt))
  print(sum(nmar[[8]])/(n*nt))
  
  nmarmissy4<-y4full
  nmarmissy4[nmar[[6]]==1]<-NA
  nmarmissy5<-y5full
  nmarmissy5[nmar[[7]]==1]<-NA
  nmarmissy6<-y6full
  nmarmissy6[nmar[[8]]==1]<-NA

  #create dataset for dynr
  ddata=matrix(NA,ncol=14,nrow=0)
  
  #NMAR data
  for (i in 1:n){
    temp=cbind(rep(i,nt),seq(1:nt),nmarmissx1[i,],nmarmissx2[i,],eta1full[i,],eta2full[i,],
               z1full[i,],z2full[i,],nmarmissy1[i,],nmarmissy2[i,],nmarmissy3[i,],
               nmarmissy4[i,],nmarmissy5[i,],nmarmissy6[i,])
    ddata=rbind(ddata,temp)
  }
  
  colnames(ddata)=c("ID","Time","x1","x2","eta1","eta2","z1","z2","y1","y2","y3","y4","y5","y6")
  ddata=as.data.frame(ddata)
  
  ddata[,3]=as.factor(ddata[,3]) 
  #save(ddata, file = "nmardata.Rdata")
}

# create missing indicator
ddata$ry1 = ifelse(is.na(ddata$y1),1,0)
ddata$ry2 = ifelse(is.na(ddata$y2),1,0)
ddata$ry3 = ifelse(is.na(ddata$y3),1,0)
ddata$ry4 = ifelse(is.na(ddata$y4),1,0)
ddata$ry5 = ifelse(is.na(ddata$y5),1,0)
ddata$ry6 = ifelse(is.na(ddata$y6),1,0)
ddata$rx1 = ifelse(is.na(ddata$x1),1,0)
ddata$rx2 = ifelse(is.na(ddata$x2),1,0)

###############dynr
## prepare data
rawdata <- dynr.data(ddata, id="ID", time="Time",
                     observed=paste0("y",1:6), covariates=c("x1","x2"))

#Define the dynamic model
dynamics <- prep.matrixDynamics(
  values.dyn=matrix(c(.5, -.2, -.3, .5), ncol=2,byrow=TRUE),
  params.dyn=matrix(c('a1', 'b1', 'b2', 'a2'), ncol=2,byrow=TRUE), 
  values.exo = matrix(c(.3,-.3,
                        .5,-.4), ncol = 2, byrow = FALSE), 
  params.exo = matrix(c('c1','c2',
                        'd1','d2'), ncol = 2, byrow = FALSE),
  covariates = c('x1','x2'),
  isContinuousTime=FALSE)

meas <- prep.measurement(
  values.load=matrix(c(1,0,
                       2,0,
                       1,0,
                       0,1,
                       0,2,
                       0,1),ncol=2,byrow=TRUE),
  params.load=matrix(c('fixed',0,
                       'lambda21',0,
                       'lambda31',0,
                       0,'fixed',
                       0,'lambda52',
                       0,'lambda62'),
                     ncol=2,byrow=TRUE),  
  values.int = matrix(rep(3,6)),
  params.int = matrix(paste0('int',1:6)), 
  state.names=c("eta1","eta2"), 
  obs.names=paste0('y',1:6) 
)

#Note that in dynr, prep.initial sets the structure of E(eta(1|0)) and Cov(eta(1|0))
#Here, initial condition covariance matrix is fixed to a diagonal matrix of 2s. 
#Could also be freely estimated with #multiple-subject data.
#Iinitial means are fixed to a vector of zeros.
initial <- prep.initial(
  values.inistate=c(0.5, -0.5),
  params.inistate=c('fixed', 'fixed'),
  values.inicov=matrix(c(3,-1,-1,9),ncol=2), #model implied covariance matrix
  params.inicov=matrix(c('fixed','fixed','fixed','fixed'),ncol=2))
#values.inicov=matrix(c(2,0,0,2),ncol=2), 
#params.inicov=matrix(c('c11','c12','c12','c22'),ncol=2))

#T=matrix(c(0.5,-0.2,-0.3,0.5),nrow=2,byrow=TRUE)
#vecsigma = c(2,0.5,0.5,6)
#I = diag(dim(Q)[1]^2)
#P = matrix(solve(I-kronecker(T,T))%*%vecsigma,nrow=2)


#Process and measurement noise covariance matrices
mdcov <- prep.noise(
  values.latent=matrix(c(2,.5,
                         .5,6),ncol=2,byrow=TRUE), 
  params.latent=matrix(c('v11','v12',
                         'v12','v22'),ncol=2,byrow=TRUE), 
  values.observed=diag(rep(1,6),6), 
  params.observed=diag(paste0('var_e',1:6),6)
)


#Put recipes and data together to prepare the full model
model <- dynr.model(dynamics=dynamics, measurement=meas,
                    noise=mdcov, initial=initial, data=rawdata,
                    outfile="PFA.c")
model@compileLib = FALSE


if(missing == 0){ #full data
  res=dynr.cook(model,verbose=FALSE)#,seed=12345)
  if (res$exitflag < 0 | sum(res$'bad.standard.errors') > 0){
    fail = 1
  }else{
    fail = 0
  }
  result=cbind(res$'transformed.parameters',res$'standard.errors',res$'conf.intervals')
  colnames(result)=c("Estimate","Std. Error","ci.lower","ci.upper")
}else{
  if(mi == 1){ # listwise deletion
    ddata2 = ddata[!with(ddata,is.na(x1) | is.na(x2)),]
    ddata3 = NULL
    ID=unique(ddata2$ID)
    for(i in 1:length(ID)){
      tmp=ddata2[ddata2$ID==ID[i],]
      tmp$Time=1:nrow(tmp)
      ddata3=rbind(ddata3,tmp)
    }
    data_listwise <- dynr.data(ddata3, id="ID", time="Time",
                               observed=paste0("y",1:6), covariates=c("x1","x2"))
    model_listwise <- dynr.model(dynamics=dynamics, measurement=meas,
                                 noise=mdcov, initial=initial, data=data_listwise,
                                 outfile="PFAlistwise.c")
    model_listwise@compileLib = FALSE
    res=dynr.cook(model_listwise,verbose=FALSE)#,seed=12345)
    if (res$exitflag < 0 | sum(res$'bad.standard.errors') > 0){
      fail = 1
    }else{
      fail = 0
    }
    result=cbind(res$'transformed.parameters',res$'standard.errors',res$'conf.intervals')
    colnames(result)=c("Estimate","Std. Error","ci.lower","ci.upper")
  }
  if(mi == 2){ # full MI
    res=dynrmi(model, which.aux=c("z1","z2","ry1","ry2","ry3","ry4","ry5","ry6","rx1","rx2"), 
               which.lag=c("x1","x2","y1","y2","y3","y4","y5","y6"), lag=1,
               m=5, iter=30, 
               imp.obs=TRUE, imp.exo=TRUE,
               diag = FALSE,
               verbose=FALSE, seed=12345)
    fail = res$fail
    result=res$estimation.result
  }
  if(mi == 3){ # partial MI
    res=dynrmi(model, which.aux=c("z1","z2","ry1","ry2","ry3","ry4","ry5","ry6","rx1","rx2"), 
                which.lag=c("x1","x2","y1","y2","y3","y4","y5","y6"), lag=1,
                m=5, iter=30, 
                imp.obs=FALSE, imp.exo=TRUE,
                diag = FALSE,
                verbose=FALSE, seed=12345)
    fail = res$fail
    result=res$estimation.result
  }
  if(mi == 4){ 
  #Define the dynamic model
  dynamics_nocov <- prep.matrixDynamics(
      values.dyn=matrix(c(.5, -.2, -.3, .5), ncol=2,byrow=TRUE),
      params.dyn=matrix(c('a1', 'b1', 'b2', 'a2'), ncol=2,byrow=TRUE), 
      isContinuousTime=FALSE)
  #Put recipes and data together to prepare the full model
  model_nocov <- dynr.model(dynamics=dynamics_nocov, measurement=meas,
                              noise=mdcov, initial=initial, data=rawdata,
                              outfile="PFAnocov.c")
  model_nocov@compileLib = FALSE
  # get factor scores
  res_tmp=dynr.cook(model_nocov,verbose=FALSE)
  fsdata=data.frame(t(res_tmp$eta_smooth_final)) #res$eta_filtered
  colnames(fsdata)=c("fs1","fs2")
  newdata = cbind(ddata, fsdata)
  data_addfs <- dynr.data(newdata, id="ID", time="Time",
                          observed=paste0("y",1:6), covariates=c("x1","x2"))
  #Put recipes and data together to prepare the full model
  model_addfs <- dynr.model(dynamics=dynamics, measurement=meas,
                            noise=mdcov, initial=initial, data=data_addfs,
                            outfile="PFAaddfs.c")
  model_addfs@compileLib = FALSE
  if(mi==4){ # MI-FS
  res=dynrmi(model_addfs, which.aux=c("z1","z2","fs1","fs2","ry1","ry2","ry3","ry4","ry5","ry6","rx1","rx2"), 
              which.lag=c("x1","x2","y1","y2","y3","y4","y5","y6","fs1","fs2"), lag=1,
              m=5, iter=30,
              imp.obs=TRUE, imp.exo=TRUE,
              diag = FALSE,
              verbose=FALSE, seed=12345)
  }
  fail = res$fail
  result=res$estimation.result
  }
}


warnings()

if(fail == 0){
write.csv(result,paste0("result_ind_missing",missing,"_mi",mi,"_n",n,"T",T,"_",r1,".csv"))
}

