#------------------------------Final Exam-------------------------------
# Based on the R code in my third assignment (a3q2.R)
# R code of Question1 d) e) f) g) and Question2 are altered as well as adjusted as below

#-------------------------------Part 1----------------------------------
#--------------d) EM Algorithm for Gamma Mixture Model------------------

GammaMixEM<-function(Y,clusters,tol,maxits,showits=TRUE){
  # Initialize log-likelihood
  log.likehood<-0
  # Initialize iteration count
  iter.count<-0
  
  # Initialize starting values of theta
  prob<-matrix(c(0.1,0.9),2,1)
  alpha<-matrix(c(0.1,0.5),2,1)
  beta<-matrix(c(0.1,0.5),2,1)
  
  # Initialize probability of cluster membership for each observation i
  n<-nrow(Y)
  wi<-matrix(1,n,clusters)
  
  # Initialize convergence
  converged<-FALSE
  # Show iterations if showits == true  
  if (showits)     
    cat(paste("Iterations of EM:", "\n"))
  while(!converged & iter.count<maxits){ 
    prob.old<-prob
    log.likehood.old<-log.likehood
    wi.old<-wi
    
#------------------------------E-Step----------------------------------
#--------------------Compute responsibilities wi-----------------------
    for (j in 1:clusters){
      wi[,j]<-prob[j,]*dgamma(Y,shape=alpha[j,],scale=beta[j,],log=FALSE) 
    }
    wi<-wi/rowSums(wi)
    
#-------------------------------M-Step----------------------------------
#---------------------Update the value of theta-------------------------

    for (j in 1:clusters){
      C<-log((t(wi[,j])%*%Y)/sum(wi[,j]))-((t(wi[,j])%*%log(Y))/sum(wi[,j]))
      
      # Construct the first-order derivative of log-likehood function about alpha
      f<-log(alpha[j,])-digamma(alpha[j,])-C
      # Construct the second-order derivative of log-likehood function about alpha
      df<-1/alpha[j,]-trigamma(alpha[j,])
      # Update alpha by Newton-Raphson Algorithm
      alpha[j,]<-alpha[j,]-f/df
      
      # Update beta based on first-order derivative of log-likehood function about beta
      beta[j,]<-(t(wi[,j])%*%Y)/(alpha[j,]*sum(wi[,j]))
      # Update p
      prob[j,]<-sum(wi[,j])/n
      # Update log-likehood
      log.likehood[j]<-(-0.5)*sum(wi[,j]*dgamma(Y,shape=alpha[j,],scale=beta[j,]),log=FALSE)
    }
    
    # Update the value of log-likehood for convergence
    log.likehood<-sum(log.likehood)
    
    # Create cluster membership
    cluster<-which(round(wi)==1,arr.ind=TRUE)
    # Order accoring to row rather than cluster
    cluster<-cluster[order(cluster[,1]),2]
    
    ## Compare old to new for convergence
    parameter.old<-c(prob.old,log.likehood.old)
    parameter.new<-c(prob,log.likehood)
    iter.count<-iter.count+1
    if (showits&iter.count==1|iter.count%%5==0)
      cat(paste(format(iter.count),"...","\n",sep=""))
    converged<-min(abs(parameter.old-parameter.new))<tol
  }

  out<-list(prob=prob,alpha=alpha,beta=beta,wi=wi,cluster=cluster,log.likehood=log.likehood)
}




#----------------------------------Part 2----------------------------------------
#----------------e) MLE of Gamma Mixture Model by E-M Algorithm------------------

load("/Users/Desktop/dat.Robj")
Y<-data.matrix(dat)

## Run and get theta
theta<-GammaMixEM(Y,clusters=2,tol=1e-8,maxits=1500,showits=TRUE)

## Numerical Result
alpha1<-theta$alpha[1,]
beta1<-theta$beta[1,]
alpha2<-theta$alpha[2,]
beta2<-theta$beta[2,]
p<-theta$prob[1,]

## Output the M.L.E. of theta
outMLE<-c(alpha1,alpha2,beta1,beta2,p)
outMLE<-matrix(outMLE,1,5)
colnames(outMLE)<-c("alpha1","alpha2","beta1","beta2","p")
rownames(outMLE)<-c("MLE of Theta")

## Graphic Result
hist(Y,breaks=100,freq=FALSE,main="Gamma Mixture Model",xlab='Y',ylab='Density')
x<-seq(0,15,length=1000)
Gamma1<-dgamma(x,shape=alpha1,scale=beta1,log=FALSE)
Gamma2<-dgamma(x,shape=alpha2,scale=beta2,log=FALSE)
GammaMix<-p*Gamma1+(1-p)*Gamma2
lines(x,GammaMix,col="red",lty=1,lwd=3)
legend("topright",legend="Mixture of Gamma",lty=1,col="red",lwd=3,merge=TRUE)




#-----------------------------------Part 3-----------------------------------
#------------------f) Generate random variates from fy-----------------------
GammaMixGen<-function(n,prob,alpha,beta){
  
  # Generate a random number from Bernoulli(p)
  j<-sample.int(length(prob),n,replace=TRUE,prob=prob)
  
  # Generate random variable from fy
  rgamma(n,alpha[j],beta[j])
}



#-----------------------------------Part 4-----------------------------------
#-------------g) Construct a 95% parametric bootstrap CI for p---------------
# Set numbers of Y(s)
S<-1000
# Set numbers of random variables yi(b) of each Y(s)
B<-1000
N<-S
genY<-rep(0,B)
# Get each phat(s) from each Y(s)
phat<-rep(0,S)
z<-rep(0,S)

# Standard Derivation of p in Bernoulli Distribution
sd.p<-sqrt(p*(1-p)/N)

## Apply parametric bootstrap by Gamma Mixture Model
for(s in 1:S){
  # Generate random variates from fy by function in part3
  genY<-GammaMixGen(B,c(p,1-p),c(alpha1,alpha2),c(beta1,beta2))
  
  # Generate thetahat from genY by E-M algorithm
  thetahat<-GammaMixEM(data.matrix(genY),clusters=2,tol=1e-8,maxits=1500)
  # Pick up phat from thetahat
  phat[s]<-thetahat$prob[1,]
  
  # Construct pivot z
  z[s]<-(phat[s]-p)/sqrt(phat[s]*(1-phat[s])/N)
}


## Compute a/2 & 1-a/2 percentile using z
q<-quantile(z,c(0.025,0.975))
ql<-min(q)
qu<-max(q)

# Output the pivot of p
outZ<-c(ql,qu)
outZ<-matrix(outZ,1,2)
colnames(outZ)<-c("2.5%-qL","97.5%-qU")
rownames(outZ)<-"Value of Z"


## Compute 100(1-a)% CI for p
pl<-p-qu*sd.p
pu<-p-ql*sd.p

# Output the CI of p
outPboot<-c(p,sd.p,pl,pu)
outPboot<-matrix(outPboot,1,4)
colnames(outPboot)<-c("Mean-p","Standard Error-p","CI-lower","CI-upper")
rownames(outPboot)<-"Value"




#--------------------------------Part 5----------------------------------
#------------Question 2) Graphicial and Numerical Summaries--------------

# Graphicial Summaries to CEO
plot(x,0*x,main="Risk Aversion of Clients",xlab='Proprietary Metric Value',ylab='Density',col="white",ylim=c(0,0.40))
matlines(x,cbind(Gamma1,Gamma2,GammaMix),lty=c(1,1,2),col=c("blue","red","black"),lwd=c(3,3,4))
legend("topright",legend=c("Risk Averter ( 1st cluster)","Risk Seeker (2nd cluster)","Representative Sample"),lty=c(1,1,2),col=c("blue","red","black"),lwd=c(3,3,4),merge = TRUE)

# Calculate numerical result to CEO
mode1<-(alpha1-1)*beta1
mean1<-alpha1*beta1
var1<-alpha1*(beta1^2)
mode2<-(alpha2-1)*beta2
mean2<-alpha2*beta2
var2<-alpha2*(beta2^2)

# Output the information about Risk Aversion
outRA<-c("Gamma( 1.787, 1.527 )","Gamma(18.450, 0.772)",p,(1-p),mode1,mode2,mean1,mean2,var1,var2)
outRA<-matrix(outRA,2,5)
colnames(outRA)<-c("Model","Probability","Mode","Mean","Variance")
rownames(outRA)<-c("Risk Averter","Risk Seeker")

