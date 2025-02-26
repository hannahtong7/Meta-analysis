install.packages("metafor",dependencies=TRUE,repos="https://cloud.r-project.org")

###FUNNEL PLOT###

#1. Simulate Data
#Set the slope and the intercept
slope<--1
intercept<-0
#Generate the predictor variable from a random normal distribution
predictor<-rnorm(n=100,mean=10,sd=10)
#Generate a response using the linear predictor
response<-intercept+slope*predictor+rnorm(n=100,mean=0,sd=40)
plot(predictor,response)
#Run with a lm 
model<-lm(response~predictor)
summary(model)

#2. Simulate many datasets 
store<-matrix(nrow=200,ncol=4)
#We need to create somewhere to store our data
for(x in 1:200){
  #we will simulate 200 different datasets 
  
  samplesize<-ceiling(exp(rnorm(1,4.5,1.5)))+3
  #we're using this code to select sample sizes at random from a log normal distribution, so that small sample sizes are common and large sample sizes are rare. And n is always > 3.                  
  
  
  predictor<-rnorm(n=samplesize,mean=10,sd=10)
  response<-intercept+predictor*slope+rnorm(n=samplesize,0,40)
  #predictor and response are just as we used before, except now n is given by samplesize rather than n = 100
  
  model<-lm(response~predictor)
  #the same linear model as we ran before
  
  store[x,]<-c(samplesize,summary(model)$coefficients[2,1:2],summary(model)$coefficients[2,4])
  #here we extract the model outputs we want and store them in our store matrix
  
  
}
store<-as.data.frame(store)
names(store)<-c("n","slope","standard.error","p.value")

#3. Producing and interpreting a funnel plot
#Have alread generated 200 datasets all with the same true slope
#X-axis is effect size and y-axis is the inverse of the sampling variance (e.g. sample size or standard error)
par(mfrow=c(1,2))
#First plot uses sampling size 
plot(store$slope,store$n,xlab="Slope",ylab="Sample size")
#Second plot uses precision 
plot(store$slope,(1/store$standard.error),xlab="Slope",ylab="Precision, (1/se)")
#Colour the slope estimates which are significant (P<0.05)
#Indicate the slope that was used for the simulation with a vertical dashed line
sigslope<-which(store$p.value<0.05)
par(mfrow=c(1,2))
plot(store$slope,store$n,xlab="Slope",ylab="Sample size")
points(store$slope[sigslope],store$n[sigslope],pch=16,col="red")
abline(v=slope,lty=2)
plot(store$slope,(1/store$standard.error),xlab="Slope",ylab="Precision, (1/se)")
points(store$slope[sigslope],(1/store$standard.error[sigslope]),pch=16,col="red")
abline(v=slope,lty=2)

###BASIC META-ANALYSIS###

#1. Estimating mean effect size ignoring sampling variance 
#Often the aim of a meta-analysis is to estimate the mean effect size
#Finding the mean effect size
model2<-lm(slope~1,data=store)
#The ~1 tells the model to fit an intercept only.
summary(model2)

#2. Estimating the mean effect size using metafor 
library(metafor)
meta<-rma(yi=slope,sei=standard.error,data=store)
meta
#Below is a summary of the key elements of the output:
#The mean estimate and the standard error of the mean and 95% confidence interval
#tau^2 - an estimate of the variance in true effect sizes among studies
#i^2 statistic - tells us what proportion of the total variance (sampling variance + among study heterogeneity in effect size) is due to among study heterogeneity in effect size, see here for more explanation
#Test for heterogeneity - this provides a test of whether there is significant heterogeneity in effect sizes
#Plot to generate a funnel plot and forest plot 
funnel(meta)
forest(meta,cex.lab=0.8,cex.axis=0.8,addfit=TRUE,shade="zebra")

###A META-ANALYSIS WITH MODERATORS AND RANDOM TERMS###

#1. Simulating a new meta-dataset 
#During the previous simulated data above, used the same slope estimate (effect size) every time 
#Now datasets will be generated such that the slope estimates (effect sizes) vary as a function of another variable.
#For example -  let’s imagine that the slopes we are generating correspond to the effect of temperature on phenology (days/degree C), and that this slope becomes more negative with latitude.
#Adjust the simulation loop and add a new step to generate the latitudinal relationship between temperature and the slope. 
#Below latitude predicts the slope e. species from further north advance timings more in response to temperature
latitude<-runif(100,0,90)
#we will randomly sample a latitude from 0,90 degree North
slope<-0+latitude*-0.1+rnorm(100,0,3)
par(mfrow=c(1,1))
plot(latitude,slope)
#Add in a random effect - slopes vary among 20 species
#Slopes 1-10 will be for species 1, slopes 11-20 species 2 and so on
store2<-matrix(nrow=200,ncol=7)
#We need to create somewhere to store our data. We'll call this one store2 to distinguish it from the previous one. This time we also want to save the latitude and species that the slope estimate comes from. We will aslo save a unique ID for each observation - we can use this later to include a residual random effect

species<-rep(1:20,each=10)
specieseffect<-rep(rnorm(20,0,2),each=10)
#we will use this to generate our 20 species random effects

for(x in 1:200){
  #we will simulate 200 different datasets 
  
  latitude<-runif(1,0,90)
  
  slope<-0+specieseffect[x]+latitude*-0.1+rnorm(1,0,3)
  
  samplesize<-ceiling(exp(rnorm(1,4.5,1.5)))
  #we're using this code to select sample sizes at random from a log normal distribution, so that small sample sizes are common and large sample sizes are rare                    
  
  if(samplesize>3){
    #we included this if function so that we don't run an analyses on datasets that are too small
    predictor<-rnorm(n=samplesize,mean=10,sd=10)
    response<-intercept+predictor*slope+rnorm(n=samplesize,0,40)
    
    model<-lm(response~predictor)
    #the same linear model as we ran before
    
    store2[x,]<-c(samplesize,summary(model)$coefficients[2,1:2],summary(model)$coefficients[2,4],latitude,species[x],x)
    #here we extract the model outputs we want and store them in our store matrix
    
    
  }}
store2<-as.data.frame(store2)
names(store2)<-c("n","slope","standard.error","p.value","latitude","species","ID")

#2. Funnel plot and simple meta-analysis 
plot(store2$slope,(1/store2$standard.error),xlab="Slope",ylab="Precision, (1/se)")
meta2<-rma(yi=slope,sei=standard.error,data=store2)
meta2

#3. A meta-analysis controlling for latitude 
#Meta-analysis models can be run like linear or LMM and include fixed and random effects.
#Fixed effects can be reffered to as moderators - all we need to add to the rma function is the model formula to include latitude as a covariate.
meta3<-rma(yi=slope,sei=standard.error,mods=~latitude,data=store2)
meta3
#Look at the latitude slope, and tau value. 
#These should be similar to the latitudinal slope we simulated and the residual standard error in slopes that we simulated, respectively.

#4. A meta-analysis with random terms 
#Use the rma.mv function and can add a random term by adding the argument term random=~1|yourterm
#Included the slope variance as the square of the standard error.
#Also included an observation level random effect (to estimate residual variance)
store2$se2<-store2$standard.error^2
store3<-store2[-which(is.na(store2$slope)==TRUE),]
#this function won't run with NAs, so we remove those rows
meta4<-rma.mv(yi=slope,V=se2,mods=~latitude,random=~1|species/ID,data=store3)
meta4

###CONFRONTING A REAL DATASET###
birdbroods<-read.csv("/Users/hannahtong/Documents/temp/meta/birdbroods.csv",sep=",",header=TRUE)
plot(birdbroods$slope,(1/birdbroods$slope.SE),xlab="Slope",ylab="Precision, (1/se)")
birdbroods$se2<-birdbroods$slope.SE^2
meta5<-rma.mv(yi=slope,V=se2,random=~1|Species/id.pop,data=birdbroods)
meta5
forest(meta5,cex.lab=0.8,cex.axis=0.8,addfit=TRUE,shade="zebra",order="obs")
#Here the effect sizes are ordered 



