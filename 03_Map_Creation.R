### Script that takes covariates rasterstack and earth engine extracted covariates along with umap calculated values#
#and creates 300 bootstraped modelling iterations using Cubist algorithm and saves performance and uncertainty##


library(raster)
library(Cubist)
library(sp)
library(snow)
library(rgdal)
library(clhs)

goof <- function(observed,
                 predicted,
                 coefficient=c('R2','concordance','MSE','RMSE','bias','MSEc','RMSEc','RPD','RPIQ'),
                 plot=TRUE,...){
  
  if(any(!coefficient%in%c('R2','concordance','MSE','RMSE','bias','MSEc','RMSEc','RPD','RPIQ'))) stop('Please choose a valid coefficient')
  
  # Coefficient of determination
  rLM <- lm(predicted ~ observed)
  R2 <- as.matrix(summary(rLM)$adj.r.squared)
  
  # Standard error of prediction ^2
  SEP2 <- mean((observed - predicted)^2)
  
  # Standard error of prediction
  SEP <- sqrt(SEP2)
  
  #Bias
  bias <- mean(predicted) - mean(observed)
  
  # residual  variance
  SEP2c <- sum(((predicted - bias - observed)^2) / length(observed))
  SEPc <- sqrt(SEP2c)
  
  # ratio of performance to deviation
  RPD <- sd(observed) / SEP
  
  # Ratio of performance to interquartile distance
  IQ <- c(quantile(observed))[3] - c(quantile(observed))[2]
  RPIQ <- IQ / SEP
  
  # Concordance
  mx <- mean(observed)
  my <- mean(predicted)
  s2x <- var(observed)
  s2y <- var(predicted)
  sxy <- mean((observed-mx) * (predicted-my))
  ccc <- 2 * sxy / (s2x + s2y + (mx - my)^2)
  
  if(plot){
    plot(observed, 
         predicted,
         ylim=c(min(c(observed,predicted)),max(c(observed,predicted))),
         xlim = c(min(c(observed,predicted)),max(c(observed,predicted))),
         asp=1,
         ...)
    
    abline(a = 0, b = 1, col = "brown4")
  }
  coefs_tmp <- data.frame(R2=R2, concordance=ccc, MSE=SEP2, RMSE=SEP, bias=bias, 
                          MSEc=SEP2c,RMSEc=SEPc, RPD=RPD, RPIQ=RPIQ, row.names=NULL)
  
  
  gf <- data.frame(coefs_tmp[,coefficient])
  
  gf
}


path_Covars<-'../../COVARS_NSW/'
target <- 'UMAP1'
type <- 'Bacteria'
iterations=300

#Locations
dat<- read.csv(dir(pattern ='.csv',full.names = T), header=T) # train data

variables <- c('clay','cec','soc','ph','Rainfall_sum','pha','amp','ET','EVI','aspect','hillshade','slope','B1','B2','B3','B4','B5','B6','K','Th','U')

dat <- dat[,c(target,variables)]



#model
beginCluster(20, type= "SOCK")
set.seed(123)

sampling <- clhs(dat,size=ceiling(nrow(dat)*0.75),iter=15000)

training <- dat[sampling,]
validation <-dat[-sampling,]

clhs_samples <- lapply(1:iterations,function(x){
  set.seed(x)
  
  sample(1:nrow(training),size =nrow(training)*.9)
})

uncertainty <- data.frame(matrix(NA,nrow = iterations,ncol=9))
preds <- data.frame(matrix(NA,nrow = iterations,ncol=nrow(validation)))

stack_NSW <- stack(x = dir(path_Covars,pattern = '.tif',full.names = T))

names(stack_NSW) <- c('clay','cec','soc','ph','Rainfall_sum','pha','amp','ET','LE','PET','PLE','ET_QC','EVI','aspect','hillshade','slope','B1','B2','B3','B4','B5','B6','K','Th','U')


stack_NSW <- stack_NSW[[variables]]

for(i in c(1:iterations)){
  
  samples_to_model <- clhs_samples[[i]]  
  
  Model <- cubist(x = training[samples_to_model,variables], y = training[samples_to_model,target], cubistControl(rules = 5,unbiased = T, label = target), committees = 20)
  
  final_Map <- clusterR(stack_NSW, predict, args=list(Model),filename=paste("cubistMap",i,target,type,".tif",sep='_'),format="GTiff",progress="text",overwrite=T)
  
  saveRDS(Model,file = paste('Model',i,type,target,'.RDS',sep='_'))
  
  uncertainty[i,] <- goof(validation[,target],predict(Model,validation),plot = F)
  preds[i,] <- predict(Model,validation)
  
  
  
  
  
}

colnames(uncertainty) <- c('R2','concordance','MSE','RMSE','bias','MSEc','RMSEc','RPD','RPIQ')

saveRDS(uncertainty,file = 'Uncertainty_models.RDS')
saveRDS(preds,file = 'preds.RDS')
saveRDS(validation,file = 'obs.RDS')

endCluster()         