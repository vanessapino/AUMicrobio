###Script that calculates performance and variable importance of modelled umap maps##

##########
require(vegan)
require(phyloseq)
require(ggplot2)
require(fossil)
require(naturalsort)
require(reshape2)
require(rgdal)
require(plyr)
require(raster)
load('Themes/Vanessa_theme.RData')


#####
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
#####

type <- 'Bacteria'
index <- 'UMAP1'
folder <- 'UMAP1'


getUncertainty <- function(type,index){
  readRDS(paste0('PredictionsARTEMIS/',type,'/',folder,'/Uncertainty_models.RDS'))
}
selectPreds <- function(Uncertainty){
  UP <- quantile(Uncertainty$R2,probs=c(0.75),names = FALSE,na.rm=T)
  DOWN <- quantile(Uncertainty$R2,probs=c(0.25),names = FALSE,na.rm=T)
  return(which(Uncertainty$R2>DOWN&Uncertainty$R2<UP))
}

getObs <- function(type,index){
  tmp <- readRDS(paste0('PredictionsARTEMIS/',type,'/',folder,'/obs.RDS'))
  tmp[,index]
}
getMeanPred <- function(type,index,selectModel){
  tmp <- readRDS(paste0('PredictionsARTEMIS/',type,'/',folder,'/preds.RDS'))
  tmp <- tmp[selectModel,]
  colMeans(tmp)
}

getMeanPredTraining <- function(type,index,selectModel){
  require(Cubist)
  tmp <- naturalsort(dir(paste0('PredictionsARTEMIS/',type,'/',folder),pattern = 'Model',full.names = T))
  tmp <- tmp[selectModel]
  samplesData <- read.csv(dir(paste0('PredictionsARTEMIS/',type,'/',folder),pattern = 'Extract',full.names = T))
  Inputvariables <- c('clay','cec','soc','ph','Rainfall_sum','pha','amp','ET','EVI',
                      'aspect','hillshade','slope',
                      'B1','B2','B3','B4','B5','B6','K','Th','U')
  samplesData <- samplesData[,Inputvariables]
  
  Preds_all_models <- do.call(cbind,lapply(tmp,function(x){
    tmp_model <- readRDS(x)
    tmp_pred <- predict(tmp_model,samplesData)
  }))
  return(rowMeans(Preds_all_models))
}
getObsTraining <- function(type,index){
  samplesData <- read.csv(dir(paste0('PredictionsARTEMIS/',type,'/',folder),pattern = 'Extract',full.names = T))
  return(samplesData[,index])
}
createStack <- function(type,index,selectModel){
  tmp_stack <- stack(naturalsort(dir(paste0('PredictionsARTEMIS/',type,'/',folder),pattern = 'cubistMap',full.names = T)))
  tmp_stack <- tmp_stack[[selectModel]]
  return(tmp_stack)
}

#This functions run on Calc() function together with the Stack
getUpPred <- function(x) quantile(x,probs=0.95,type=7,names = FALSE,na.rm = TRUE)
getLoPred <- function(x) quantile(x,probs=0.05,type=7,names = FALSE,na.rm = TRUE)
getSD <- function(x) sd(x,na.rm = TRUE)


getModelsConditionsCount <- function(type,index,selectModel){
  tmp <- naturalsort(dir(paste0('PredictionsARTEMIS/',type,'/',folder),pattern = 'Model',full.names = T))
  tmp <- tmp[selectModel]
  variables_all_models <- do.call(cbind,lapply(tmp,function(x){
    tmp_model <- readRDS(x)$usage
  }))
  
  COUNTS <-  count(unlist(sapply(1:nrow(variables_all_models),function(x) {
    var_names <- as.character(variables_all_models[x,c(3,6,3:length(tmp)*3)])
    count_per_var <- as.numeric(variables_all_models[x,c(3,6,3:length(tmp)*3)-2])
    vector_freq_vars <- rep(var_names,count_per_var)
    return(vector_freq_vars)
  })
  )
  )
}
getModelsTipsCount <- function(type,index,selectModel){
  tmp <- naturalsort(dir(paste0('PredictionsARTEMIS/',type,'/',folder),pattern = 'Model',full.names = T))
  tmp <- tmp[selectModel]
  variables_all_models <- do.call(cbind,lapply(tmp,function(x){
    tmp_model <- readRDS(x)$usage
  }))
  
  COUNTS <-  count(unlist(sapply(1:nrow(variables_all_models),function(x) {
    var_names <- as.character(variables_all_models[x,c(3,6,3:length(tmp)*3)])
    count_per_var <- as.numeric(variables_all_models[x,c(3,6,3:length(tmp)*3)-1])
    vector_freq_vars <- rep(var_names,count_per_var)
    return(vector_freq_vars)
  })
  )
  )
}
createPerformancePlot <- function(type,index,obsTrain,predTrain,obsVal,predVal,coef='concordance'){
  
  dataVal <- data.frame(Observed=obsVal,PredsValidation=predVal)
  dataTrain <- data.frame(Observed=obsTrain,PredsTraining=predTrain)
  
  min_value <- range(c(dataVal$PredsValidation,dataVal$Observed,dataTrain$Observed,dataTrain$PredsTraining))[1]
  max_value <- range(c(dataVal$PredsValidation,dataVal$Observed,dataTrain$Observed,dataTrain$PredsTraining))[2]
  
  perf_val <- round(goof(dataVal$Observed,dataVal$PredsValidation,coefficient = coef),2)
  perf_train <- round(goof(dataTrain$Observed,dataTrain$PredsTraining,coefficient = coef),2)
  
  ggplot(dataTrain,aes(Observed,PredsTraining))+
    geom_point(colour='blue',size=2)+
    geom_point(data = dataVal,aes(Observed,PredsValidation),colour='red',size=3)+
    coord_equal()+
    geom_abline(slope=1,aes(linetype=2))+
    lims(x = c(min_value, max_value), y = c(min_value, max_value))+
    ylab('Predicted')+
    ggtitle(paste0(type,' - ',index,' Model performace; Train ',coef,' = ', perf_train,'; Val ',coef,' = ', perf_val))
}



#Start#
#Get the uncertainty to choose the right models
Uncertainty <- getUncertainty(type,index)
# View(Uncertainty)
# View(Uncertainty[selectPreds(Uncertainty),])

# length(selectPreds(Uncertainty))
# goof(getObs(type,index),getMeanPred(type,index,selectPreds(Uncertainty)))
# goof(getObsTraining(type,index),getMeanPredTraining(type,index,selectPreds(Uncertainty)))

#made a plot of the performance
plot_performance <- createPerformancePlot(type,
                                          index,
                                          obsTrain = getObsTraining(type,index),
                                          predTrain = getMeanPredTraining(type,index,selectPreds(Uncertainty)),
                                          obsVal=getObs(type,index),
                                          predVal = getMeanPred(type,index,selectPreds(Uncertainty)))

#create a stack of rasters from the selected models (based again on their respective undertainty)
Stack <- createStack(type,index,selectPreds(Uncertainty))
meanPreds <- mean(Stack)

# plot(meanPreds)
UncertPreds <- calc(Stack,getSD)
# plot(UncertPreds)

#write rasters
writeRaster(meanPreds,
            filename = paste0('PredictionsARTEMIS/',type,'/FromR/Preds_',folder,'_',index,'.tif'),
            format="GTiff",
            overwrite=TRUE)

writeRaster(UncertPreds,
            filename = paste0('PredictionsARTEMIS/',type,'/FromR/Uncert_',folder,'_',index,'.tif'),
            format="GTiff",
            overwrite=TRUE)



a <- getModelsConditionsCount(type,index,selectPreds(Uncertainty))
b <- getModelsTipsCount(type,index,selectPreds(Uncertainty))

#write model details 
write.csv(a,file =paste0('PredictionsARTEMIS/',type,'/FromR/Splits_',type,'_',folder,'.csv'))
write.csv(b,file =paste0('PredictionsARTEMIS/',type,'/FromR/Tips_',type,'_',folder,'.csv'))


#save performance plots


ggsave(plot_performance,filename = paste0('PredictionsARTEMIS/',type,'/FromR/Performance_',type,'_',folder,'.png'),height = 7,width = 10)
