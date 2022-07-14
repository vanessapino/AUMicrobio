#Read all the prediction rasters save a stack as an rds file

####Select depending on which computer are you####

# for mario
basepath_ <- 'C:/Users/mfaj1435/OneDrive - The University of Sydney (Staff)/Molec_Eco/'

#for pinto
# basepath_ <- 'C:/Users/vpin5530/OneDrive - The University of Sydney (Staff)/Papers/My_articles/2nd_beta/' 

#####

# get all predictions folders

folders_ <- c('Bacteria/UMAP1',
              'Bacteria/UMAP2',
              'Bacteria/UMAP3',
              'ITS/UMAP1',
              'ITS/UMAP2',
              'ITS/UMAP3')


lapply(folders_, function(folder){
 # browser()
  name_pred <- gsub('/','',folder,perl = T)
  predictions_ <- raster::stack(dir(paste0('R:/PRJ-Biodiversity/Transects/BiodiversityMaps/PredictionsARTEMIS/',folder),
                                    pattern = 'cubistMap', full.names = T))
  raster::writeRaster(predictions_[[1]],'check.tif',overwrite=TRUE)
  uncertainty <- raster::calc(predictions_, sd)
  raster::writeRaster(uncertainty,paste0(basepath_,'/plots_uncertainty/',name_pred,'.tif'))
})
