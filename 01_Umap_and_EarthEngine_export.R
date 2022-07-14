##UMAP calculation and export to earth engine#
require(phyloseq)
require(uwot)
require(ggplot2)
require(DESeq2)


Metadata <- readRDS('filtered_NSW_Dataset.rds')
Details <- get_variable(sample_data(Metadata$Bacteria@sam_data))
Details2 <- Metadata$responses

###16s analysis with SV####

Taxonomy <- readRDS('varianceStabilised16sOtuTable.RDS')

OrderSamples <- match(Details$SampleID,sample_names(Taxonomy))

FilteredTaxonomy <-prune_samples(sample_names(Taxonomy)[OrderSamples],Taxonomy)

FilteredTaxonomy <- filter_taxa(FilteredTaxonomy,function(x)sum(x>2)>1,TRUE)

#check
identical(sort(sample_names(FilteredTaxonomy)),sort(Details$SampleID))

Details <- Details[match(Details$SampleID,sample_names(FilteredTaxonomy)),]

#check
rownames(FilteredTaxonomy)[1:5]
rownames(Details)[1:5]

microbes <- FilteredTaxonomy@.Data

set.seed(123)
coordsUMAP <- uwot::tumap(microbes,n_components =3,pca = 50)

Attributes <- Details

Scores <-data.frame(ID=Attributes,
                    UMAP1=coordsUMAP[,1],
                    UMAP2=coordsUMAP[,2],
                    UMAP3=coordsUMAP[,3])



ggplot(Scores,aes(UMAP1,UMAP2,colour=ID.Transect))+
  geom_point()+
  geom_text(aes(UMAP1,UMAP2,label=ID.Site))




Scores[Scores$ID.Transect=='NS'&as.numeric(Scores$ID.Site)%in%20:22,'time_start'] <- lubridate::as_date('2013-05-22')
Scores[Scores$ID.Transect=='NS'&as.numeric(Scores$ID.Site)%in%c(15:19,23:26),'time_start'] <- lubridate::as_date('2013-03-06')
Scores[Scores$ID.Transect=='NS'&as.numeric(Scores$ID.Site)%in%0:8,'time_start'] <- lubridate::as_date('2013-03-20')
Scores[Scores$ID.Transect=='NS'&as.numeric(Scores$ID.Site)%in%9:14,'time_start'] <- lubridate::as_date('2013-04-06')
Scores[Scores$ID.Transect=='EW'&as.numeric(Scores$ID.Site)%in%27:35,'time_start'] <- lubridate::as_date('2014-03-12')
Scores[Scores$ID.Transect=='EW'&as.numeric(Scores$ID.Site)%in%36:48,'time_start'] <- lubridate::as_date('2014-04-03')
colnames(Scores)[10:11] <- c('longitude','latitude')


order_Details2 <- match(Scores$ID.index,gsub('SAMPLE_','',Details2$index))

toExport <- data.frame(Scores,Details2[order_Details2,])

colnames(toExport)
variablesToExport <- colnames(toExport)[c(13,14,15,16,1,5,6,7,8,9,10,11,30,38,39,45)]

toExport <- toExport[,variablesToExport]

View(toExport)

write.csv(toExport,file = 'Results/ExtractForEarthEngine.csv',row.names = F)


