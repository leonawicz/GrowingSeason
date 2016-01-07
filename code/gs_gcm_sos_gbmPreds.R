setwd("/atlas_scratch/mfleonawicz/projects/GrowingSeason/workspaces")

library(gbm)
library(maptools)
library(raster)
library(data.table)
library(dplyr)
library(tidyr)
library(parallel)

source("../code/gs_functions.R")

#d.gbm1, d.gbm.preds, d.gbm.data
load("appdata_gbm.RData")

files <- list.files(pattern="^doytdd_qmap_.*.RData")
models <- unlist(strsplit( sapply(strsplit(files, "doytdd_qmap_"), "[", 2), "\\.RData"))
yrs.narr <- 1979:2010
yrs.hist <- 1958:2005
yrs.proj <- 2006:2100
load(files[1])

# Alaska ecoregion level two mask
shpDir <- "/atlas_scratch/mfleonawicz/projects/DataExtraction/data/shapefiles"
eco_shp <- shapefile(file.path(shpDir, "AK_ecoregions/akecoregions.shp"))
eco_shp <- spTransform(eco_shp, CRS(projection(doytdd0[[1]])))
eco_shp <- unionSpatialPolygons(eco_shp, eco_shp@data$LEVEL_2)
ecomask <- rasterize(eco_shp, doytdd0[[1]])
regions <- names(eco_shp)

# Space-time predictions
raster_preds <- function(gbm.table, x.rasterlist, regions, ecomask, yrs=NULL){
    cells <- ncell(x.rasterlist[[1]])
    yr.id <- 1:nlayers(x.rasterlist[[1]])
    if(is.null(yrs)) yrs <- yr.id
    x <- lapply(x.rasterlist, function(x) as.numeric(getValues(x)))
    reg.id <- getValues(ecomask)
    reg.id[is.na(reg.id)] <- length(regions) + 1
    
    get_id_cols <- function(data, id, newdata){
        grp <- as.character(groups(data))
        for(i in 1:length(grp)) newdata <- filter_(newdata, .dots=list(paste0(grp[i], "==\'", data[[grp]][1], "\'")))
        newdata[[id]]
    }
    
    d <- data.table(Index=1:cells, Year_ID=rep(yr.id, each=cells), Region=c(regions, NA)[reg.id],
        DOY_TDD05=x[[1]], DOY_TDD10=x[[2]], DOY_TDD15=x[[3]], DOY_TDD20=x[[4]]) %>% group_by(Region)
    d.preds <- gbm.table %>% do(
        Year_ID=get_id_cols(., id="Year_ID", newdata=d),
        Index=get_id_cols(., id="Index", newdata=d),
        Pred=as.integer(round(get_preds(., model=GBM1, newdata=d, n.trees=BI, type.err="cv")))) %>% unnest
    d <- data.table(left_join(d, d.preds, copy=T)) %>% arrange(Year_ID, Index)
    x <- setValues(x.rasterlist[[1]], d$Pred)
    names(x) <- paste0("SOS_", yrs)
    x
}

system.time({
for(i in 1:length(models)){
    model <- models[i]
    load(files[i])
    b.narr <- raster_preds(d.gbm1, doytdd0, regions, ecomask, yrs.narr)
    print("1. NARR")
    b.hist <- raster_preds(d.gbm1, doytdd1, regions, ecomask, yrs.hist)
    print("2. GCM Historical")
    b.hist.qm <- raster_preds(d.gbm1, doytdd1.mapped, regions, ecomask, yrs.hist)
    print("3. GCM Historical Mapped")
    b.rcp60 <- raster_preds(d.gbm1, doytdd.rcp60, regions, ecomask, yrs.proj)
    print("4. GCM RCP 6.0")
    b.rcp60.qm <- raster_preds(d.gbm1, doytdd.rcp60.mapped, regions, ecomask, yrs.proj)
    print("5. GCM RCP 6.0 Mapped")
    b.rcp85 <- raster_preds(d.gbm1, doytdd.rcp85, regions, ecomask, yrs.proj)
    print("6. GCM RCP 8.5")
    b.rcp85.qm <- raster_preds(d.gbm1, doytdd.rcp85.mapped, regions, ecomask, yrs.proj)
    print("7. GCM RCP 8.5 Mapped")
    save(b.narr, b.hist, b.hist.qm, b.rcp60, b.rcp60.qm, b.rcp85, b.rcp85.qm, file=paste0("sos_gbm_preds_", model, ".RData"))
    print(paste("######## Workspace", i, "of 3 saved. ########"))
}
})
