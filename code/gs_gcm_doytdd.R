# @knitr setup
lapply(c("raster", "dplyr", "purrr", "parallel"), library, character.only=TRUE)
rasterOptions(tmpdir="/atlas_scratch/mfleonawicz/raster_tmp", chunksize=10e10, maxmemory=10e11)
setwd("/atlas_scratch/mfleonawicz/projects/GrowingSeason/data")

# @knitr functions
prep_files <- function(files){
    labs <- sapply(strsplit(basename(files), "_"), function(x) gsub(".tif", "", paste0(x[c(8,6,5)], collapse="_")))
    models <- sapply(strsplit(labs, "_"), "[", 3)
    list(files=split(files, models), labs=split(labs, models), models=unique(models))
}

get_DOYTDDpct <- function(p, files, clim, r, pct=c(0.05, 0.1, 0.15, 0.2)){
    x <- brick(files[p]) %>% rotate %>% projectRaster(r) %>% crop(r) %>% mask(r)
    gc()
    x[x<0] <- 0
    gc()
    x <- calc(x, cumsum)
    gc()
    xlist <- vector("list", length(pct))
    for(i in 1:length(pct)){
        xlist[[i]] <- overlay(x, clim, fun=function(x, y){ ifelse(x < y*pct[i], 1, 0) }) %>% calc(sum) %>% `+`(1)
        gc()
    }
    rm(x)
    gc()
    names(xlist) <- pct
    xlist
}

qmap <- function(i, x0=NULL, x1, deltas, map.by="value", non.negative=TRUE){
    xi <- x1[[i]]
    if(map.by=="value") a <- t(mapply(function(x0, x1, deltas) if(any(is.na(x1))) rep(NA, length(x1)) else { e <- round(ecdf(x0)(x1)*100); e[e==0] <- 1; deltas[e] },
        as.data.frame(t(getValues(x0[[i]]))), as.data.frame(t(getValues(xi))), as.data.frame(t(getValues(deltas[[i]])))))
    if(map.by=="quantile") a <- t(mapply(function(x1, deltas) if(any(is.na(x1))) rep(NA, length(x1)) else deltas[round(ecdf(x1)(x1)*100)],
        as.data.frame(t(getValues(xi))), as.data.frame(t(getValues(deltas[[i]])))))
    xi2 <- setValues(xi, a)
    xi <- xi - xi2
    if(non.negative) xi[xi<0] <- 0
    xi
}

# @knitr setup2
hist.info <- prep_files(list.files("ar5_daily_tas/historical", pattern="\\.tif$", full=TRUE))
rcp60.info <- prep_files(list.files("ar5_daily_tas/rcp60", pattern="\\.tif$", full=TRUE))
rcp85.info <- prep_files(list.files("ar5_daily_tas/rcp85", pattern="\\.tif$", full=TRUE))
r <- raster("sos_1982_2010.tif")
clim <- raster("clim_tdd_1979_2010.tif") %>% projectRaster(r) %>% crop(r) %>% mask(r)
doytdd0 <- lapply(list.files(pattern="^pct.*.tif$", full=TRUE), function(x, r) brick(x) %>% projectRaster(r) %>% crop(r) %>% mask(r), r=r)
doytdd0.qtiles <- lapply(doytdd0, function(x) calc(x, function(x, ...) quantile(x, seq(0.01, 1, by=0.01), ...), na.rm=T))

# @knitr processing
for(k in 1:length(hist.info$files)){
    doytdd1 <- mclapply(1:length(hist.info$files[[k]]), get_DOYTDDpct, files=hist.info$files[[k]], clim=clim, r=r, mc.cores=16)
    names(doytdd1) <- hist.info$labs[[k]]
    doytdd1 <- lapply(zip_n(doytdd1), stack)
    doytdd1.qtiles <- lapply(doytdd1, function(x) calc(x, function(x, ...) quantile(x, seq(0.01, 1, by=0.01), ...), na.rm=T))
    doytdd.deltas <- lapply(1:length(doytdd0), function(i, x, y) y[[i]] - x[[i]], x=doytdd0.qtiles, y=doytdd1.qtiles) # deltas
    doytdd1.mapped <- mclapply(1:length(doytdd0), qmap, x1=doytdd1, deltas=doytdd.deltas, map.by="quantile", mc.cores=length(doytdd0)) # historical mapped
    gc()
    doytdd.rcp60 <- mclapply(1:length(rcp60.info$files[[k]]), get_DOYTDDpct, files=rcp60.info$files[[k]], clim=clim, r=r, mc.cores=20)
    names(doytdd.rcp60) <- rcp60.info$labs[[k]]
    doytdd.rcp60 <- lapply(zip_n(doytdd.rcp60), stack)
    doytdd.rcp60.mapped <- mclapply(1:length(doytdd0), qmap, x0=doytdd1, x1=doytdd.rcp60, deltas=doytdd.deltas, mc.cores=length(doytdd0)) # rcp60 mapped
    gc()
    doytdd.rcp85 <- mclapply(1:length(rcp85.info$files[[k]]), get_DOYTDDpct, files=rcp85.info$files[[k]], clim=clim, r=r, mc.cores=20)
    names(doytdd.rcp85) <- rcp85.info$labs[[k]]
    doytdd.rcp85 <- lapply(zip_n(doytdd.rcp85), stack)
    doytdd.rcp85.mapped <- mclapply(1:length(doytdd0), qmap, x0=doytdd1, x1=doytdd.rcp85, deltas=doytdd.deltas, mc.cores=length(doytdd0)) # rcp85 mapped
    names(doytdd0) <- names(doytdd1.mapped) <- names(doytdd.rcp60.mapped) <- names(doytdd.rcp85.mapped) <- names(doytdd1)
    save(doytdd0, doytdd1, doytdd1.mapped, doytdd.rcp60, doytdd.rcp60.mapped, doytdd.rcp85, doytdd.rcp85.mapped, file=paste0("../workspaces/doytdd_qmap_", hist.info$models[k], ".RData"))
}
