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
  x <- readAll(brick(files[p]))
  y <- raster(ext=extent(0, 360, -90, 90), res=2.5)
  x <- resample(x, y)
  gc()
  x <- rotate(x) %>% projectRaster(r) %>% crop(r) %>% mask(r)
  gc()
  x[x<0] <- 0
  x <- calc(x, cumsum)
  gc()
  xlist <- purrr::map(pct, ~overlay(x, clim, fun=function(x, y){ ifelse(x < y*.x, 1, 0) }) %>% calc(sum) %>% `+`(1))
  x <- subset(x, 1)
  x <- setValues(x, runif(ncell(x), -0.5, 0.5))
  xlist <- purrr::map(xlist, ~.x + x)
  rm(x)
  gc()
  names(xlist) <- pct
  xlist
}

qmap <- function(i, x0=NULL, x1, deltas, map.by="value", non.negative=TRUE){
  do_adjust <- function(x) (x - min(x) + 0.01)/(max(x) - min(x) + 0.01)
  xi <- x1[[i]]
  if(map.by=="value") a <- t(mapply(function(x0, x1, deltas) if(any(is.na(x1))) rep(NA, length(x1)) else { e <- round(do_adjust(ecdf(x0)(x1))*100); e[e==0] <- 1; deltas[e] },
      as.data.frame(t(getValues(x0[[i]]))), as.data.frame(t(getValues(xi))), as.data.frame(t(getValues(deltas[[i]])))))
  if(map.by=="quantile") a <- t(mapply(function(x1, deltas) if(any(is.na(x1))) rep(NA, length(x1)) else deltas[round(do_adjust(ecdf(x1)(x1))*100)],
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
sos <- readAll(brick("sos_1982_2010.tif"))
r <- calc(sos, mean)
clim <- raster("clim_tdd_1979_2010.tif") %>% projectRaster(r) %>% crop(r) %>% mask(r)
doytdd0 <- mclapply(list.files(pattern="^pct.*.tif$", full=TRUE), function(x, r){
  x <- brick(x) 
  x[x<=1] <- NA
  x %>% projectRaster(r) %>% crop(r) %>% mask(r)
  }, r=r, mc.cores=4)
doytdd0.qtiles <- mclapply(doytdd0, function(x) calc(x, function(x, ...) quantile(x, seq(0, 1, length=100), ...), na.rm=T), mc.cores=4)
n <- length(doytdd0)

# @knitr processing
for(k in 1:length(hist.info$files)){
    doytdd1 <- mclapply(1:length(hist.info$files[[k]]), get_DOYTDDpct, files=hist.info$files[[k]], clim=clim, r=r, mc.cores=24)
    names(doytdd1) <- hist.info$labs[[k]]
    doytdd1 <- map(transpose(doytdd1), stack)
    doytdd1.qtiles <- mclapply(doytdd1, function(x) calc(x, function(x, ...) quantile(x, seq(0, 1, length=100), ...), na.rm=T), mc.cores=n)
    doytdd.deltas <- mclapply(1:n, function(i, x, y) y[[i]] - x[[i]], x=doytdd0.qtiles, y=doytdd1.qtiles, mc.cores=n) # deltas
    doytdd1.mapped <- mclapply(1:n, qmap, x1=doytdd1, deltas=doytdd.deltas, map.by="quantile", mc.cores=n) # historical mapped
    gc()
    doytdd.rcp60 <- mclapply(1:length(rcp60.info$files[[k]]), get_DOYTDDpct, files=rcp60.info$files[[k]], clim=clim, r=r, mc.cores=24)
    names(doytdd.rcp60) <- rcp60.info$labs[[k]]
    doytdd.rcp60 <- lapply(transpose(doytdd.rcp60), stack)
    doytdd.rcp60.mapped <- mclapply(1:length(doytdd0), qmap, x0=doytdd1, x1=doytdd.rcp60, deltas=doytdd.deltas, mc.cores=length(doytdd0)) # rcp60 mapped
    gc()
    doytdd.rcp85 <- mclapply(1:length(rcp85.info$files[[k]]), get_DOYTDDpct, files=rcp85.info$files[[k]], clim=clim, r=r, mc.cores=20)
    names(doytdd.rcp85) <- rcp85.info$labs[[k]]
    doytdd.rcp85 <- lapply(transpose(doytdd.rcp85), stack)
    doytdd.rcp85.mapped <- mclapply(1:length(doytdd0), qmap, x0=doytdd1, x1=doytdd.rcp85, deltas=doytdd.deltas, mc.cores=length(doytdd0)) # rcp85 mapped
    names(doytdd0) <- names(doytdd1.mapped) <- names(doytdd.rcp60.mapped) <- names(doytdd.rcp85.mapped) <- names(doytdd1)
    save(doytdd0, doytdd1, doytdd1.mapped, doytdd.rcp60, doytdd.rcp60.mapped, doytdd.rcp85, doytdd.rcp85.mapped, file=paste0("../workspaces/doytdd_qmap_", hist.info$models[k], ".RData"))
}
