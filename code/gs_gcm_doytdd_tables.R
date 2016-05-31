setwd("/workspace/UA/mfleonawicz/projects/GrowingSeason/workspaces")
pkgs <- list("rgdal", "raster", "maptools", "ggplot2", "data.table", "dplyr", "tidyr", "parallel")
lapply(pkgs, function(x) library(x, character.only=T))

rasterOptions(tmpdir="/atlas_scratch/mfleonawicz/raster_tmp", chunksize=10e10, maxmemory=10e11)
source("../code/gs_functions.R")

files <- list.files(pattern="^doytdd_qmap_.*.RData")
models <- unlist(strsplit( sapply(strsplit(files, "doytdd_qmap_"), "[", 2), "\\.RData"))
yrs.narr <- 1979:2010
yrs.hist <- 1958:2005
yrs.proj <- 2006:2100

sos <- readAll(brick("../data/sos_1982_2010.tif"))
r <- calc(sos, mean)
shpDir <- "/atlas_scratch/mfleonawicz/projects/DataExtraction/data/shapefiles"
eco_shp <- shapefile(file.path(shpDir, "AK_ecoregions/akecoregions.shp")) %>% spTransform(CRS(projection(r)))
eco_shp <- unionSpatialPolygons(eco_shp, eco_shp@data$LEVEL_2)
regions <- names(eco_shp)
dem <- raster("/Data/Base_Data/GIS/GIS_Data/Raster/DEMs/PRISM_2km_DEM/AKCanada_2km_DEM_mosaic.tif") %>% resample(r)
n <- 4

rcp_merge <- function(historical, rcp60, rcp85, model, cells){
  lev <- c("Historical", "RCP 6.0", "RCP 8.5")
  id <- names(rcp60)
  historical <- filter(historical, Cell %in% cells) %>% mutate(RCP=lev[1], Model=model)
  rcp60 <- filter(rcp60, Cell %in% cells) %>% mutate(RCP=lev[2], Model=model)
  rcp85 <- filter(rcp85, Cell %in% cells) %>% mutate(RCP=lev[3], Model=model)
  bind_rows(historical, rcp60, rcp85) %>% mutate(RCP=factor(RCP, levels=lev)) %>% select_(.dots=as.list(c("RCP", "Model", id)))
}

d.tdd <- vector("list", length(files))
system.time({
for(i in seq_along(files)){
  load(files[i])
  if(i==1){
    d.narr <- make_TDD_dt(doytdd0, extractBy=eco_shp, y=sos, keep.y=TRUE, years=yrs.narr, dem=dem, mc.cores=n) %>% filter(!(is.na(DOY_TDD05) | is.na(DOY_TDD10) | is.na(DOY_TDD15) | is.na(DOY_TDD20)))
    d.narr <- mutate(d.narr, RCP="Historical", Model="NARR") %>% select_(.dots=as.list(c("RCP", "Model", names(d.narr))))
  }
  doytdd1.mapped <- make_TDD_dt(doytdd1.mapped, extractBy=eco_shp, y=r, years=yrs.hist, dem=dem, mc.cores=n) %>% filter(!(is.na(DOY_TDD05) | is.na(DOY_TDD10) | is.na(DOY_TDD15) | is.na(DOY_TDD20)))
  doytdd.rcp60.mapped <- make_TDD_dt(doytdd.rcp60.mapped, extractBy=eco_shp, y=r, years=yrs.proj, dem=dem, mc.cores=n) %>% filter(!(is.na(DOY_TDD05) | is.na(DOY_TDD10) | is.na(DOY_TDD15) | is.na(DOY_TDD20)))
  doytdd.rcp85.mapped <- make_TDD_dt(doytdd.rcp85.mapped, extractBy=eco_shp, y=r, years=yrs.proj, dem=dem, mc.cores=n) %>% filter(!(is.na(DOY_TDD05) | is.na(DOY_TDD10) | is.na(DOY_TDD15) | is.na(DOY_TDD20)))
  cells <- intersect(d.narr$Cell, doytdd1.mapped$Cell) %>% intersect(doytdd.rcp60.mapped$Cell) %>% intersect(doytdd.rcp85.mapped$Cell)
  d.narr <- filter(d.narr, Cell %in% cells)
  d.tdd[[i]] <- rcp_merge(doytdd1.mapped, doytdd.rcp60.mapped, doytdd.rcp85.mapped, models[i], cells)
}
})

cells <- sort(unique(unlist(purrr::map(d.tdd, ~.x$Cell))))
d.tdd <- bind_rows(d.narr, bind_rows(d.tdd) %>% filter(Cell %in% cells))
save(d.tdd, file="tdd_table.RData")
