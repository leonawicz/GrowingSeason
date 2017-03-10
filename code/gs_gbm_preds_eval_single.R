setwd("/atlas_scratch/mfleonawicz/projects/GrowingSeason/workspaces")
pkgs <- list("rgdal", "raster", "maptools", "ggplot2", "data.table", "dplyr", "tidyr", "parallel", "gbm")
lapply(pkgs, function(x) library(x, character.only=T))

rasterOptions(tmpdir="/atlas_scratch/mfleonawicz/raster_tmp", chunksize=10e10, maxmemory=10e11)
load("data.RData") # d, d.stats, d.stats2, d.hm, sos, ecomask, yrs, cbpal
source("../code/gs_functions.R")
sos <- readAll(brick("../data/sos_1982_2010.tif"))
r <- calc(sos, mean)

d <- dcast(d, Region + Year + Obs + x + y + SOS ~ Threshold, value.var="TDD") %>% tbl_df
d <- mutate(d, Cell=cellFromXY(sos, cbind(x,y))) %>% select(-x, -y)
setnames(d, c("Region", "Year", "Obs", "SOS", paste0("DOY_TDD", c("05", 10, 15, 20)), "Cell"))

cells <- (d %>% group_by(Cell) %>% summarise(Region=unique(Region), n=n()) %>% filter(n==29))$Cell
d <- filter(d, Cell %in% cells)

set.seed(564)
use_ak <- FALSE

if(use_ak){
  d_ak <- copy(d)
  d_ak$Region <- "Alaska"
  d <- bind_rows(d, d_ak)
  n.trees <- #100
    #2400
    (0.5)*c(3000, 2500, 4500, 6500, 3000, 1500, 3000, 6000, 1500)
  shrink <- c(0.5, c(0.03808634, 0.01393819, 0.11464891, 0.11070950, 0.04268367, 0.01483978, 0.08393778, 0.20000000, 0.06569468)) # first one is a test for AK
  frac <- c(0.1, c(rep(0.1, 5), .5, .5, .1, .5)) # first one test
  regions <- sort(unique(d$Region))
  predictors <- paste0("DOY_TDD", c(15, 15, 20, 10, "05", 10, "05", 20, "05", "05"))
  file_suffix <- "_withAK"

} else{
  n.trees <- #100
    #2400
    (0.5)*c(3000, 2500, 4500, 6500, 3000, 1500, 3000, 6000, 1500)
  shrink <- #rep(0.2, 9) #
    c(0.029872457, 0.014395137, 0.133199817, 0.119391448, 0.040367520, 0.005136854, 0.022135554, 0.200000000, 0.011665086)
  frac <- c(rep(0.1, 5), .5, .5, .1, .5)
  regions <- sort(unique(d$Region))
  predictors <- paste0("DOY_TDD", c(15, 20, 10, "05", 10, "05", 20, "05", "05"))
  file_suffix <- ""
}

# build gbm models
get_cv_err <- function(.){
  k <- 200
  b <- .$CV[[1]]
  m <- .$GBM1[[1]]
  n <- m$n.trees
  d <- data.table(Trees=1:n, `Training Error`=m$train.error, `CV Error`=m$cv.error)
  be <- d$`CV Error`[b]
  d <- melt(d, id="Trees")
  d <- mutate(d, `Number of Trees`=b, Error=be)
  setnames(d, c("Number of Trees","Type of Error","Error", "Optimal_Trees", "Optimal_Error"))
  nest(d, `Number of Trees`, `Type of Error`, `Error`)
}

gbm_explore <- function(i, data, n.trees, shrinkage, frac, years=sort(unique(data$Year)), by.year=TRUE, agg=FALSE, predictors){
  n <- if(by.year) length(years) else 1
  out <- vector("list", n)
  data <- filter(data, Year %in% years)
  gc()
  if(length(frac)==1) frac <- rep(frac, length(unique(d$Region)))
  for(j in seq_along(out)){
    if(by.year){ d <- filter(data, Year==years[j]) } else { d <- data; rm(data); gc() }
    #d <- group_by(d, Region, Year) %>% sample_frac(frac) %>% group_by(Region)
    d <- d %>% split(.$Region) %>% purrr::map2(frac, ~group_by(.x, Year) %>% sample_frac(.y)) %>% bind_rows %>% data.table %>% group_by(Region)
    if(agg && !by.year) d <- group_by(d, Region, Year) %>% summarise(SOS=mean(SOS),
                                                                     DOY_TDD05=mean(DOY_TDD05),
                                                                     DOY_TDD10=mean(DOY_TDD10),
                                                                     DOY_TDD15=mean(DOY_TDD15),
                                                                     DOY_TDD20=mean(DOY_TDD20)
    ) %>% group_by(Region)
    d.train <- sample_frac(d, 0.5)
    d.test <- setdiff(d, d.train)
    form <- paste("SOS ~", predictors)
    d.gbm <- d.train %>% split(.$Region)
    d.gbm <- d.gbm %>%
      purrr::map2(seq_along(form), ~gbm(as.formula(form[.y]), data=.x,
        distribution="gaussian", bag.fraction=0.5, cv.folds=5, train.fraction=1,
        interaction.depth=1, n.minobsinnode=5, n.trees=n.trees[.y], shrinkage=shrinkage[.y], verbose=FALSE, keep.data=TRUE, n.cores=1))

    d.gbm <- d.train %>% group_by %>% select(Region) %>% distinct(Region) %>% mutate(GBM1=d.gbm) %>% group_by(Region)
    regions <- d.gbm$Region
    #print(paste("regions[1]=='Alaska':", regions[1]=="Alaska"))
    if(regions[1]=="Alaska"){
      regions <- regions[-1]
      #print(regions)
      d.gbm <- d.gbm %>% group_by %>% slice(rep(1, 9)) %>% mutate(Region=regions) %>% group_by(Region)
      d.train <- filter(d.train, Region!="Alaska")
      d.test <- filter(d.test, Region!="Alaska")
    }
    d.tdd <- select(d.train, -Obs, -Year, -Cell) %>% melt(id.vars=c("Region", "SOS"), variable.name="Var", value.name="Val") %>% group_by(Region, Var)

    rm(d.train); gc()
    d.bi <- d.gbm %>% do(BI=get_bi(., model=.$GBM1, plotDir, saveplot=F))
    d.gbm <- data.table(suppressMessages(left_join(d.gbm, d.bi))) %>% group_by(Region)
    d.gbm <- d.gbm %>% mutate(CV=purrr::map(BI, ~.x$CV)) %>% group_by(Region)
    d.ri <- d.gbm %>% do(RI=get_ri(., model=.$GBM1, n.trees=.$BI, plotDir, saveplot=F))
    d.gbm <- data.table(suppressMessages(left_join(d.gbm, d.ri))) %>% group_by(Region)
    d.preds <- d.gbm %>%
      do(Predicted=get_preds(., model=.$GBM1, newdata=d.test, n.trees=.$BI, type.err="cv", grp=as.character(groups(d.gbm)))) %>% group_by(Region)
    d.err <- d.gbm %>% group_by(Region) %>% do(Error=get_cv_err(.))
    d.gbm <- suppressMessages(left_join(d.gbm, d.err)) %>% group_by(Region)
    rm(d.bi, d.ri, d.err); gc()

    d.pd <- d.gbm %>% do(PD=get_pd(., source_data=d.tdd, x="Val", y="SOS", outDir=NULL, model=.$GBM1, vars=1, order.by.ri=TRUE,
                                   suffix=.$Region, saveplot=FALSE, grp=as.character(groups(d.gbm)))) %>% group_by(Region)

    d.gbm <- data.table(left_join(d.gbm, d.pd)) %>% group_by(Region)
    rm(d.tdd); gc()

    d.test$Predicted <- unnest(d.preds)$Predicted
    d.test$Run <- i
    d.test <- d.test %>% select(Region, Year, Run, SOS, Predicted) %>% setnames(c("Region", "Year", "Run", "Observed", "Predicted")) %>%
      melt(id.vars=c("Region", "Year", "Run"), value.name="SOS") %>% data.table %>% setnames(c("Region", "Year", "Run", "Source", "SOS")) %>%
      group_by(Region, Year, Run, Source) %>% summarise(SOS=round(mean(SOS)))
    if(by.year){
      d.coef <- mutate(d.coef, Year=years[j]) %>% select(Region, Year, intercept, slope)
      d.gbm <- mutate(d.gbm, Year=years[j]) %>% select(Region, Year, GBM1, RI, BI, CV, Error, PD)
    }
    out[[j]] <- list(GBM=d.gbm, data=d.test)
    rm(d.gbm, d.test); gc()
    if(i==1 & n > 1) print(j)
  }
  if(n==1){
    out <- out[[1]]
  } else {
    out <- list(GBM=rbindlist(out %>% purrr::map(~.x$GBM)) %>% group_by, data=bind_rows(out %>% purrr::map(~.x$data)), LM=bind_rows(out %>% purrr::map(~.x$LM)))
  }
  out
}

# Run
system.time( dlist <- mclapply(1:32, gbm_explore, d, n.trees=n.trees, shrinkage=shrink, frac=frac,
                               by.year=FALSE, agg=FALSE, predictors=predictors, mc.cores=32) )
gbm_explore(1, d, n.trees=n.trees, shrinkage=shrink, frac=frac, by.year=FALSE, agg=FALSE, predictors=predictors)

# Extract tables of models, CV optimal trees, predictions, corrections, etc.
gbm.out <- lapply(dlist, "[[", 1)
ri.out <-  purrr::map2(gbm.out, seq_along(gbm.out), ~select(.x, Region, RI) %>% mutate(Run=.y) %>% unnest) %>% bind_rows %>% filter(Method=="CV")
cv.out <-  purrr::map(gbm.out, ~select(.x, Region, CV)) %>% bind_rows %>% unnest
err.out <-  purrr::map2(gbm.out, seq_along(gbm.out), ~select(.x, Region, Error) %>% mutate(Run=.y)) %>% bind_rows %>% unnest
pd.out <-  purrr::map2(gbm.out, seq_along(gbm.out), ~group_by(.x) %>% select(PD) %>% mutate(Run=.y)) %>% bind_rows %>% unnest
d.out <- rbindlist(lapply(dlist, "[[", 2))
#lm.out <- lapply(dlist, "[[", 3)
d.out <- group_by(d.out, Region, Year, Source) %>% summarise(SOS=mean(SOS)) %>% bind_rows(filter(d.out))

dir.create("singlePred_outputs", showWarnings=FALSE)
save(ri.out, cv.out, pd.out, d.out, file=paste0("singlePred_outputs/final_gbm_summary_tables", file_suffix, ".RData"))

# use this one
gbm_prediction_maps <- function(newdata, d, r, year.method="match", lm.pars=NULL, output="maps", simplify=TRUE, agg=FALSE, agg.frac=0.05, n.cores=32){
  yrs <- sort(unique(newdata$Year))
  if("Year" %in% names(d[[1]])) gbm.yrs <- sort(unique(d[[1]]$Year))
  grp <- if("Year" %in% names(d[[1]])) list("Region", "Year") else "Region"
  d <- purrr::map(d, ~group_by_(.x, .dots=grp))
  newdata <- if("Year" %in% names(d[[1]])) arrange(newdata, Year, Region, Obs) else arrange(newdata, Region, Year, Obs)
  if(agg) newdata <- group_by(newdata, Region, Year) %>% sample_frac(agg.frac) %>%
    summarise(Cell=Region[1],
              DOY_TDD05=mean(DOY_TDD05),
              DOY_TDD10=mean(DOY_TDD10),
              DOY_TDD15=mean(DOY_TDD15),
              DOY_TDD20=mean(DOY_TDD20)) %>% group_by(Region) # for point models
  if(year.method=="match"){
    d <- mclapply(seq_along(d), function(i, x) x[[i]] %>%
                    do(Predicted=as.integer(round(
                      get_preds(., model=.$GBM1, newdata=newdata, n.trees=.$BI, type.err="cv", grp=as.character(groups(d[[1]])))))) %>%
                    mutate(Run=i), x=d, mc.cores=n.cores)
  } else if(year.method=="random"){
    par_preds <- function(i, x) x[[i]] %>%
      do(Predicted=as.integer(round(
        get_preds(., model=.$GBM1, newdata=newdata[[i]], n.trees=.$BI, type.err="cv", grp=as.character(groups(d[[1]])))))) %>% mutate(Run=i)
    newdata <- purrr::map(seq_along(d), ~mutate(newdata, YearBackup=Year) %>% group_by(YearBackup) %>%
                            mutate(Year=sample(gbm.yrs, 1, replace=TRUE)) %>% group_by)
    d <- mclapply(seq_along(d), par_preds, x=d, mc.cores=n.cores)
    newdata <- mutate(newdata[[1]], Year=YearBackup) %>% select(-YearBackup)
  }
  d <- d %>% purrr::map(~unnest(.x) %>% mutate(Cell=newdata$Cell, Year=as.integer(newdata$Year)) %>%
                          nest(Cell, Predicted) %>% group_by(Region, Year, Run))
  if(!is.null(lm.pars)){
    d <- d %>% purrr::map2(lm.pars, ~suppressMessages(left_join(.x, .y, copy=T)))
    d <- d %>% purrr::map(~unnest(.x) %>% mutate(`Bias corrected`=(Predicted-intercept)/slope) %>%
                            nest(Cell, Predicted, `Bias corrected`) %>% group_by(Region, Year, Run))
  }
  if(output=="table" || agg) return(bind_rows(d)) # for point models
  gc()
  setPred <- function(d, r, values){
    d2 <- d %>% group_by %>% unnest %>% select_(.dots=list(paste0("`", values, "`"), "Cell"))
    r[d2$Cell] <- d2[[values]]
    r
  }
  r <- raster(r)
  r <- setValues(r, NA)
  values <- if(is.null(lm.pars)) "Predicted" else c("Predicted", "Bias corrected")
  b <- vector("list", length(values))
  for(k in seq_along(values)){
    b[[k]] <- mclapply(d, function(x) brick(purrr::map(x %>% split(.$Year), ~setPred(.x, r, values[k]))), mc.cores=n.cores)
    gc()
  }
  names(b) <- values
  if(simplify & length(b)==1) b <- b[[1]]
  b
}

# Run

load("tdd_table.RData")
rcps <- unique(d.tdd$RCP)[-1]
models <- unique(d.tdd$Model)[-1]
pred.maps <- purrr::map(vector("list", length(rcps)), ~vector("list", length(models)) %>% setNames(models)) %>% setNames(rcps)
set.seed(1)

system.time({
  for(i in seq_along(rcps)){
    for(j in seq_along(models)){
      d.gcm <- d.tdd %>% data.table %>% filter(RCP %in% c("Historical", rcps[i]) & Model==models[j]) %>%
        select(-RCP, -Model, -x, -y, -Elev, -SOS) %>% tbl_df
      pred.maps[[i]][[j]] <- gbm_prediction_maps(d.gcm, gbm.out[1:10], r, agg=FALSE, n.cores=10)
      print(j)
    }
  }
})
saveRDS(pred.maps, file=paste0("singlePred_outputs/gbm_pred_rasters", file_suffix, ".rds"))

setwd("/atlas_scratch/mfleonawicz/projects/GrowingSeason/workspaces")
library(raster)
library(maptools)
library(data.table)
library(dplyr)
library(purrr)

load("data.RData") # d, d.stats, d.stats2, d.hm, sos, ecomask, yrs, cbpal
shpDir <- "/atlas_scratch/mfleonawicz/projects/DataExtraction/data/shapefiles"
eco_shp <- shapefile(file.path(shpDir, "AK_ecoregions/akecoregions.shp")) %>% spTransform(CRS(projection(sos)))
eco_shp <- unionSpatialPolygons(eco_shp, eco_shp@data$LEVEL_2)

#pred.maps <- readRDS(paste0("singlePred_outputs/gbm_pred_rasters", file_suffix, ".rds"))
rcps <- c("RCP 6.0", "RCP 8.5")
gcms <- c("GFDL-CM3", "IPSL-CM5A-LR", "MRI-CGCM3")

extract_to_dt <- function(x, y, fun, rcp, gcm, run, ...){
  raster::extract(x, y, fun, ...) %>% t %>% data.table %>% setnames(names(y)) %>% mutate(Run=run) %>%
    melt(id.vars="Run", variable.name="Region", value.name="SOS") %>%
    mutate(Year=rep(as.integer(substr(names(x), 2, 5)), length(names(y))),
           RCP=factor(rcp, levels=c("RCP 6.0", "RCP 8.5")),
           Model=gcm, Source="Projected") %>%
    select(RCP, Model, Region, Year, Source, SOS, Run)
}

d.proj <- vector("list", length(rcps)*length(gcms)*length(pred.maps[[1]][[1]]))
idx <- 1
for(i in seq_along(rcps)){
  for(j in seq_along(gcms)){
    for(k in seq_along(pred.maps[[i]][[j]])){
      d.proj[[idx]] <- extract_to_dt(pred.maps[[i]][[j]][[k]], eco_shp, mean, rcps[i], gcms[j], run=k, na.rm=TRUE)
      idx <- idx + 1
      print(paste("j =", j))
    }
  }
}
d.proj <- bind_rows(d.proj)
saveRDS(d.proj, paste0("singlePred_outputs/sos_projections", file_suffix, ".rds"))
######################################################################################## END
set.seed(1)
n <- 10
outfiles <- paste0("run", 1:n, ".rds")
system.time({
  for(i in seq_along(rcps)){
    for(j in seq_along(models)){
      x <- d.tdd %>% data.table %>% filter(RCP==rcps[i] & Model==models[j]) %>% select(-RCP, -Model, -x, -y, -Elev, -SOS)
      x <- gbm_prediction_maps(x, gbm.out[sample(seq_along(gbm.out), n)], r, year.method="random", n.cores=min(n, 10))
      dir.create(mapDir <- file.path("/atlas_scratch/mfleonawicz/projects/GrowingSeason/workspaces", rcps[i], models[j]), recursive=TRUE, showWarnings=FALSE)
      for(k in 1:n){
        xx <- x[[k]]
        saveRDS(xx, file.path(mapDir, outfiles[k]))
      }
    }
  }
})

#d.tdd.tmp <- d.tdd %>% data.table %>% select(-x, -y, -Elev, -SOS) %>% filter(Model!="NARR")

#f <- function(x) purrr::map(x, ~select(.x, -RCP, -Model) %>% gbm_prediction_maps(gbm.out[1:16], r, n.cores=16))
#system.time({ pred.maps <- d.tdd.tmp %>% split(.$RCP) %>% purrr::map(~.x %>% split(.$Model) %>% f) })

shpDir <- "/atlas_scratch/mfleonawicz/projects/DataExtraction/data/shapefiles"
eco_shp <- shapefile(file.path(shpDir, "AK_ecoregions/akecoregions.shp")) %>% spTransform(CRS(projection(sos)))
eco_shp <- unionSpatialPolygons(eco_shp, eco_shp@data$LEVEL_2)

extract_to_dt <- function(x, y, fun, rcp, gcm, n.cores=32){
  x <- x[[rcp]][[gcm]]
  f <- function(i, ...){
    raster::extract(x[[i]], y, fun, na.rm=TRUE) %>% t %>% data.table %>% setnames(names(y)) %>% mutate(Run=i) %>% melt(id.vars="Run", variable.name="Region", value.name="SOS") %>%
      mutate(Year=as.integer(substr(names(x[[i]]), 2, 5)), RCP=factor(rcp, levels=c("RCP 6.0", "RCP 8.5")), Model=gcm, Source="Projected") %>% select(RCP, Model, Region, Year, Source, SOS, Run)
  }
  mclapply(seq_along(x), f, x=x, y=y, fun=fun, rcp=rcp, gcm=gcm, mc.cores=n.cores) %>% bind_rows
}

# for point models
extract_to_dt <- function(x, y, fun, rcp, gcm, n.cores=32){
  x <- x[[rcp]][[gcm]]
  f <- function(i, ...){
    x <- if(is.data.frame(x[[i]])) rename(unnest(x[[i]]), SOS=Predicted) %>% select(-Cell) else raster::extract(x[[i]], y, fun, na.rm=TRUE) %>% t %>% data.table %>% setnames(names(y)) %>% mutate(Run=i) %>%
      melt(id.vars="Run", variable.name="Region", value.name="SOS") %>% mutate(Year=as.integer(substr(names(x[[i]]), 2, 5)))
    x <- mutate(x, RCP=factor(rcp, levels=c("RCP 6.0", "RCP 8.5")), Model=gcm, Source="Projected") %>% select(RCP, Model, Region, Year, Source, SOS, Run)
    x
  }
  if(is.data.frame(x)) x <- list(x)
  mclapply(seq_along(x), f, x=x, y=y, fun=fun, rcp=rcp, gcm=gcm, mc.cores=n.cores) %>% bind_rows
}

d.hoy <- d.preds %>% select(-SOS) %>% mutate(Source="Predicted HOY") %>% rename(SOS=Predicted)
d.hoy <- d.preds %>% select(-Predicted) %>% distinct %>% mutate(Source="Global observed") %>% bind_rows(d.hoy)
d.ts <- bind_rows(filter(d.out, Source!="Bias corrected"), d.hoy) %>% mutate(Source=factor(Source, levels=c("Observed", "Predicted", "Global observed", "Predicted HOY")))

d.proj1 <- extract_to_dt(pred.maps, eco_shp, mean, "RCP 6.0", "GFDL-CM3")
d.proj2 <- extract_to_dt(pred.maps, eco_shp, mean, "RCP 6.0", "IPSL-CM5A-LR")
d.proj3 <- extract_to_dt(pred.maps, eco_shp, mean, "RCP 6.0", "MRI-CGCM3")
d.proj4 <- extract_to_dt(pred.maps, eco_shp, mean, "RCP 8.5", "GFDL-CM3")
d.proj5 <- extract_to_dt(pred.maps, eco_shp, mean, "RCP 8.5", "IPSL-CM5A-LR")
d.proj6 <- extract_to_dt(pred.maps, eco_shp, mean, "RCP 8.5", "MRI-CGCM3")

d.proj.mean <- d.proj %>% group_by(Region, Year, Source) %>% summarise(SOS=mean(SOS), Run=1)
d.smooth <- filter(d.ts, is.na(Run) & Source=="Predicted") %>% bind_rows(filter(d.proj, Year > 2010)) %>% group_by(Region, Year) %>% summarise(SOS=mean(SOS), Source="Trend", Run=1)

dir.create(plotDir <- file.path("../plots/gbm/models"), recursive=T, showWarnings=F)

# historical predictions over observations time series
clrs <- c("peru", "royalblue", "#8B451320", "#00008B20", "#FF303020")
png(file.path(plotDir, paste0("gbm_TSpreds_byRegion_plusRCP60_01pctSample.png")), width=3200, height=1600, res=200)
ggplot(filter(d.out, !is.na(Run) & Source!="Bias corrected"), aes(x=Year, y=SOS, colour=Source, group=interaction(Source, Run))) +
  scale_color_manual(values=clrs) +
  geom_line(size=1, alpha=0.25) +
  #geom_line(data=filter(d.out, Source=="Bias corrected" & is.na(Run)), colour="#B8860B", size=1, linetype=1) +
  geom_line(data=filter(d.out, Source=="Observed" & is.na(Run)), colour="black", size=2) +
  geom_line(data=filter(d.out, Source=="Predicted" & is.na(Run)), colour="black", size=2) +
  geom_line(data=filter(d.out, Source=="Predicted" & is.na(Run)), colour=clrs[2], size=1) +
  geom_line(data=filter(d.out, Source=="Observed" & is.na(Run)), colour=clrs[1], size=1) +
  geom_line(data=filter(d.out, Source=="Predicted" & is.na(Run)), colour=clrs[2], size=1, linetype=2) +
  geom_line(data=d.proj1 %>% mutate(Source="Projected: GFDL"), size=1) +
  geom_line(data=d.proj2 %>% mutate(Source="Projected: IPSL"), size=1) +
  geom_line(data=d.proj3 %>% mutate(Source="Projected: MRI"), size=1) +
  #geom_smooth(data=d.smooth) +
  theme_bw() + theme(legend.position="bottom") + ggtitle("Observed and modeled start of growing season") +
  #scale_x_continuous(breaks=c(1982,1990,2000,2010)) +
  guides(fill=guide_legend(override.aes=list(alpha=1)), colour=guide_legend(override.aes=list(alpha=1))) +
  facet_wrap(~Region, ncol=3, scales="fixed")
dev.off()





# Prep specific maps for plotting script
pred.maps.gbm.samples <- pred.maps[[1]] %>% purrr::map(~calc(.x, mean)) %>% stack
pred.maps.gbm1 <- pred.maps %>% purrr::map(~.x[[32]])

s1 <- calc(pred.maps.gbm1$Predicted - sos, mean) #stack(calc(pred.maps.gbm1$Predicted - sos, mean), calc(pred.maps.gbm1$`Bias corrected` - sos, mean))
sMean <- pred.maps[[1]] %>% purrr::map(~.x-sos) #list(pred.maps[[1]] %>% purrr::map(~.x-sos), pred.maps[[2]] %>% purrr::map(~.x-sos))
sMean <- Reduce("+", sMean)/length(sMean) #sMean %>% purrr::map(~Reduce("+", .x)/length(.x))
sMean2002 <- subset(sMean, match(2002, yrs)) #sMean %>% purrr::map(~subset(.x, match(2002, yrs))) %>% stack
sMean <- calc(sMean, mean) #sMean %>% purrr::map(~calc(.x, mean)) %>% stack
names(s1) <- names(sMean) <- names(sMean2002) <- "Prediction_deltas" #names(s1) <- names(sMean) <- names(sMean2002) <- c("Prediction_deltas", "Bias_corrected_deltas")
#s1 <- subset(s1, 1)
#sMean <- subset(sMean, 1)
#sMean2002 <- subset(sMean2002, 1)

set.seed(1)
predict_hold_out_year <- function(gbm.out, d, regions, years){
  inner_fun <- function(gbm.out, d, region, year){
    d.test <- d %>% filter(Region==region & Year==year) %>% group_by(Region, Year)
    rep.gbm <- sample(seq_along(gbm.out), length(gbm.out), replace=T)
    d.preds <- purrr::map2(gbm.out, rep.gbm,
                           ~filter(.x, Region==region & Year!=year) %>% sample_n(1) %>% group_by(Region) %>% do(Predicted=get_preds(., model=GBM1, newdata=d.test, n.trees=BI, type.err="cv") %>% mean)) %>%
      bind_rows %>% unnest %>% mutate(Region=region, Year=year)
    d.test <- suppressMessages(summarise(d.test, SOS=mean(SOS)) %>% left_join(d.preds, copy=T))
    print(paste0(region, ": ", year))
    d.test
  }
  purrr::map(regions, ~purrr::map2(.x, years, ~inner_fun(gbm.out, d, .x, .y)) %>% bind_rows) %>% bind_rows
}

d.preds <- predict_hold_out_year(gbm.out, d, unique(d$Region), yrs)
save(ri.out, cv.out, d.out, file="gbm_preds_eval_TEST.RData")

err.out2 <- unnest(err.out)
clrs <- c("black", "royalblue")
dir.create(plotDir <- file.path("../plots/gbm/models"), recursive=T, showWarnings=F)
png(file.path(plotDir, paste0("gbm_error.png")), width=3200, height=1600, res=200)
ggplot(err.out2 %>% filter(Region=="Arctic Tundra" & Run %in% 1:2), aes(`Number of Trees`, Error, group=interaction(Region, Year, Run, `Type of Error`), colour=`Type of Error`)) + geom_line(size=1) +
  scale_colour_manual("", values=clrs) + ggtitle("predictive error by GBM trees") +
  geom_point(data=err.out %>% filter(Region=="Arctic Tundra" & Run %in% 1:2), aes(x=Optimal_Trees, y=Optimal_Error, group=NULL, colour=NULL), size=1, colour="black") +
  geom_text(data=err.out %>% filter(Region=="Arctic Tundra" & Run %in% 1:2), aes(group=NULL, colour=NULL), colour="black") +
  theme_gray(base_size=16) + theme(legend.position="bottom", legend.box="horizontal", strip.background=element_blank()) +
  facet_wrap(~Region, ncol=3, scales="free")
dev.off()
