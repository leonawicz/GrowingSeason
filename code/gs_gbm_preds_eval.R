setwd("/workspace/UA/mfleonawicz/projects/GrowingSeason/workspaces")
pkgs <- list("rgdal", "raster", "maptools", "ggplot2", "data.table", "dplyr", "tidyr", "parallel", "gbm")
lapply(pkgs, function(x) library(x, character.only=T))

rasterOptions(tmpdir="/atlas_scratch/mfleonawicz/raster_tmp", chunksize=10e10, maxmemory=10e11)
load("data.RData") # d, d.stats, d.stats2, d.hm, sos, ecomask, yrs, cbpal
source("../code/gs_functions.R")
sos <- readAll(brick("../data/sos_1982_2010.tif"))
r <- calc(sos, mean)
#dem <- raster("/Data/Base_Data/GIS/GIS_Data/Raster/DEMs/PRISM_2km_DEM/AKCanada_2km_DEM_mosaic.tif") %>% resample(sos)

d <- dcast(d, Region + Year + Obs + x + y + SOS ~ Threshold, value.var="TDD")
#d <- mutate(d, Elev=raster::extract(dem, cbind(x, y)))
#d <- mutate(d, Cell=cellFromXY(sos, cbind(x,y)), Elev=raster::extract(dem, cbind(x, y)))
d <- mutate(d, Cell=cellFromXY(sos, cbind(x,y))) %>% select(-x, -y)
setnames(d, c("Region", "Year", "Obs", "SOS", paste0("DOY_TDD", c("05", 10, 15, 20)), "Cell"))
#setnames(d, c("Region", "Year", "Obs", "x", "y", "SOS", paste0("DOY_TDD", c("05", 10, 15, 20)), "Elev"))
#setnames(d, c("Region", "Year", "Obs", "x", "y", "SOS", paste0("DOY_TDD", c("05", 10, 15, 20)), "Cell", "Elev"))

cells <- (d %>% group_by(Cell) %>% summarise(Region=unique(Region), n=n()) %>% filter(n==29))$Cell
d <- filter(d, Cell %in% cells)

set.seed(564)
n.trees=(0.25)*c(3000, 2500, 4500, 6500, 3000, 1500, 3000, 6000, 1500)

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

gbm_explore <- function(i, data, n.trees, frac, years=sort(unique(data$Year)), by.year=TRUE){
  n <- if(by.year) length(years) else 1
  out <- vector("list", n)
  data <- filter(data, Year %in% years)
  gc()
  for(j in seq_along(out)){
    if(by.year){ d <- filter(data, Year==years[j]) } else { d <- data; rm(data); gc() }
    d <- group_by(d, Region, Year) %>% sample_frac(frac) %>% group_by(Region)
    d.train <- sample_frac(d, 0.9)
    d.test <- setdiff(d, d.train)
    d.gbm <- d.train %>% split(.$Region) %>% purrr::map2(n.trees, ~gbm(SOS ~ DOY_TDD05 + DOY_TDD10 + DOY_TDD15 + DOY_TDD20, data=.x,
      distribution="gaussian", bag.fraction=0.5, cv.folds=5, train.fraction=1,
      interaction.depth=1, n.minobsinnode=5, n.trees=.y, shrinkage=0.2, verbose=FALSE, keep.data=FALSE, n.cores=1))
    d.gbm <- d.train %>% group_by %>% select(Region) %>% distinct(Region) %>% mutate(GBM1=d.gbm) %>% group_by(Region)
    rm(d.train); gc()
    d.bi <- d.gbm %>% do(BI=get_bi(., model=GBM1, plotDir, saveplot=F))
    d.gbm <- data.table(suppressMessages(left_join(d.gbm, d.bi))) %>% group_by(Region)
    d.gbm <- d.gbm %>% mutate(CV=purrr::map(BI, ~.x$CV)) %>% group_by(Region)
    d.ri <- d.gbm %>% do(RI=get_ri(., model=GBM1, n.trees=BI, plotDir, saveplot=F))
    d.gbm <- data.table(suppressMessages(left_join(d.gbm, d.ri))) %>% group_by(Region)
    d.preds <- d.gbm %>% do(Predicted=get_preds(., model=GBM1, newdata=d.test, n.trees=BI, type.err="cv")) %>% group_by(Region)
    d.err <- d.gbm %>% group_by(Region) %>% do(Error=get_cv_err(.))
    d.gbm <- suppressMessages(left_join(d.gbm, d.err)) %>% group_by(Region)
    rm(d.bi, d.ri, d.err); gc()
    d.test$Predicted <- unnest(d.preds)$Predicted
    d.bias <- d.test %>% nest(-Region) %>% group_by(Region) %>% do(Region=.$Region, Predicted=.$Predicted[[1]], Coef=purrr::map2(.$SOS, .$Predicted, ~lm(.y ~ .x)$coefficients))
    d.bias <- d.bias %>% do(Region=.$Region, Predicted=.$Predicted, Coef=.$Coef, `Bias corrected`=(.$Predicted-.$Coef[[1]]["(Intercept)"])/.$Coef[[1]][".x"])
    d.coef <- d.bias %>% select(Region, Coef) %>% unnest(Region) %>% mutate(Coef=purrr::map(.$Coef, ~data.frame(t(.x[[1]])) %>% setnames(c("intercept", "slope")))) %>% unnest
    d.bias <- d.bias %>% select(-Coef) %>% unnest(Region) %>% unnest
    d.test <- suppressMessages(left_join(d.test, d.bias, copy=T))
    rm(d.bias); gc()
    d.test$Run <- i
    d.test <- d.test %>% select(Region, Year, Run, SOS, Predicted, `Bias corrected`) %>% setnames(c("Region", "Year", "Run", "Observed", "Predicted", "Bias corrected")) %>%
      melt(id.vars=c("Region", "Year", "Run"), value.name="SOS") %>% data.table %>% setnames(c("Region", "Year", "Run", "Source", "SOS")) %>%
      group_by(Region, Year, Run, Source) %>% summarise(SOS=round(mean(SOS)))
    if(by.year){
      d.coef <- mutate(d.coef, Year=years[j]) %>% select(Region, Year, intercept, slope)
      d.gbm <- mutate(d.gbm, Year=years[j]) %>% select(Region, Year, GBM1, RI, BI, CV, Error)
    }
    out[[j]] <- list(GBM=d.gbm, data=d.test, LM=d.coef)
    rm(d.gbm, d.test, d.coef); gc()
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
system.time( dlist <- mclapply(1:32, gbm_explore, d, n.trees=n.trees, frac=0.25, by.year=T, mc.cores=32) )

# Extract tables of models, CV optimal trees, predictions, corrections, etc.
gbm.out <- lapply(dlist, "[[", 1)
ri.out <-  purrr::map2(gbm.out, seq_along(gbm.out), ~select(.x, Region, Year, RI) %>% mutate(Run=.y) %>% unnest) %>% bind_rows %>% filter(Method=="CV")
cv.out <-  purrr::map(gbm.out, ~select(.x, Region, Year, CV)) %>% bind_rows
err.out <-  purrr::map2(gbm.out, seq_along(gbm.out), ~select(.x, Region, Year, Error) %>% mutate(Run=.y)) %>% bind_rows %>% unnest
d.out <- rbindlist(lapply(dlist, "[[", 2))
lm.out <- lapply(dlist, "[[", 3)
d.out <- group_by(d.out, Region, Year, Source) %>% summarise(SOS=mean(SOS)) %>% bind_rows(filter(d.out))

rm(dlist)
gc()
save(gbm.out, file="/atlas_scratch/mfleonawicz/projects/GrowingSeason/workspaces/gbm_models_all_huge.RData") 

# Spatial predictions
gbm_prediction_maps <- function(d, newdata, r, lm.pars=NULL, output="maps", n.cores=32){
  yrs <- sort(unique(newdata$Year))
  grp <- if("Year" %in% names(d[[1]])) list("Region", "Year") else list("Region")
  if("Year" %in% names(d[[1]])) newdata <- arrange(newdata, Year, Region, Obs)
  d <- mclapply(seq_along(d), function(i, x) x[[i]] %>% group_by_(.dots=grp) %>% do(Predicted=get_preds(., model=GBM1, newdata=newdata, n.trees=BI, type.err="cv")) %>% mutate(Run=i), x=d, mc.cores=n.cores)
  d <- d %>% purrr::map(~unnest(.x) %>% mutate(Cell=newdata$Cell, Year=as.integer(newdata$Year)) %>% nest(Cell, Predicted) %>% group_by(Region, Year, Run))
  if(!is.null(lm.pars)){
    d <- d %>% purrr::map2(lm.pars, ~suppressMessages(left_join(.x, .y, copy=T)))
    d <- d %>% purrr::map(~unnest(.x) %>% mutate(`Bias corrected`=(Predicted-intercept)/slope) %>% nest(Cell, Predicted, `Bias corrected`) %>% group_by(Region, Year, Run))
  }
  if(output=="table") return(bind_rows(d))
  gc()
  setPred <- function(i, r, d, values, yrs){
    r.pred <- raster(r)
    setValues(r.pred, NA)
    d2 <- select_(d, .dots=list(paste0("`", values, "`"), "Cell")) %>% filter(Year==yrs[i]) %>% unnest
    r.pred[d2$Cell] <- d2[[values]]
    r.pred
  }
  values <- if(is.null(lm.pars)) "Predicted" else c("Predicted", "Bias corrected")
  b <- vector("list", length(values))
  for(k in seq_along(values)){
    b[[k]] <- mclapply(d, function(x) brick(lapply(seq_along(yrs), setPred, r=r, d=x, values=values[k], yrs=yrs)), mc.cores=n.cores)
    gc()
  }
  names(b) <- values
  b
}

gbm_prediction_maps <- function(newdata, d, r, year.method="match", lm.pars=NULL, output="maps", simplify=TRUE, n.cores=32){
  yrs <- sort(unique(newdata$Year))
  gbm.yrs <- sort(unique(d[[1]]$Year))
  grp <- if("Year" %in% names(d[[1]])) list("Region", "Year") else list("Region")
  if("Year" %in% names(d[[1]])) newdata <- arrange(newdata, Year, Region, Obs)
  if(year.method=="match"){
    d <- mclapply(seq_along(d), function(i, x) x[[i]] %>% group_by_(.dots=grp) %>% do(Predicted=as.integer(round(get_preds(., model=GBM1, newdata=newdata, n.trees=BI, type.err="cv")))) %>% mutate(Run=i), x=d, mc.cores=n.cores)
  } else if(year.method=="random"){
    par_preds <- function(i, x) x[[i]] %>% group_by_(.dots=grp) %>% do(Predicted=as.integer(round(get_preds(., model=GBM1, newdata=newdata[[i]], n.trees=BI, type.err="cv")))) %>% mutate(Run=i)
    newdata <- purrr::map(seq_along(d), ~mutate(newdata, YearBackup=Year) %>% group_by(YearBackup) %>% mutate(Year=sample(gbm.yrs, 1, replace=TRUE)) %>% group_by)
    d <- mclapply(seq_along(d), par_preds, x=d, mc.cores=n.cores)
    newdata <- mutate(newdata[[1]], Year=YearBackup) %>% select(-YearBackup)
  }
  d <- d %>% purrr::map(~unnest(.x) %>% mutate(Cell=newdata$Cell, Year=as.integer(newdata$Year)) %>% nest(Cell, Predicted) %>% group_by(Region, Year, Run))
  if(!is.null(lm.pars)){
    d <- d %>% purrr::map2(lm.pars, ~suppressMessages(left_join(.x, .y, copy=T)))
    d <- d %>% purrr::map(~unnest(.x) %>% mutate(`Bias corrected`=(Predicted-intercept)/slope) %>% nest(Cell, Predicted, `Bias corrected`) %>% group_by(Region, Year, Run))
  }
  if(output=="table") return(bind_rows(d))
  gc()
  setPred <- function(d, r, values){
    d2 <- select_(d, .dots=list(paste0("`", values, "`"), "Cell")) %>% unnest
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
#pred.maps <- gbm_prediction_maps(d, gbm.out, r, n.cores=16)

load("tdd_table.RData")
rcps <- unique(d.tdd$RCP)[-1]
models <- unique(d.tdd$Model)[-1]
pred.maps <- purrr::map(vector("list", length(rcps)), ~vector("list", length(models)) %>% setNames(models)) %>% setNames(rcps)
set.seed(1)
#system.time({
#pred.maps <- d.tdd.tmp %>% split(.$RCP) %>%
#  purrr::map(~.x %>% split(.$Model) %>%
#    purrr::map(~select_(.x, .dots=list("-RCP", "-Model")) %>% gbm_prediction_maps(gbm.out[1:16], r, year.method="random", n.cores=16)))
#})
system.time({
for(i in seq_along(rcps)){
  for(j in seq_along(models)){
    d.gcm <- d.tdd %>% data.table %>% filter(RCP==rcps[i] & Model==models[j]) %>% select(-RCP, -Model, -x, -y, -Elev, -SOS)
    pred.maps[[i]][[j]] <- gbm_prediction_maps(d.gcm, gbm.out[1:10], r, year.method="random", n.cores=10)
  }
}
})
saveRDS(pred.maps, file="/atlas_scratch/mfleonawicz/projects/GrowingSeason/workspaces/gbm_preds_gcm_all.RData") 

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

d.hoy <- d.preds %>% select(-SOS) %>% mutate(Source="Predicted HOY") %>% rename(SOS=Predicted)
d.hoy <- d.preds %>% select(-Predicted) %>% distinct %>% mutate(Source="Global observed") %>% bind_rows(d.hoy)
d.ts <- bind_rows(filter(d.out, Source!="Bias corrected"), d.hoy) %>% mutate(Source=factor(Source, levels=c("Observed", "Predicted", "Global observed", "Predicted HOY")))

d.proj <- extract_to_dt(pred.maps, eco_shp, mean, "RCP 6.0", "GFDL-CM3")

d.proj.mean <- d.proj %>% group_by(Region, Year, Source) %>% summarise(SOS=mean(SOS), Run=1)
d.smooth <- filter(d.ts, is.na(Run) & Source=="Predicted") %>% bind_rows(filter(d.proj, Year > 2010)) %>% group_by(Region, Year) %>% summarise(SOS=mean(SOS), Source="Trend", Run=1)

dir.create(plotDir <- file.path("../plots/gbm/models"), recursive=T, showWarnings=F)

# historical predictions over observations time series
clrs <- c("peru", "royalblue", "black", "red")
png(file.path(plotDir, paste0("gbm_TSpreds_byRegion_plusGCMs.png")), width=3200, height=1600, res=200)
ggplot(filter(d.out, !is.na(Run) & Source!="Bias corrected"), aes(x=Year, y=SOS, colour=Source, group=interaction(Source, Run))) +
  scale_color_manual(values=clrs) + 
  geom_line(size=1, alpha=0.25) +
  #geom_line(data=filter(d.out, Source=="Bias corrected" & is.na(Run)), colour="#B8860B", size=1, linetype=1) +
  geom_line(data=filter(d.out, Source=="Observed" & is.na(Run)), colour="black", size=2) + 
  geom_line(data=filter(d.out, Source=="Predicted" & is.na(Run)), colour="black", size=2) +
  geom_line(data=filter(d.out, Source=="Predicted" & is.na(Run)), colour=clrs[2], size=1) +
  geom_line(data=filter(d.out, Source=="Observed" & is.na(Run)), colour=clrs[1], size=1) + 
  geom_line(data=filter(d.out, Source=="Predicted" & is.na(Run)), colour=clrs[2], size=1, linetype=2) +
  geom_line(data=d.proj, colour="#00000030", size=1) + 
  geom_smooth(data=d.smooth) +
  theme_bw() + theme(legend.position="bottom") + ggtitle("Observed and modeled start of growing season") +
  #scale_x_continuous(breaks=c(1982,1990,2000,2010)) +
  guides(fill=guide_legend(override.aes=list(alpha=1)), colour=guide_legend(override.aes=list(alpha=1))) +
  facet_wrap(~Region, ncol=3, scales="free")
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
save(ri.out, cv.out, d.out, s1, sMean, sMean2002, d.preds, file="gbm_preds_eval.RData")

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
