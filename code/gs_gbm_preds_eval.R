setwd("/workspace/UA/mfleonawicz/projects/GrowingSeason/workspaces")
pkgs <- list("rgdal", "raster", "maptools", "ggplot2", "data.table", "dplyr", "tidyr", "parallel", "gbm")
lapply(pkgs, function(x) library(x, character.only=T))

rasterOptions(tmpdir="/atlas_scratch/mfleonawicz/raster_tmp", chunksize=10e10, maxmemory=10e11)
load("data.RData") # d, d.stats, d.stats2, d.hm, sos, ecomask, yrs, cbpal
source("../code/gs_functions.R")
dem <- raster("/Data/Base_Data/GIS/GIS_Data/Raster/DEMs/PRISM_2km_DEM/AKCanada_2km_DEM_mosaic.tif") %>% resample(sos)
sos <- readAll(brick("../data/sos_1982_2010.tif"))

d <- dcast(d, Region + Year + Obs + x + y + SOS ~ Threshold, value.var="TDD")
d <- mutate(d, Elev=raster::extract(dem, cbind(x, y)))
setnames(d, c("Region", "Year", "Obs", "x", "y", "SOS", paste0("DOY_TDD", c("05", 10, 15, 20)), "Elev"))

set.seed(564)
n.trees=round(0.25*c(15000, 25000, 40000, 30000, 15000, 5000, 25000, 40000, 5000))

# build gbm models
gbm_explore <- function(i, data, n.trees, frac, years=sort(unique(d$Year)), by.year=TRUE){
  n <- if(by.year) length(years) else 1
  out <- vector("list", n)
  data <- filter(data, Year %in% years)
  gc()
  for(j in seq_along(out)){
  if(by.year){ d <- filter(data, Year==years[j]) } else { d <- data; rm(data); gc() }
  d <- group_by(d, Region, Year) %>% sample_frac(frac) %>% group_by(Region)
  d.train <- sample_frac(d, 0.9)
  d.test <- setdiff(d, d.train)
  d.gbm <- d.train %>% split(.$Region) %>% purrr::map2(n.trees, ~gbm(SOS ~ DOY_TDD05 + DOY_TDD10 + DOY_TDD15 + DOY_TDD20 + x + y + Elev, data=.x,
    distribution="gaussian", bag.fraction=0.5, cv.folds=5, train.fraction=1,
    interaction.depth=1, n.minobsinnode=5, n.trees=.y, shrinkage=0.1, verbose=FALSE, keep.data=FALSE, n.cores=1))
  d.gbm <- d.train %>% group_by %>% select(Region) %>% distinct(Region) %>% mutate(GBM1=d.gbm) %>% group_by(Region)
  rm(d.train); gc()
  d.bi <- d.gbm %>% do(BI=get_bi(., model=GBM1, plotDir, suffix=Region, saveplot=F))
  d.gbm <- data.table(suppressMessages(left_join(d.gbm, d.bi))) %>% group_by(Region)
  d.preds <- d.gbm %>% do(Predicted=get_preds(., model=GBM1, newdata=d.test, n.trees=BI, type.err="cv")) %>% group_by(Region)
  d.gbm <- d.gbm %>% mutate(CV=purrr::map(BI, ~.x$CV))
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
    d.gbm <- mutate(d.gbm, Year=years[j]) %>% select(Region, Year, GBM1, BI, CV)
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
system.time( dlist <- mclapply(1:32, gbm_explore, d, n.trees=n.trees, frac=0.5, by.year=T, mc.cores=32) )

# Extract tables of models, CV optimal trees, predictions, corrections, etc.
gbm.out <- lapply(dlist, "[[", 1)
cv.out <-  purrr::map(gbm.out, ~select(.x, Region, Year, CV)) %>% bind_rows
d.out <- rbindlist(lapply(dlist, "[[", 2))
lm.out <- lapply(dlist, "[[", 3)
d.out <- group_by(d.out, Region, Year, Source) %>% summarise(SOS=mean(SOS)) %>% bind_rows(filter(d.out))

# Spatial predictions
gbm_prediction_maps <- function(d, newdata, r, lm.pars=NULL, output="maps", n.cores=32){
  yrs <- sort(unique(newdata$Year))
  grp <- if("Year" %in% names(d[[1]])) list("Region", "Year") else list("Region")
  if("Year" %in% names(d[[1]])) newdata <- arrange(newdata, Year, Region, Obs)
  d <- mclapply(seq_along(d), function(i, x) x[[i]] %>% group_by_(.dots=grp) %>% do(Predicted=get_preds(., model=GBM1, newdata=newdata, n.trees=BI, type.err="cv")) %>% mutate(Run=i), x=d, mc.cores=n.cores)
  d <- d %>% purrr::map(~unnest(.x) %>% mutate(Cell=cellFromXY(r, select(newdata, x, y)), Year=as.integer(newdata$Year)) %>% nest(Cell, Predicted) %>% group_by(Region, Year, Run))
  if(!is.null(lm.pars)){
    d <- d %>% purrr::map2(lm.pars, ~suppressMessages(left_join(.x, .y, copy=T)))
    d <- d %>% purrr::map(~unnest(.x) %>% mutate(`Bias corrected`=(Predicted-intercept)/slope) %>% nest(Cell, Predicted, `Bias corrected`) %>% group_by(Region, Year, Run))
  }
  if(output=="table") return(bind_rows(d))
  gc()
  setPred <- function(i, r, d, values, yrs){
    r.pred <- raster(r)
    setValues(r.pred, NA)
    d2 <- filter(d, Year==yrs[i]) %>% select_(.dots=list(paste0("`", values, "`"), "Cell")) %>% unnest
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

# Run
pred.maps <- gbm_prediction_maps(gbm.out, d, subset(sos, 1), lm.out, n.cores=4)

save(cv.out, d.out, pred.maps, file="gbm_preds_eval_big.RData") # save all raw predictions maps in case needed

# Prep specific maps for plotting script
pred.maps.gbm.samples <- pred.maps[[1]] %>% purrr::map(~calc(.x, mean)) %>% stack
pred.maps.gbm1 <- pred.maps %>% purrr::map(~.x[[32]])

s1 <- stack(calc(pred.maps.gbm1$Predicted - sos, mean), calc(pred.maps.gbm1$`Bias corrected` - sos, mean))
sMean <- list(pred.maps[[1]] %>% purrr::map(~.x-sos), pred.maps[[2]] %>% purrr::map(~.x-sos))
sMean <- sMean %>% purrr::map(~Reduce("+", .x)/length(.x))
sMean2002 <- sMean %>% purrr::map(~subset(.x, match(2002, yrs))) %>% stack
sMean <- sMean %>% purrr::map(~calc(.x, mean)) %>% stack
names(s1) <- names(sMean) <- names(sMean2002) <- c("Prediction_deltas", "Bias_corrected_deltas")
s1 <- subset(s1, 1)
sMean <- subset(sMean, 1)
sMean2002 <- subset(sMean2002, 1)

save(cv.out, d.out, s1, sMean, sMean2002, file="gbm_preds_eval.RData")
