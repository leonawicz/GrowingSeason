#############################################
#### Growing season thaw degree days GBM ####
#############################################

#### Script author:  Matthew Leonawicz ####
#### Maintainted by: Matthew Leonawicz ####
#### Last updated:   02/29/2015        ####

# @knitr setup
setwd("/workspace/UA/mfleonawicz/projects/GrowingSeason/workspaces")

library(rasterVis)
library(maptools)
library(ggplot2)
library(data.table)
library(dplyr)
library(tidyr)
library(parallel)
library(gbm)

rasterOptions(tmpdir="/atlas_scratch/mfleonawicz/raster_tmp", chunksize=10e10, maxmemory=10e11)
load("data.RData") # d, d.stats, d.stats2, d.hm, sos, ecomask, yrs, cbpal
source("../code/gs_functions.R")
sw.shp <- shapefile("/big_scratch/mfleonawicz/Alf_Files_20121129/statewide_shape/Alaska_Albers_ESRI.shp")
dem <- raster("/Data/Base_Data/GIS/GIS_Data/Raster/DEMs/PRISM_2km_DEM/AKCanada_2km_DEM_mosaic.tif") %>% resample(sos)
dir.create(plotDir <- file.path("../plots/gbm/models"), recursive=T, showWarnings=F)
sos <- readAll(brick("../data/sos_1982_2010.tif"))

d <- dcast(d, Region + Year + Obs + x + y + SOS ~ Threshold, value.var="TDD")
d <- mutate(d, Elev=raster::extract(dem, cbind(x, y)))
setnames(d, c("Region", "Year", "Obs", "x", "y", "SOS", paste0("DOY_TDD", c("05", 10, 15, 20)), "Elev"))

set.seed(564)
d.sub <- group_by(d, Region, Year) %>% sample_frac(0.50) %>% setorder(Region, Year) %>% group_by(Region)

##############################################
set.seed(564)
n.trees=round(0.25*c(15000, 25000, 40000, 30000, 15000, 5000, 25000, 40000, 5000))

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

system.time( dlist <- mclapply(1:32, gbm_explore, d, n.trees=n.trees, frac=0.5, by.year=T, mc.cores=32) )
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

pred.maps <- gbm_prediction_maps(gbm.out, d, subset(sos, 1), lm.out, n.cores=4)
save(cv.out, d.out, pred.maps, file="gbm_preds_eval_big.RData")

rTheme <- function(region=colorRampPalette(c("darkred", "firebrick1", "white", "royalblue", "darkblue"))(19), ...){
    theme <- custom.theme.2(region=region, ...)
    theme$strip.background$col <- theme$strip.shingle$col <- theme$strip.border$col <- theme$panel.background$col <- theme$axis.line$col <- "transparent"
    theme
}

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

rng <- range(c(s1[], sMean[]), na.rm=T)
rng <- range(rng, -rng)
at.vals <- seq(rng[1], rng[2], length=19)
colkey <- list(at=at.vals, labels=list(labels=round(at.vals), at=at.vals))

png(file.path(plotDir, paste0("gbm_pred_deltas_gbm1.png")), width=3200, height=1600, res=200)
p <- levelplot(s1, maxpixels=ncell(s1), main=paste("1982-2010 mean start of growing season deltas: single GBM"), par.settings=rTheme, xlab=NULL, ylab=NULL, scales=list(draw=F), contour=F, margin=F, at=at.vals, colorkey=colkey) +
        latticeExtra::layer(sp.polygons(sw.shp, fill='transparent', alpha=0.3))
print(p)
dev.off()

png(file.path(plotDir, paste0("gbm_pred_deltas_gbmMean.png")), width=3200, height=1600, res=200)
p <- levelplot(sMean, maxpixels=ncell(sMean), main=paste("1982-2010 mean start of growing season deltas: mean GBM"), par.settings=rTheme, xlab=NULL, ylab=NULL, scales=list(draw=F), contour=F, margin=F, at=at.vals, colorkey=colkey) +
        latticeExtra::layer(sp.polygons(sw.shp, fill='transparent', alpha=0.3))
print(p)
dev.off()

rng <- range(sMean2002[], na.rm=T)
rng <- range(rng, -rng)
at.vals <- seq(rng[1], rng[2], length=19)
colkey <- list(at=at.vals, labels=list(labels=round(at.vals), at=at.vals))
png(file.path(plotDir, paste0("gbm_pred_deltas_gbmMean2002.png")), width=3200, height=1600, res=200)
p <- levelplot(sMean2002, maxpixels=ncell(sMean2002), main=paste("2002 start of growing season deltas: mean GBM"), par.settings=rTheme, xlab=NULL, ylab=NULL, scales=list(draw=F), contour=F, margin=F, at=at.vals, colorkey=colkey) +
        latticeExtra::layer(sp.polygons(sw.shp, fill='transparent', alpha=0.3))
print(p)
dev.off()

arctic <- Which(ecomask==3)
arctic[arctic==0] <- NA
arctic <- trim(arctic)

rTheme <- function(region=brewer.pal(9, "Blues"), ...){
    theme <- custom.theme.2(region=region, ...)
    theme$strip.background$col <- theme$strip.shingle$col <- theme$strip.border$col <- theme$panel.background$col <- theme$axis.line$col <- "transparent"
    theme
}

png(file.path(plotDir, paste0("gbm_pred_byModel.png")), width=3200, height=1600, res=200)
p <- levelplot(pred.maps.gbm.samples, maxpixels=ncell(pred.maps.gbm.samples), main=paste("1982-2010 mean start of growing season GBM comparisons"), par.settings=rTheme, xlab=NULL, ylab=NULL, scales=list(draw=F), contour=F, margin=F) +#, at=at.vals, colorkey=colkey)
        latticeExtra::layer(sp.polygons(sw.shp, fill='transparent', alpha=0.3))
print(p)
dev.off()

sos2002 <- subset(sos, 21)
s <- stack(sos2002, sos2002 + sMean2002)
names(s) <- c("Observed", "Predicted")
png(file.path(plotDir, paste0("gbm_pred_gbmMean2002.png")), width=3200, height=1600, res=200)
p <- levelplot(s, maxpixels=ncell(s), main=paste("2002 start of growing season: mean GBM"), par.settings=rTheme, xlab=NULL, ylab=NULL, scales=list(draw=F), contour=F, margin=F) +
        latticeExtra::layer(sp.polygons(sw.shp, fill='transparent', alpha=0.3))
print(p)
dev.off()

s <- crop(s, arctic) %>% mask(arctic)
png(file.path(plotDir, paste0("gbm_pred_gbmMean2002_ArcticTundra.png")), width=3200, height=1600, res=200)
p <- levelplot(s, maxpixels=ncell(s), main=paste("2002 start of Arctic Tundra growing season: mean GBM"), par.settings=rTheme, xlab=NULL, ylab=NULL, scales=list(draw=F), contour=F, margin=F) +
        latticeExtra::layer(sp.polygons(sw.shp, fill='transparent', alpha=0.3))
print(p)
dev.off()
##############################################

system.time({
d.gbm <- d.sub %>% do(GBM1=gbm(SOS ~ DOY_TDD05 + DOY_TDD10 + DOY_TDD15 + DOY_TDD20 + x + y + Elev, data=.,
    distribution="gaussian", bag.fraction=0.5, cv.folds=5, train.fraction=1,
    #weights=sqrt(SOS),
    interaction.depth=1, n.minobsinnode=20, n.trees=5000, shrinkage=0.01,
    keep.data=F,
    verbose=FALSE, n.cores=1)) %>% group_by(Region)
})

# Error curves
d.bi <- d.gbm %>% do(BI=get_bi(., model=GBM1, plotDir, suffix=Region, saveplot=T))
d.gbm <- data.table(left_join(d.gbm, d.bi)) %>% group_by(Region)

# Relative influence
d.ri <- d.gbm %>% do(RI=get_ri(., model=GBM1, n.trees=BI, plotDir, suffix=Region, saveplot=T))
d.gbm <- data.table(left_join(d.gbm, d.ri)) %>% group_by(Region)

png(file.path(plotDir, paste0("gbm_RI_byRegion.png")), width=1600, height=1600, res=200)
    tmp <- select(d.gbm, RI) %>% unnest %>% data.table %>% filter(Method=="CV") %>% group_by(Region) %>%
        mutate(barorder=as.numeric(strsplit(paste(RI, collapse=","), ",")[[1]][1])) %>% mutate(Region=factor(Region, levels=unique(Region)[order(unique(barorder))]))
    g <- ggplot(tmp, aes(Region, RI, fill=Predictor)) + geom_bar(stat="identity", position="stack") +
    scale_fill_manual(values=cbpal[-2]) + theme_bw() + theme(legend.position="bottom") + coord_flip() +
    ggtitle(paste("Predictor relative influence on start of growing season"))
    print(g)
dev.off()

# Partial dependence
#d.gbm %>% do(get_pd2(., model=GBM1, plotDir, suffix=Region, saveplot=T))
#gbm1 <- (d.gbm %>% filter(Region=="Arctic Tundra"))$GBM1[[1]]
#get_pd(source_data=d.tdd %>% filter(Region=="Arctic Tundra"), x=DOY_TDD, y=SOS, outDir=plotDir, model=gbm1, vars=1:4, order.by.ri=TRUE)
d.tdd <- select(d.sub, -Obs, -Year) %>% melt(id.vars=c("Region", "SOS"), variable.name="Var", value.name="Val") %>% group_by(Region, Var)
d.pd <- d.gbm %>% do(PD=get_pd(., source_data=d.tdd, x=Val, y=SOS, outDir=plotDir, model=GBM1, vars=1:4, order.by.ri=TRUE, suffix=Region, saveplot=T)) %>% group_by(Region)
d.gbm <- data.table(left_join(d.gbm, d.pd)) %>% group_by(Region)

d.pd.test <- d.pd %>% group_by %>% select(-Region) %>% unnest %>% mutate(Var=sapply(strsplit(Var, ":"), "[", 1), Var=factor(Var, levels=sort(unique(Var)))) %>% group_by(Region, Var)
png(file.path(plotDir, paste0("gbm_PD_test1.png")), width=3200, height=1600, res=200)
ggplot(d.pd.test %>% filter(Region %in% c("Arctic Tundra", "Bering Tundra")) %>% mutate(Ymin=min(Prob)), aes(x=Val)) +
    facet_wrap(Region ~ Var, ncol=4, scales="fixed") + labs(x="DOY_TDD", y="SOS") +
    geom_ribbon(aes(ymin=Ymin, ymax=Prob), fill="orange") + geom_line(aes(x=x, y=y), size=1)
dev.off()

png(file.path(plotDir, paste0("gbm_PD_test2.png")), width=3200, height=1600, res=200)
    ggplot(d.pd.test %>% filter(Region %in% c("Arctic Tundra", "Bering Tundra")) %>% mutate(Ymin=min(Prob)), aes(x=Val)) +
    facet_wrap(~ Var, ncol=2) + labs(x="DOY_TDD", y="SOS") +
    geom_ribbon(aes_string(ymin="Ymin", ymax="Prob", colour="Region", fill="Region"), alpha=0.5) + geom_line(aes_string(x="x", y="y", colour="Region"), size=1)
dev.off()

# Time series
d.preds <- d.gbm %>% do(Pred=get_preds(., model=GBM1, newdata=d.sub, n.trees=BI, type.err="cv")) %>% group_by(Region)

# Prep to test exchangeability of pairs of regional GBMs
region.pair <- c("Arctic Tundra", "Coastal Rainforests")
d.sub.tmp <- copy(d.sub)
d.sub.tmp2 <- d.sub.tmp <- filter(d.sub.tmp, Region %in% region.pair)
d.gbm.tmp2 <- d.gbm.tmp <- copy(d.gbm)
d.gbm.tmp <- filter(d.gbm.tmp, Region %in% region.pair)
d.gbm.tmp2 <- filter(d.gbm.tmp, Region %in% region.pair)
d.gbm.tmp2$Region <- d.gbm.tmp$Region[2:1]
d.preds.tmp <- d.gbm.tmp %>% do(Pred=get_preds(., model=GBM1, newdata=d.sub.tmp, n.trees=BI, type.err="cv")) %>% group_by(Region)
d.preds.tmp2 <- d.gbm.tmp2 %>% do(Pred=get_preds(., model=GBM1, newdata=d.sub.tmp, n.trees=BI, type.err="cv")) %>% group_by(Region)

d.gbm1 <- copy(d.gbm) # for shiny app
d.gbm.data <- copy(d.sub) # for shiny app

d.gbm <- data.table(left_join(d.gbm, d.preds)) %>% group_by(Region)

# Time series regional average predictions
d.sub$Pred <- unnest(d.preds)$Pred
d.sub %>% select(Region, Year, SOS, Pred) %>% setnames(c("Region", "Year", "Obs", "Pred")) %>%
    melt(id.vars=c("Region", "Year"), value.name="SOS") %>% data.table %>% setnames(c("Region", "Year", "Source", "SOS")) %>%
    group_by(Region, Year, Source) %>% summarise(SOS=round(mean(SOS))) -> d.sub2

d.gbm.preds <- copy(d.sub2) # for shiny app

png(file.path(plotDir, paste0("gbm_TSpreds_byRegion.png")), width=3200, height=1600, res=200)
ggplot(d.sub2, aes(x=Year, y=SOS, colour=Source)) + scale_color_manual(values=c("black", "red")) + geom_line(size=1) + geom_point() +
    theme_bw() + theme(legend.position="bottom") + ggtitle("Observed and modeled start of growing season") +
    scale_x_continuous(breaks=c(1982,1990,2000,2010)) +
    facet_wrap(~Region, ncol=3)
dev.off()
png(file.path(plotDir, paste0("gbm_TSpreds_byRegionSD.png")), width=3200, height=1600, res=200)
ggplot(d.sub2 %>% group_by(Region, Source) %>% summarise(SD=sd(SOS)), aes(x=Region, y=SD, fill=Source)) + scale_fill_manual(values=c("black", "red")) + geom_bar(stat="identity", position="dodge") +
    theme_bw() + theme(legend.position="bottom") + ggtitle("Observed and modeled start of growing season period SD")
dev.off()

# save workspace for shiny app
save(d.gbm1, d.gbm.preds, d.gbm.data, file="appdata_gbm.RData")

# Time series predictions GBM exchangeability
d.sub.tmp$Pred <- unnest(d.preds.tmp)$Pred
d.sub.tmp2$Pred <- unnest(d.preds.tmp2)$Pred
d.sub.tmp %>% select(Region, Year, SOS, Pred) %>% setnames(c("Region", "Year", "Obs", "Pred")) %>%
    melt(id.vars=c("Region", "Year"), value.name="SOS") %>% data.table %>% setnames(c("Region", "Year", "Source", "SOS")) %>%
    group_by(Region, Year, Source) %>% summarise(SOS=round(mean(SOS))) -> d.sub2.tmp
d.sub.tmp2 %>% select(Region, Year, SOS, Pred) %>% setnames(c("Region", "Year", "Obs", "Pred")) %>%
    melt(id.vars=c("Region", "Year"), value.name="SOS") %>% data.table %>% setnames(c("Region", "Year", "Source", "SOS")) %>%
    group_by(Region, Year, Source) %>% summarise(SOS=round(mean(SOS))) -> d.sub2.tmp2
d.sub2.tmp2$Region <- paste(d.sub2.tmp2$Region, "(swapped GBM)")
d.sub2.tmp <- rbind(d.sub2.tmp, d.sub2.tmp2)
d.sub2.tmp <- mutate(d.sub2.tmp, Region=factor(Region, levels=unique(Region)[c(1,2,3,4)]))

png(file.path(plotDir, paste0("gbm_TSpreds_TestRegionalModels.png")), width=3200, height=1600, res=200)
ggplot(d.sub2.tmp, aes(x=Year, y=SOS, colour=Source)) + scale_color_manual(values=c("black", "red")) + geom_line(size=1) + geom_point() +
    theme_bw() + theme(legend.position="bottom") + ggtitle("Observed and modeled start of growing season when swapping regional GBM predictive models") +
    scale_x_continuous(breaks=c(1982,1990,2000,2010)) +
    facet_wrap(~Region, ncol=2)
dev.off()

# Spatial predictions
r.sos.mean <- calc(sos, mean, na.rm=TRUE)
r.sos.mean[is.nan(r.sos.mean)] <- NA
d.sub[, Cell:=cellFromXY(r.sos.mean, d.sub[, c("x","y"), with=FALSE])]

setPred <- function(i, r, d, yrs){
    r.pred <- raster(r)
    setValues(r.pred, NA)
    d %>% filter(Year==yrs[i]) %>% select(Pred, Cell) -> d2
    r.pred[d2$Cell] <- d2$Pred
    r.pred
}

r.pred <- calc(stack(lapply(1:length(yrs), setPred, r=r.sos.mean, d=d.sub, yrs=yrs)), mean, na.rm=TRUE)
r.pred[is.nan(r.pred)] <- NA

s <- stack(trim(mask(r.sos.mean, r.pred)), trim(r.pred), trim(r.pred-mask(r.sos.mean, r.pred)))
names(s) <- c("Observed", "Modeled", "Difference")

# Theme settings
revRasterTheme <- function(pch = 19, cex = 0.7, region=terrain.colors(30), ...){
    theme <- custom.theme.2(pch = pch, cex = cex, region = region, ...)
    theme$strip.background$col <- theme$strip.shingle$col <- theme$strip.border$col <- "transparent"
    theme$add.line$lwd = 0.4
    theme
}

png(file.path(plotDir, paste0(region,"_gbm4_predsMap1.png")), width=1600, height=1600, res=200)
p <- levelplot(subset(s, 1:2), maxpixels=ncell(r.sos.mean), main=paste("1982-2010 mean start of", region, "growing season"), par.settings=revRasterTheme, contour=F, margin=F)# +#, at=at.vals, colorkey=colkey)
        #layer(sp.polygons(eco_shp, fill='transparent', alpha=0.3))
print(p)
dev.off()

png(file.path(plotDir, paste0(region,"_gbm4_predsMap2.png")), width=1600, height=1600, res=200)
p <- levelplot(subset(s, 3), maxpixels=ncell(r.sos.mean), main=paste("1982-2010 pred. - obs. mean start of", region, "growing season"), par.settings=revRasterTheme, contour=F, margin=F)# +#, at=at.vals, colorkey=colkey)
        #layer(sp.polygons(eco_shp, fill='transparent', alpha=0.3))
print(p)
dev.off()


