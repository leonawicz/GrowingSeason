setwd("/atlas_scratch/mfleonawicz/projects/GrowingSeason/workspaces")
#setwd("C:/github/GrowingSeason/workspaces")

library(maptools)
library(rasterVis)
library(data.table)
library(dplyr)

source("../code/gs_functions.R")

files <- list.files(pattern="^sos_gbm_preds_.*.RData")
models <- unlist(strsplit( sapply(strsplit(files, "sos_gbm_preds_"), "[", 2), "\\.RData"))
yrs.sos <- 1982:2010
yrs.narr <- 1979:2010
yrs.hist <- 1958:2005
yrs.proj <- 2006:2100
load(files[1])
b.sos <- brick("../data/sos/sos_1982_2010.tif")
names(b.sos) <- paste0("SOS_", yrs.sos)

# Alaska ecoregion level two mask
shpDir <- "/atlas_scratch/mfleonawicz/projects/DataExtraction/data/shapefiles"
#shpDir <- "C:/github/DataExtraction/data/shapefiles"
eco_shp <- shapefile(file.path(shpDir, "AK_ecoregions/akecoregions.shp"))
eco_shp <- spTransform(eco_shp, CRS(projection(b.narr)))
eco_shp <- unionSpatialPolygons(eco_shp, eco_shp@data$LEVEL_2)
ecomask <- rasterize(eco_shp, b.narr)
regions <- names(eco_shp)

myFun <- function(x, ...) c(Mean=mean(x,...), SD=sd(x,...), quantile(x, prob=c(0.05, 0.95),...))

get_brick_stats <- function(x, shp, years, scenario, model, map, ...){
    stopifnot(any(class(x) %in% c("RasterBrick", "RasterStack")))
    ids <- sapply(slot(shp, "polygons"), function(x) slot(x, "ID"))
    x <- extract(x, shp, ...)
    if(is.list(x)) stop("Unequal length results returned by polygon. Check summary function.")
    n <- nrow(x)/length(ids)
    x <- lapply(1:n, function(i, x) x[rownames(x)==unique(rownames(x))[i],], x=x)
    x <- lapply(x, function(x, years) data.table(Year=years, Stat=unique(rownames(x)), t(x)), years=years)
    x <- rbindlist(x)
    setnames(x, c("Year", "Stat", ids))
    x <- melt(x, id.vars=c("Year", "Stat"), variable.name="Region", value.name="SOS")
    x <- mutate(x, Scenario=scenario, Model=model, Map=map)
    setcolorder(x, c("Scenario", "Model", "Map", "Region", "Year", "Stat", "SOS"))
    x
}

get_stat_maps <- function(x, ...) calc(x, myFun, ...)

system.time({
maplist.sos <- vector("list", length(models))
for(i in 1:length(models)){
    model <- models[i]
    load(files[i])
    if(i==1) d <- get_brick_stats(x=b.sos, shp=eco_shp, years=yrs.sos, scenario="Historical", model="SOS", map="Baseline", fun=myFun, na.rm=TRUE)
    if(i==1) d <- bind_rows(d, get_brick_stats(x=b.narr, shp=eco_shp, years=yrs.narr, scenario="Historical", model="NARR", map="Baseline", fun=myFun, na.rm=TRUE))
    d <- bind_rows(d, get_brick_stats(x=b.hist, shp=eco_shp, years=yrs.hist, scenario="Historical", model=model, map="Original", fun=myFun, na.rm=TRUE))
    d <- bind_rows(d, get_brick_stats(x=b.hist.qm, shp=eco_shp, years=yrs.hist, scenario="Historical", model=model, map="Quantile-mapped", fun=myFun, na.rm=TRUE))
    d <- bind_rows(d, get_brick_stats(x=b.rcp60, shp=eco_shp, years=yrs.proj, scenario="RCP 6.0", model=model, map="Original", fun=myFun, na.rm=TRUE))
    d <- bind_rows(d, get_brick_stats(x=b.rcp60.qm, shp=eco_shp, years=yrs.proj, scenario="RCP 6.0", model=model, map="Quantile-mapped", fun=myFun, na.rm=TRUE))
    d <- bind_rows(d, get_brick_stats(x=b.rcp85, shp=eco_shp, years=yrs.proj, scenario="RCP 8.5", model=model, map="Original", fun=myFun, na.rm=TRUE))
    d <- bind_rows(d, get_brick_stats(x=b.rcp85.qm, shp=eco_shp, years=yrs.proj, scenario="RCP 8.5", model=model, map="Quantile-mapped", fun=myFun, na.rm=TRUE))
    
    if(i==1) x0 <- get_stat_maps(x=b.sos, na.rm=TRUE)
    x1narr <- get_stat_maps(x=b.narr, na.rm=TRUE)
    x1 <- get_stat_maps(x=b.hist, na.rm=TRUE)
    x1m <- get_stat_maps(x=b.hist.qm, na.rm=TRUE)
    x1narr.deltas <- x1narr-x0
    x1.deltas <- x1-x0#lapply(1:length(x0), function(i, x, y) y[[i]]-x[[i]], x=x0, y=x1)
    x1m.deltas <- x1m-x0#lapply(1:length(x0), function(i, x, y) y[[i]]-x[[i]], x=x0, y=x1m)
    #names(x1.deltas) <- names(x1m.deltas) <- names(doytdd0)
    maplist.sos[[i]] <- list(diffx1narrx0=x1narr.deltas, diffx1x0=x1.deltas, diffx1mx0=x1m.deltas)
}
})

names(maplist.sos) <- models
d.sos <- d
save(d.sos, maplist.sos, ecomask, eco_shp, file="sos_final.RData")

# test plot: single pixel
library(ggplot2)
clrs <- c("black", "darkorange", "purple")
ind <- 100000
x0 <- as.numeric(b.narr[ind])
x1 <- as.numeric(b.hist[ind])
x1m <- as.numeric(b.hist.qm[ind])
d <- data.table(Source=factor(rep(c("NARR", "GCM", "QMAP"), times=c(length(x0), length(x1), length(x1m))), levels=c("NARR", "GCM", "QMAP")), SOS=c(x0, x1, x1m))

ggplot(d, aes(x=SOS, fill=Source, colour=Source)) + geom_line(stat="density", size=1) +
    xlim(50,150) + scale_colour_manual(values=clrs) + theme(legend.position="bottom") +
    labs(title="Sample pixel SOS time series distribution\nNARR 1979-2010 and GCM 1958-2005") 

###################################################################

setwd("C:/github/GrowingSeason/workspaces")
library(data.table)
library(dplyr)
library(rasterVis)
library(ggplot2)

load("sos_final.RData")

sos.means <- d.sos %>% filter(Model=="SOS") %>% group_by(Region, Stat) %>% summarise(SOS_Mean=mean(SOS))
d2.sos <- left_join(d.sos, sos.means) %>% mutate(`SOS Deltas`=SOS-SOS_Mean)

# save workspace for shiny app
d.sos.qmap <- d2.sos
save(d.sos.qmap, maplist.sos, file="C:/github/shiny-apps/gs/appdata_qmap_sos.RData")

clrs <- c("black", "darkorange", "purple")
scens <- c("Historical", "RCP 6.0")
mods <- c("SOS", "NARR")
spstat <- "Mean"
regions <- names(eco_shp)

dsub <- d.sos %>% filter(Scenario %in% scens & Model %in% mods)
dsub2 <- d2.sos %>% filter(Scenario %in% scens & Model %in% mods)

# Time series and histogram DOY TDD and deltas
tsDir <- "../plots/sos_qmap/timeseries"
denDir <- "../plots/sos_qmap/distributions"
png(paste0(tsDir,"/ts_full_", spstat, "_sos_byRegion.png"), width=3200, height=2400, res=200)
ggplot(dsub %>% filter(Stat==spstat), aes(x=Year, y=SOS, group=interaction(Model, Map), colour=Model)) +
    geom_line(data=dsub %>% filter(Stat==spstat & Model=="SOS"), size=1) + geom_line() +
    scale_colour_manual(values=clrs) + theme(legend.position="bottom") + facet_wrap(~ Region, scales="fixed") +
    labs(y="SOS", title=paste("SOS baseline and", mods[2], scens[scens!="Historical"][1], "original and qm-based regional", spstat, "SOS predictions"))
dev.off()

png(paste0(denDir, "/density_historical_", spstat, "_doytdd5pct20pct_byRegion.png"), width=3200, height=2400, res=200)
ggplot(dsub %>% filter(Pct %in% c("0.05", "0.2") & Stat==spstat & Scenario=="Historical"), aes(x=`DOY TDD`, group=interaction(Pct, Map), linetype=Pct, colour=Map)) + geom_line(stat="density", size=1) +
    scale_colour_manual(values=clrs) + theme(legend.position="bottom", legend.box="horizontal") + facet_wrap(~ Region, scales="free_y") +
    labs(x="DOY TDD 5% and 20%", title=paste("NARR baseline and historical", mods[2], "original and quantile-mapped regional", spstat, "DOY TDD 5% and 20%"))
dev.off()

png(paste0(tsDir, "/ts_deltas_historical_", spstat, "_doytdd10pct_byRegion.png"), width=3200, height=2400, res=200)
ggplot(dsub2 %>% filter(Pct=="0.1" & Stat==spstat & Year %in% 1979:2010), aes(x=Year, y=`DOY TDD Deltas`, group=interaction(Model, Map), colour=Map)) +
    geom_line(data=dsub2 %>% filter(Pct=="0.1" & Stat==spstat & Model=="NARR"), size=1) + geom_line() +
    scale_colour_manual(values=clrs) + theme(legend.position="bottom") + facet_wrap(~ Region, scales="free_y") +
    labs(y=expression(DOY[GCM]-mean(DOY[NARR])), title=paste(mods[2], "original and quantile-mapped regional", spstat, "DOY TDD 10% deltas"))
dev.off()

png(paste0(denDir, "/density_deltas_historical_", spstat, "_doytdd5pct20pct_byRegion.png"), width=3200, height=2400, res=200)
ggplot(dsub2 %>% filter(Pct %in% c("0.05", "0.2") & Stat==spstat & Scenario=="Historical"), aes(x=`DOY TDD Deltas`, group=interaction(Pct, Map), linetype=Pct, colour=Map)) + geom_line(stat="density", size=1) +
    scale_colour_manual(values=clrs) + theme(legend.position="bottom", legend.box="horizontal") + facet_wrap(~ Region, scales="free_y") +
    labs(x=expression(DOY[GCM]-DOY[NARR[mean]]), title=paste(mods[2], "original and quantile-mapped regional", spstat, "DOY TDD 5% and 20% deltas"))
dev.off()

# Maps
# Theme settings
revRasterTheme <- function (pch = 19, cex = 0.7, region=colorRampPalette(c("navyblue", "gray90", "darkred"))(9), ...){   #brewer.pal(9, "RdBu")[-1]
    theme <- custom.theme.2(pch = pch, cex = cex, region = region, ...)
    theme$strip.background$col <- theme$strip.shingle$col <- theme$strip.border$col <- "transparent"
    theme$panel.background$col <- "gray20"
    theme$add.line$lwd = 0.4
    theme
}

for(i in 1:9){
for(j in 1:4){
region <- regions[i]
msk <- Which(ecomask==match(region, regions))
msk[msk==0] <- NA
stat <- c("mean", "SD", "pct05", "pct95")[j]
stat2 <- c("mean", "std. dev.", "5th percentile", "95th percentile")[j]
r1 <- subset(maplist[[mods[2]]]$diffx1x0$`0.1`, j)
r2 <- subset(maplist[[mods[2]]]$diffx1mx0$`0.1`, j)
r1[is.nan(r1)] <- NA
r2[is.nan(r2)] <- NA
s <- stack(r1, r2)
s <- mask(s, msk)
s <- trim(s)
names(s) <- c("Original", "Mapped")
rng <- range(as.numeric(cellStats(s, range)))
at.vals <- seq(-max(abs(rng)), max(abs(rng)), length=40)
#at.vals <- seq(rng[1], rng[2], length=19)
colkey <- list(at=at.vals, labels=list(labels=round(at.vals), at=at.vals))
mapDir <- "../plots/doy_tdd_qmap/maps"
png(paste0(mapDir, "/map_deltas_doytdd10pct_", gsub(" ", "", region), "_historical_", stat, ".png"), width=3200, height=1600, res=200)
p <- levelplot(s, maxpixels=ncell(s), main=paste(region, "DOY TDD 10% deltas for historical", mods[2], "original and quantile-mapped", stat2), na.col="black", par.settings=revRasterTheme, contour=T, margin=F, at=at.vals, colorkey=colkey) +
    layer(sp.polygons(eco_shp, col='gray70'))
print(p)
dev.off()
}
}

# test plot: single pixel
ind <- 100000
x0 <- as.numeric(doytdd0[[2]][ind])
x1 <- as.numeric(doytdd1[[2]][ind])
x1m <- as.numeric(doytdd1.mapped[[2]][ind])
d <- data.table(Source=rep(c("NARR", "GCM", "QMAP"), times=c(length(x0), length(x1), length(x1m))), `DOY TDD 10%`=c(x0, x1, x1m))

ggplot(d, aes(x=`DOY TDD 10%`, fill=Source, colour=Source)) + geom_line(stat="density", size=1) +
    xlim(125,185) + scale_colour_manual(values=clrs) + theme(legend.position="bottom") +
    labs(title="Sample pixel time series distribution\nNARR 1979-2010 and GCM 1958-2005") 

# test plot: multiple pixels
set.seed(29)
n <- 9
ind <- which(!is.na(subset(doytdd0[[2]], 1)[]))
ind <- sample(ind, n)
x0 <- as.numeric(doytdd0[[2]][ind])
x1 <- as.numeric(doytdd1[[2]][ind])
x1m <- as.numeric(doytdd1.mapped[[2]][ind])
d <- data.table(
    Source=factor(rep(c("NARR", "GCM", "QMAP"), times=n*c(length(x0), length(x1), length(x1m))), levels=c("NARR", "GCM", "QMAP")),
    `DOY TDD 10%`=c(x0, x1, x1m), Obs=1:n)

ggplot(d, aes(x=`DOY TDD 10%`, group=interaction(Source, Obs), fill=Source, colour=Source)) + geom_line(stat="density", size=1) +
    xlim(125,185) + scale_colour_manual(values=clrs) + theme(legend.position="bottom") + facet_wrap(~ Obs, ncol=3, scales="free") +
    labs(title="Sample pixel time series distribution\nNARR 1979-2010 and GCM 1958-2005") 

