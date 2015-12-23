setwd("/atlas_scratch/mfleonawicz/projects/GrowingSeason/workspaces")

library(maptools)
library(raster)
library(data.table)
library(dplyr)
library(parallel)

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

myFun <- function(x, ...) c(Mean=mean(x,...), SD=sd(x,...), quantile(x, prob=c(0.05, 0.95),...))

get_brick_stats <- function(i, x, shp, years, scenario, model, map, ...){
    stopifnot(any(class(x[[i]]) %in% c("RasterBrick", "RasterStack")))
    nam <- factor(names(x)[i], levels=names(x))
    ids <- sapply(slot(shp, "polygons"), function(x) slot(x, "ID"))
    x <- extract(x[[i]], shp, ...)
    if(is.list(x)) stop("Unequal length results returned by polygon. Check summary function.")
    n <- nrow(x)/length(ids)
    x <- lapply(1:n, function(i, x, rnames) x[rownames(x)==unique(rownames(x))[i],], x=x)
    x <- lapply(x, function(x, years, nam) data.table(Year=years, Pct=nam, Stat=unique(rownames(x)), t(x)), years=years, nam=nam)
    x <- rbindlist(x)
    setnames(x, c("Year", "Pct", "Stat", ids))
    x <- melt(x, id.vars=c("Year", "Pct", "Stat"), variable.name="Region", value.name="DOY TDD")
    x <- mutate(x, Scenario=scenario, Model=model, Map=map)
    setcolorder(x, c("Scenario", "Model", "Map", "Region", "Pct", "Year", "Stat", "DOY TDD"))
    x
}

get_stat_maps <- function(i, x, ...) calc(x[[i]], myFun, ...)

system.time({
maplist <- vector("list", length(models))
for(i in 1:length(models)){
    model <- models[i]
    load(files[i])
    n <- length(doytdd0)
    names(doytdd0) <- names(doytdd1.mapped) <- names(doytdd.rcp60.mapped) <- names(doytdd.rcp85.mapped) <- names(doytdd1)
    if(i==1) d <- rbindlist(mclapply(1:n, get_brick_stats, x=doytdd0, shp=eco_shp, years=yrs.narr, scenario="Historical", model="NARR", map="Baseline", fun=myFun, na.rm=TRUE, mc.cores=n))
    d <- bind_rows(d, rbindlist(mclapply(1:n, get_brick_stats, x=doytdd1, shp=eco_shp, years=yrs.hist, scenario="Historical", model=model, map="Original", fun=myFun, na.rm=TRUE, mc.cores=n)))
    d <- bind_rows(d, rbindlist(mclapply(1:n, get_brick_stats, x=doytdd1.mapped, shp=eco_shp, years=yrs.hist, scenario="Historical", model=model, map="Quantile-mapped", fun=myFun, na.rm=TRUE, mc.cores=n)))
    d <- bind_rows(d, rbindlist(mclapply(1:n, get_brick_stats, x=doytdd.rcp60, shp=eco_shp, years=yrs.proj, scenario="RCP 6.0", model=model, map="Original", fun=myFun, na.rm=TRUE, mc.cores=n)))
    d <- bind_rows(d, rbindlist(mclapply(1:n, get_brick_stats, x=doytdd.rcp60.mapped, shp=eco_shp, years=yrs.proj, scenario="RCP 6.0", model=model, map="Quantile-mapped", fun=myFun, na.rm=TRUE, mc.cores=n)))
    d <- bind_rows(d, rbindlist(mclapply(1:n, get_brick_stats, x=doytdd.rcp85, shp=eco_shp, years=yrs.proj, scenario="RCP 8.5", model=model, map="Original", fun=myFun, na.rm=TRUE, mc.cores=n)))
    d <- bind_rows(d, rbindlist(mclapply(1:n, get_brick_stats, x=doytdd.rcp85.mapped, shp=eco_shp, years=yrs.proj, scenario="RCP 8.5", model=model, map="Quantile-mapped", fun=myFun, na.rm=TRUE, mc.cores=n)))
    
    if(i==1) x0 <- mclapply(1:n, get_stat_maps, x=doytdd0, na.rm=TRUE, mc.cores=n)
    x1 <- mclapply(1:n, get_stat_maps, x=doytdd1, na.rm=TRUE, mc.cores=n)
    x1m <- mclapply(1:n, get_stat_maps, x=doytdd1.mapped, na.rm=TRUE, mc.cores=n)
    x1.deltas <- lapply(1:length(x0), function(i, x, y) y[[i]]-x[[i]], x=x0, y=x1)
    x1m.deltas <- lapply(1:length(x0), function(i, x, y) y[[i]]-x[[i]], x=x0, y=x1m)
    names(x1.deltas) <- names(x1m.deltas) <- names(doytdd0)
    maplist[[i]] <- list(diffx1x0=x1.deltas, diffx1mx0=x1m.deltas)
}
})

names(maplist) <- models

save(d, maplist, ecomask, eco_shp, file="doytdd_final.RData")
###################################################################
setwd("C:/github/GrowingSeason/workspaces")
library(data.table)
library(dplyr)
library(ggplot2)
library(rasterVis)

load("doytdd_final.RData")

clrs <- c("black", "darkorange", "purple")
scens <- c("Historical", "RCP 6.0")
mods <- c("NARR", "MRI-CGCM3")
spstat <- "SD"
regions <- names(eco_shp)

narr.means <- d %>% filter(Model=="NARR") %>% group_by(Region, Pct, Stat) %>% summarise(NARR_Mean=mean(`DOY TDD`))
d2 <- left_join(d, narr.means) %>% mutate(`DOY TDD Deltas`=`DOY TDD`-NARR_Mean)

# save workspace for shiny app
d.qmap <- d2
save(d.qmap, eco_shp, ecomask, regions, maplist, file="C:/github/shiny-apps/sos/appdata_qmap.RData")

dsub <- d %>% filter(Scenario %in% scens & Model %in% mods)
dsub2 <- d2 %>% filter(Scenario %in% scens & Model %in% mods)

# Time series and histogram DOY TDD and deltas
tsDir <- "../plots/doy_tdd_qmap/timeseries"
denDir <- "../plots/doy_tdd_qmap/distributions"
png(paste0(tsDir,"/ts_full_", spstat, "_doytdd10pct_byRegion.png"), width=3200, height=2400, res=200)
ggplot(dsub %>% filter(Pct=="0.1" & Stat==spstat), aes(x=Year, y=`DOY TDD`, group=interaction(Model, Map), colour=Map)) +
    geom_line(data=dsub %>% filter(Pct=="0.1" & Stat==spstat & Model=="NARR"), size=1) + geom_line() +
    scale_colour_manual(values=clrs) + theme(legend.position="bottom") + facet_wrap(~ Region, scales="free_y") +
    labs(y="DOY TDD 10%", title=paste("NARR baseline and", mods[2], scens[scens!="Historical"][1], "original and quantile-mapped regional", spstat, "DOY TDD 10%"))
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

