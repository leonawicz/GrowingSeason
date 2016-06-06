setwd("C:/github/GrowingSeason/workspaces")
pkgs <- list("rasterVis", "maptools", "ggplot2", "data.table", "dplyr", "tidyr")
lapply(pkgs, function(x) library(x, character.only=T))

load("data.RData") # d, d.stats, d.stats2, d.hm, sos, ecomask, yrs, cbpal
load("gbm_preds_eval_final.RData") # ri.out, cv.out, pd.out, d.out, s1, sMean, sMean2002
shpDir <- "C:/github/DataExtraction/data/shapefiles"
eco_shp <- shapefile(file.path(shpDir, "AK_ecoregions/akecoregions.shp")) %>% spTransform(CRS(projection(sos)))
eco_shp <- unionSpatialPolygons(eco_shp, eco_shp@data$LEVEL_2)
sw_shp <- shapefile(file.path(shpDir, "Political/Alaska.shp")) %>% spTransform(CRS(projection(sos)))
dir.create(plotDir <- file.path("../plots/gbm/final"), recursive=T, showWarnings=F)

# Spatial Summaries
# CV optimal trees distribution boxplots
png(file.path(plotDir, paste0("gbm_cvbi_byRegion.png")), width=3200, height=1600, res=200)
ggplot(cv.out, aes(x=Region, y=CV)) + geom_boxplot(fill="gray") + geom_point(position=position_jitter(width=0.5)) +
  theme_bw() + theme(legend.position="bottom") + labs(title="5-fold CV optimal trees distributions", y="Number of trees")
dev.off()

# Relative influence by year and region
png(file.path(plotDir, paste0("gbm_RI_byRegion.png")), width=1600, height=1600, res=200)
ggplot(ri.out %>% group_by(Region, Predictor) %>% summarise(RI=mean(RI)), aes(Region, RI, fill=Predictor)) + geom_bar(stat="identity", position="stack") +
  scale_fill_manual(values=cbpal[-2]) + theme_bw(base_size=10) + theme(legend.position="bottom") + coord_flip() +
  ggtitle(paste("Predictor relative influence on start of growing season")) #+ facet_wrap(~Region, ncol=5)
dev.off()

# gbm observed vs. fitted values
png(file.path(plotDir, paste0("gbm_ObsVsFitted_byRegion.png")), width=2000, height=2000, res=200)
ggplot(filter(d.out, is.na(Run)) %>% dcast(Region + Year + Run ~ Source, value.var="SOS"), aes(x=Predicted, y=Observed, colour=Region)) +
  geom_point(position=position_jitter(), size=2) + geom_smooth(method="lm", se=FALSE) + scale_color_manual(values=brewer.pal(9, "Set1")) +
  theme_bw() + theme(legend.position="bottom") + labs(title="Observed vs. fitted values", x="Predicted", y="Observed") #+
  #facet_wrap(~Region, ncol=3)
dev.off()

# partial dependence
pd <- pd.out %>% mutate(Var2=substr(Var, 1, 16)) %>% group_by(Region, Var2) %>% mutate(Ymin=min(Prob))
pd.mean <- summarise(pd, Mean_RI=round(mean(as.numeric(substr(Var, 17, nchar(Var)))), 1))
pd <- left_join(pd, pd.mean)
pd <- group_by(pd) %>% mutate(Var=paste0(Var2, Mean_RI)) %>% select(-Var2, -Mean_RI)

pd_merge_data <- function(x, n=1000){
  f <- function(x, rv, rx, n){
    vp <- approx(x$Val, x$Prob, xout=seq(rv[1], rv[2], length=n))
    xy <- approx(x$x, x$y, xout=seq(rx[1], rx[2], length=n))
    data.frame(Run=x$Run[1], Val=vp$x, Prob=vp$y, x=xy$x, y=xy$y)
  }
  rv <- range(x$Val)
  rx <- range(x$x)
  x0 <- select_(x, .dots=list("-Val", "-Prob", "-x", "-y")) %>% distinct
  x1 <- x %>% split(.$Run) %>% purrr::map(~f(.x, rv, rx, n)) %>% bind_rows
  left_join(x0, x1)
}  

pd.interp <- pd %>% split(paste(.$Region, .$Var)) %>% purrr::map(~pd_merge_data(.x)) %>% bind_rows
pd.ci <- pd.interp %>% group_by(Region, Var, Ymin, Val, x) %>% summarise(
  n_na=length(which(is.na(y))),
  Prob_mean=mean(Prob, na.rm=TRUE),
  Prob_025=quantile(Prob, 0.025, na.rm=TRUE),
  Prob_975=quantile(Prob, 0.975, na.rm=TRUE),
  y_mean=mean(y, na.rm=TRUE),
  y_025=ifelse(n_na > 0, NA, quantile(y, 0.025, na.rm=TRUE)),
  y_975=ifelse(n_na > 0, NA, quantile(y, 0.975, na.rm=TRUE)))
  
pd.ci <- pd.interp %>% group_by(Region, Var, Ymin, Val, x) %>% summarise(
  n_na=length(which(is.na(y))),
  Prob_mean=mean(Prob, na.rm=TRUE),
  Prob_min=min(Prob, na.rm=TRUE),
  Prob_max=max(Prob, na.rm=TRUE),
  y_mean=ifelse(n_na > 16, NA, mean(y, na.rm=TRUE)),
  y_min=ifelse(n_na > 16, NA, quantile(y, 0.025, na.rm=TRUE)),
  y_max=ifelse(n_na > 16, NA, quantile(y, 0.975, na.rm=TRUE)))

# bands
save_pdplot <- function(x, outDir){
  region <- gsub(" ", "", unique(x$Region))
  g <- ggplot(x, aes(x=Val)) + facet_wrap(~Var, ncol=2, scales="free") + labs(x="x", y="y")
  g <- g + geom_ribbon(aes(ymin=Ymin, ymax=Prob_mean), fill="orange") +
  geom_ribbon(aes(x=x, ymin=y_min, ymax=y_max), fill="black", alpha=0.2) +
  #geom_line(aes(x=x, y=y_mean), size=1) +
  geom_line(aes(x=x, y=y_min), size=1) +
  geom_line(aes(x=x, y=y_max), size=1)
  png(file.path(outDir, paste0("gbm_PD_", region, ".png")), width=2000, height=2000, res=200)
  print(g)
  dev.off()
}

pd.ci %>% split(.$Region) %>% purrr::walk(~save_pdplot(.x, outDir=plotDir))

# spaghetti
save_pdplot <- function(x, outDir){
  region <- gsub(" ", "", unique(x$Region))
  g <- ggplot(x, aes(x=Val)) + facet_wrap(~Var, ncol=2, scales="free") + labs(x="x", y="y")
  g <- g + geom_ribbon(aes(ymin=Ymin, ymax=Prob, group=Run), fill="orange", alpha=0.2) + geom_line(aes(x=x, y=y, group=Run), size=1, alpha=0.2)
  png(file.path(outDir, paste0("gbm_PD_", region, ".png")), width=2000, height=2000, res=200)
  print(g)
  dev.off()
}

pd %>% split(.$Region) %>% purrr::walk(~save_pdplot(.x, outDir=plotDir))

#################################

d.hoy <- d.preds %>% select(-SOS) %>% mutate(Source="Predicted HOY") %>% rename(SOS=Predicted)
d.hoy <- d.preds %>% select(-Predicted) %>% distinct %>% mutate(Source="Global observed") %>% bind_rows(d.hoy)
d.ts <- bind_rows(filter(d.out, Source!="Bias corrected"), d.hoy) %>% mutate(Source=factor(Source, levels=c("Observed", "Predicted", "Global observed", "Predicted HOY")))

extract_to_dt <- function(x, y, fun, rcp, gcm, run, ...){
  raster::extract(x, y, fun, ...) %>% t %>% data.table %>% setnames(names(y)) %>% mutate(Run=run) %>% melt(id.vars="Run", variable.name="Region", value.name="SOS") %>%
    mutate(Year=as.integer(substr(names(x), 5, 8)), RCP=factor(rcp, levels=c("RCP 6.0", "RCP 8.5")), Model=gcm, Source="Projected") %>% select(RCP, Model, Region, Year, Source, SOS, Run)
}

load("sos_gbm_preds_GFDL-CM3_test.RData")
d.proj <- bind_rows(extract_to_dt(b.rcp60.qm, eco_shp, mean, "RCP 6.0", "GFDL-CM3", run=1, na.rm=TRUE), extract_to_dt(b.rcp85.qm, eco_shp, mean, "RCP 8.5", "GFDL-CM3", run=2, na.rm=TRUE))
load("sos_gbm_preds_IPSL-CM5A-LR.RData")
d.proj <- bind_rows(d.proj, extract_to_dt(b.rcp60.qm, eco_shp, mean, "RCP 6.0", "IPSL-CM5A-LR", run=3, na.rm=TRUE), extract_to_dt(b.rcp85.qm, eco_shp, mean, "RCP 8.5", "IPSL-CM5A-LR", run=4, na.rm=TRUE))
load("sos_gbm_preds_MRI-CGCM3.RData")
d.proj <- bind_rows(d.proj, extract_to_dt(b.rcp60.qm, eco_shp, mean, "RCP 6.0", "MRI-CGCM3", run=5, na.rm=TRUE), extract_to_dt(b.rcp85.qm, eco_shp, mean, "RCP 8.5", "MRI-CGCM3", run=6, na.rm=TRUE))
d.proj.mean <- d.proj %>% group_by(Region, Year, Source) %>% summarise(SOS=mean(SOS), Run=1)
d.smooth <- filter(d.ts, is.na(Run) & Source=="Predicted") %>% bind_rows(filter(d.proj, Year > 2010)) %>% group_by(Region, Year) %>% summarise(SOS=mean(SOS), Source="Trend", Run=1)

load("old_sos_preds/sos_gbm_preds_GFDL-CM3.RData")
d.proj.old <- bind_rows(extract_to_dt(b.rcp60.qm, eco_shp, mean, "RCP 6.0", "GFDL-CM3", run=1, na.rm=TRUE), extract_to_dt(b.rcp85.qm, eco_shp, mean, "RCP 8.5", "GFDL-CM3", run=2, na.rm=TRUE))
load("old_sos_preds/sos_gbm_preds_IPSL-CM5A-LR.RData")
d.proj.old <- bind_rows(d.proj.old, extract_to_dt(b.rcp60.qm, eco_shp, mean, "RCP 6.0", "IPSL-CM5A-LR", run=3, na.rm=TRUE), extract_to_dt(b.rcp85.qm, eco_shp, mean, "RCP 8.5", "IPSL-CM5A-LR", run=4, na.rm=TRUE))
load("old_sos_preds/sos_gbm_preds_MRI-CGCM3.RData")
d.proj.old <- bind_rows(d.proj.old, extract_to_dt(b.rcp60.qm, eco_shp, mean, "RCP 6.0", "MRI-CGCM3", run=5, na.rm=TRUE), extract_to_dt(b.rcp85.qm, eco_shp, mean, "RCP 8.5", "MRI-CGCM3", run=6, na.rm=TRUE))
d.proj.old.mean <- d.proj.old %>% group_by(Region, Year, Source) %>% summarise(SOS=mean(SOS), Run=1)
d.smooth.old <- filter(d.ts, is.na(Run) & Source=="Predicted") %>% bind_rows(filter(d.proj.old, Year > 2010)) %>% group_by(Region, Year) %>% summarise(SOS=mean(SOS), Source="Trend", Run=1)

z <- b.rcp60.qm %>% subset(5:94) %>% stackApply(rep(1:9, each=10), mean) %>% subset(c(2,9)) %>% calc(function(x) x[2]-x[1])
rTheme <- function(region=colorRampPalette(c("darkred", "firebrick1", "white", "royalblue", "darkblue"))(19), ...){
  theme <- custom.theme.2(region=region, ...)
  theme$strip.background$col <- theme$strip.shingle$col <- theme$strip.border$col <- theme$panel.background$col <- theme$axis.line$col <- "transparent"
  theme
}

rng <- range(z[], na.rm=T)
rng <- range(rng, -rng)
at.vals <- seq(rng[1], rng[2], length=19)
colkey <- list(at=at.vals, labels=list(labels=round(at.vals), at=at.vals))

#png(file.path(plotDir, paste0(".png")), width=3200, height=1600, res=200)
p <- levelplot(z, maxpixels=ncell(z), main=paste("2090s - 2010s deltas"),
  par.settings=rTheme, xlab=NULL, ylab=NULL, scales=list(draw=F), contour=F, margin=F, at=at.vals, colorkey=colkey) +
  latticeExtra::layer(sp.polygons(sw_shp, fill='transparent', alpha=0.3))
print(p)
#dev.off()

# historical predictions over observations time series
clrs <- c("peru", "royalblue", "black", "red")
png(file.path(plotDir, paste0("gbm_TSpreds_byRegion.png")), width=3200, height=1600, res=200)
ggplot(filter(d.out, !is.na(Run) & Source!="Bias corrected"), aes(x=Year, y=SOS, colour=Source, group=interaction(Source, Run))) +
  scale_color_manual(values=clrs) + 
  geom_line(size=1, alpha=0.25) +
  #geom_line(data=filter(d.out, Source=="Bias corrected" & is.na(Run)), colour="#B8860B", size=1, linetype=1) +
  geom_line(data=filter(d.out, Source=="Observed" & is.na(Run)), colour="black", size=2) + 
  geom_line(data=filter(d.out, Source=="Predicted" & is.na(Run)), colour="black", size=2) +
  geom_line(data=filter(d.out, Source=="Predicted" & is.na(Run)), colour=clrs[2], size=1) +
  geom_line(data=filter(d.out, Source=="Observed" & is.na(Run)), colour=clrs[1], size=1) + 
  geom_line(data=filter(d.out, Source=="Predicted" & is.na(Run)), colour=clrs[2], size=1, linetype=2) +
  theme_bw() + theme(legend.position="bottom") + ggtitle("Observed and modeled start of growing season") +
  #scale_x_continuous(breaks=c(1982,1990,2000,2010)) +
  guides(fill=guide_legend(override.aes=list(alpha=1)), colour=guide_legend(override.aes=list(alpha=1))) +
  facet_wrap(~Region, ncol=3, scales="free")
dev.off()

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
  geom_line(data=d.proj.old, colour="#00000030", size=1) + 
  geom_smooth(data=d.smooth.old) +
  theme_bw() + theme(legend.position="bottom") + ggtitle("Observed and modeled start of growing season") +
  #scale_x_continuous(breaks=c(1982,1990,2000,2010)) +
  guides(fill=guide_legend(override.aes=list(alpha=1)), colour=guide_legend(override.aes=list(alpha=1))) +
  facet_wrap(~Region, ncol=3, scales="free")
dev.off()

# historical predictions over observations time series, plus predictions on withheld years
png(file.path(plotDir, paste0("gbm_TSpreds_byRegion_HOY.png")), width=3200, height=1600, res=200)
ggplot(filter(d.ts, !is.na(Run) & !(Source %in% c("Bias corrected", "Global observed", "Predicted HOY"))), aes(x=Year, y=SOS, colour=Source, group=interaction(Source, Run))) +
  scale_color_manual(values=clrs[c(3,1,2,4)]) + 
  geom_line(size=1, alpha=0.25) +
  #geom_line(data=filter(d.out, Source=="Bias corrected" & is.na(Run)), colour="#B8860B", size=1, linetype=1) +
  geom_line(data=filter(d.out, Source=="Observed" & is.na(Run)), colour="black", size=2) + 
  geom_line(data=filter(d.out, Source=="Predicted" & is.na(Run)), colour="black", size=2) +
  geom_line(data=filter(d.out, Source=="Predicted" & is.na(Run)), colour=clrs[2], size=1) +
  geom_line(data=filter(d.out, Source=="Observed" & is.na(Run)), colour=clrs[1], size=1) + 
  geom_line(data=filter(d.out, Source=="Predicted" & is.na(Run)), colour=clrs[2], size=1, linetype=2) +
  geom_point(data=d.ts %>% filter(Source=="Global observed"), aes(group=NULL), size=3) +
  geom_point(data=d.ts %>% filter(Source=="Predicted HOY"), size=1, position=position_jitter(width=0.2, height=0)) +
  theme_bw() + theme(legend.position="bottom") + ggtitle("Observed and modeled start of growing season") +
  scale_x_continuous(breaks=c(1982,1990,2000,2010)) +
  guides(fill=guide_legend(override.aes=list(alpha=1)), colour=guide_legend(override.aes=list(alpha=1))) +
  facet_wrap(~Region, ncol=3, scales="free")
dev.off()

d.sd <- filter(d.out, Source=="Observed" | (!is.na(Run) & Source!="Bias corrected")) %>% group_by(Region, Source) %>% summarise(SD=sd(SOS))
d.hoy.sd <- d.preds %>% group_by(Region) %>% summarise(`Global observered`=sd(unique(SOS)), `Predicted HOY`=sd(Predicted)) %>% melt(id.vars="Region", variable.name="Source", value.name="SD")
d.sd <- bind_rows(d.sd, d.hoy.sd)

# predictions and observations historical inter-annual variability
png(file.path(plotDir, paste0("gbm_TSpreds_byRegionSD.png")), width=3200, height=1600, res=200)
ggplot(d.sd, aes(x=Region, y=SD, fill=Source)) +
  scale_fill_manual(values=clrs) + geom_bar(stat="identity", position="dodge") +
  theme_bw() + theme(legend.position="bottom") + ggtitle("Observed and modeled start of growing season historical inter-annual SD")
dev.off()

# Maps
# Settings
rTheme <- function(region=colorRampPalette(c("darkred", "firebrick1", "white", "royalblue", "darkblue"))(19), ...){
  theme <- custom.theme.2(region=region, ...)
  theme$strip.background$col <- theme$strip.shingle$col <- theme$strip.border$col <- theme$panel.background$col <- theme$axis.line$col <- "transparent"
  theme
}

rng <- range(c(s1[], sMean[]), na.rm=T)
rng <- range(rng, -rng)
at.vals <- seq(rng[1], rng[2], length=19)
colkey <- list(at=at.vals, labels=list(labels=round(at.vals), at=at.vals))

# historical mean prediction deltas map for one example gbm model
png(file.path(plotDir, paste0("gbm_pred_deltas_gbm1.png")), width=3200, height=1600, res=200)
p <- levelplot(s1, maxpixels=ncell(s1), main=paste("1982-2010 mean start of growing season deltas: single GBM"),
  par.settings=rTheme, xlab=NULL, ylab=NULL, scales=list(draw=F), contour=F, margin=F, at=at.vals, colorkey=colkey) +
  latticeExtra::layer(sp.polygons(sw_shp, fill='transparent', alpha=0.3))
print(p)
dev.off()

# historical mean prediction deltas map averaged over multiple gbm models
png(file.path(plotDir, paste0("gbm_pred_deltas_gbmMean.png")), width=3200, height=1600, res=200)
p <- levelplot(sMean, maxpixels=ncell(sMean), main=paste("1982-2010 mean start of growing season deltas: mean GBM"),
  par.settings=rTheme, xlab=NULL, ylab=NULL, scales=list(draw=F), contour=F, margin=F, at=at.vals, colorkey=colkey) +
  latticeExtra::layer(sp.polygons(sw_shp, fill='transparent', alpha=0.3))
print(p)
dev.off()

rng <- range(sMean2002[], na.rm=T)
rng <- range(rng, -rng)
at.vals <- seq(rng[1], rng[2], length=19)
colkey <- list(at=at.vals, labels=list(labels=round(at.vals), at=at.vals))

# prediction deltas map for 2002 averaged over multiple gbm models
png(file.path(plotDir, paste0("gbm_pred_deltas_gbmMean2002.png")), width=3200, height=1600, res=200)
p <- levelplot(sMean2002, maxpixels=ncell(sMean2002), main=paste("2002 start of growing season deltas: mean GBM"),
  par.settings=rTheme, xlab=NULL, ylab=NULL, scales=list(draw=F), contour=F, margin=F, at=at.vals, colorkey=colkey) +
  latticeExtra::layer(sp.polygons(sw_shp, fill='transparent', alpha=0.3))
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

# historical mean predictions maps by multiple gbm models for comparison
png(file.path(plotDir, paste0("gbm_pred_byModel.png")), width=3200, height=1600, res=200)
p <- levelplot(pred.maps.gbm.samples, maxpixels=ncell(pred.maps.gbm.samples), main=paste("1982-2010 mean start of growing season GBM comparisons"),
  par.settings=rTheme, xlab=NULL, ylab=NULL, scales=list(draw=F), contour=F, margin=F) +
  latticeExtra::layer(sp.polygons(sw_shp, fill='transparent', alpha=0.3))
print(p)
dev.off()

sos2002 <- subset(sos, 21)
s <- stack(sos2002, sos2002 + sMean2002)
names(s) <- c("Observed", "Predicted")

# predictions map for 2002 averaged over multiple gbm models
png(file.path(plotDir, paste0("gbm_pred_gbmMean2002.png")), width=3200, height=1600, res=200)
p <- levelplot(s, maxpixels=ncell(s), main=paste("2002 start of growing season: mean GBM"),
  par.settings=rTheme, xlab=NULL, ylab=NULL, scales=list(draw=F), contour=F, margin=F) +
  latticeExtra::layer(sp.polygons(sw_shp, fill='transparent', alpha=0.3))
print(p)
dev.off()

s <- crop(s, arctic) %>% mask(arctic)

# predictions map for 2002 averaged over multiple gbm models masked to Arctic tundra ecoregion
png(file.path(plotDir, paste0("gbm_pred_gbmMean2002_ArcticTundra.png")), width=3200, height=1600, res=200)
p <- levelplot(s, maxpixels=ncell(s), main=paste("2002 start of Arctic Tundra growing season: mean GBM"),
  par.settings=rTheme, xlab=NULL, ylab=NULL, scales=list(draw=F), contour=F, margin=F) +
  latticeExtra::layer(sp.polygons(sw_shp, fill='transparent', alpha=0.3))
print(p)
dev.off()
