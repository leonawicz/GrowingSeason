setwd("C:/github/GrowingSeason/workspaces")
pkgs <- list("rasterVis", "maptools", "ggplot2", "data.table", "dplyr", "tidyr")
lapply(pkgs, function(x) library(x, character.only=T))

load("data.RData") # d, d.stats, d.stats2, d.hm, sos, ecomask, yrs, cbpal
load("gbm_preds_eval.RData") # ri.out, cv.out, d.out, s1, sMean, sMean2002
shpDir <- "C:/github/DataExtraction/data/shapefiles"
eco_shp <- shapefile(file.path(shpDir, "AK_ecoregions/akecoregions.shp")) %>% spTransform(CRS(projection(sos))) %>%
  unionSpatialPolygons(eco_shp@data$LEVEL_2)
sw_shp <- shapefile(file.path(shpDir, "Political/Alaska.shp")) %>% spTransform(CRS(projection(sos)))
dir.create(plotDir <- file.path("../plots/gbm/models"), recursive=T, showWarnings=F)

# Spatial Summaries
# CV optimal trees distribution boxplots
png(file.path(plotDir, paste0("gbm_cvbi_byRegion.png")), width=3200, height=1600, res=200)
ggplot(cv.out, aes(x=Region, y=CV)) + geom_boxplot(fill="gray") + geom_point(position=position_jitter(width=0.5)) +
  theme_bw() + theme(legend.position="bottom") + labs(title="5-fold CV optimal trees distributions", y="Number of trees")
dev.off()

# Relative influence by year and region
png(file.path(plotDir, paste0("gbm_RI_byRegion.png")), width=1600, height=1600, res=200)
ggplot(ri.out %>% group_by(Region, Year, Predictor) %>% summarise(RI=mean(RI)), aes(factor(Year), RI, fill=Predictor)) + geom_bar(stat="identity", position="stack") +
  scale_fill_manual(values=cbpal[-2]) + theme_bw(base_size=10) + theme(legend.position="bottom") + coord_flip() +
  ggtitle(paste("Predictor relative influence on start of growing season")) +
  facet_wrap(~Region, ncol=5)
dev.off()

# gbm observed vs. fitted values
png(file.path(plotDir, paste0("gbm_ObsVsFitted_byRegion.png")), width=2000, height=2000, res=200)
ggplot(filter(d.out, is.na(Run)) %>% dcast(Region + Year + Run ~ Source, value.var="SOS"), aes(x=Predicted, y=Observed, colour=Region)) +
  geom_point(position=position_jitter(), size=2) + geom_smooth(method="lm", se=FALSE) + scale_color_manual(values=brewer.pal(9, "Set1")) +
  theme_bw() + theme(legend.position="bottom") + labs(title="Observed vs. fitted values", x="Predicted", y="Observed") #+
  #facet_wrap(~Region, ncol=3)
dev.off()

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
  geom_point(data=d.preds %>% select(-Predicted) %>% distinct, aes(y=SOS, colour=NULL, group=NULL), colour="black", size=3) +
  geom_point(data=d.preds, aes(y=Predicted, colour=NULL, group=NULL), colour="red", size=1, position=position_jitter(width=0.2, height=0)) +
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
