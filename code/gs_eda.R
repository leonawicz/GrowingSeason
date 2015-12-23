#############################################
#### Growing season thaw degree days EDA ####
#############################################

#### Script author:  Matthew Leonawicz ####
#### Maintainted by: Matthew Leonawicz ####
#### Last updated:   06/25/2015        ####

# @knitr setup
setwd("C:/github/GrowingSeason/workspaces")
plotDir <- "../plots"

library(rasterVis)
library(maptools)
library(ggplot2)
library(data.table)
library(dplyr)

yrs <- 1982:2010
cbpal <- c("#8B4500", "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Start of season map
sos <- brick("../data/sos_1982_2010.tif")

# Threshold thaw degree days maps
tdd05 <- dropLayer(brick("../data/pct05_tdd_spring_1979_2010.tif"), 1:3)
tdd10 <- dropLayer(brick("../data/pct10_tdd_spring_1979_2010.tif"), 1:3)
tdd15 <- dropLayer(brick("../data/pct15_tdd_spring_1979_2010.tif"), 1:3)
tdd20 <- dropLayer(brick("../data/pct20_tdd_spring_1979_2010.tif"), 1:3)

# Alaska ecoregion level two mask
shpDir <- "C:/github/DataExtraction/data/shapefiles"
eco_shp <- shapefile(file.path(shpDir, "AK_ecoregions/akecoregions.shp"))
eco_shp <- spTransform(eco_shp, CRS(projection(sos)))
eco_shp <- unionSpatialPolygons(eco_shp, eco_shp@data$LEVEL_2)
eco_IDs <- sapply(slot(eco_shp, "polygons"), function(x) slot(x, "ID"))
ecomask <- rasterize(eco_shp, sos)

# @knitr data_prep
make_TDD_df <- function(d, extractBy, y){
    d[d <= 1] <- NA
    d %>% projectRaster(y) %>% resample(y, method="bilinear") %>% mask(y) %>% extract(extractBy, cellnumbers=TRUE) -> d
    s <- extract(y, extractBy)
    d <- rbindlist(lapply(1:length(d),
        function(i, tdd, sos, years, eco, shp, mask){
            xy <- xyFromCell(mask, as.numeric(tdd[[i]][,1]))
            data.table(Region=eco[i], Year=rep(years, each=nrow(tdd[[i]])), SOS=as.numeric(sos[[i]]), TDD=as.numeric(tdd[[i]][,-1]), x=xy[,1], y=xy[,2])
        }, tdd=d, sos=s, years=yrs, eco=eco_IDs, shp=extractBy, mask=y))
    d[!is.na(SOS) & !is.na(TDD),]
}

d <- lapply(list(tdd05, tdd10, tdd15, tdd20), make_TDD_df, extractBy=eco_shp, y=sos)
lapply(1:length(d), function(i, x, pct) x[[i]][, Threshold := pct[i]], x=d, pct=paste0(c("05",10,15,20), "pct"))
d <- rbindlist(d)
tmp <- with(d, paste(Region, Threshold, Year))
tmp <- as.numeric(lapply(split(tmp, tmp), length))
setorder(d, Region, Threshold, Year)
d[, Obs := unlist(lapply(tmp, function(x) 1:x))]
setcolorder(d, c("Region", "Year", "Threshold", "Obs", "x", "y", "SOS", "TDD"))

d %>% group_by(Threshold, Region) %>%
    summarise(Min=min(TDD, na.rm=T),
              Pct05=quantile(TDD, 0.05, na.rm=T),
              Pct10=quantile(TDD, 0.10, na.rm=T),
              Pct25=quantile(TDD, 0.25, na.rm=T),
              Pct50=quantile(TDD, 0.5, na.rm=T),
              Pct75=quantile(TDD, 0.75, na.rm=T),
              Pct90=quantile(TDD, 0.90, na.rm=T),
              Pct95=quantile(TDD, 0.95, na.rm=T),
              Max=max(TDD, na.rm=T), Mean=mean(TDD, na.rm=T), SD=sd(TDD, na.rm=T)) -> d.stats
              
d %>% group_by(Threshold, Year) %>%
    summarise(Min=min(TDD, na.rm=T),
              Pct05=quantile(TDD, 0.05, na.rm=T),
              Pct10=quantile(TDD, 0.10, na.rm=T),
              Pct25=quantile(TDD, 0.25, na.rm=T),
              Pct50=quantile(TDD, 0.5, na.rm=T),
              Pct75=quantile(TDD, 0.75, na.rm=T),
              Pct90=quantile(TDD, 0.90, na.rm=T),
              Pct95=quantile(TDD, 0.95, na.rm=T),
              Max=max(TDD, na.rm=T), Mean=mean(TDD, na.rm=T), SD=sd(TDD, na.rm=T)) -> d.stats2

d %>% group_by(Region, Threshold, Year) %>% summarise(SD=sd(TDD, na.rm=T)) -> d.hm
save(d, d.stats, d.stats2, d.hm, sos, eco_shp, ecomask, yrs, cbpal, file="data.RData")

# @knitr setup2
setwd("C:/github/GrowingSeason/workspaces")
plotDir <- "../plots"

library(rasterVis)
library(maptools)
library(ggplot2)
library(data.table)
library(dplyr)

load("data.RData")

# @knitr plot_setup
# Plot ecoregions
econames <- c("AK Range Trans.", "Aleutian Meadows", "Arctic Tundra", "Bering Taiga", "Bering Tundra", "Coastal Mtn. Trans.", "Coastal Rainforests", "Intermontane Boreal", "Pacific Mtn. Trans.")
n.reg <- length(econames)
at.vals <- 0:n.reg
colkey <- list(at=at.vals, labels=list(labels=econames, at=at.vals + 0.5))

# rasterVis theme settings
revRasterTheme <- function (pch = 19, cex = 0.7, region=cbpal, ...){
    theme <- custom.theme.2(pch = pch, cex = cex, region = region, ...)
    theme$strip.background$col <- theme$strip.shingle$col <- theme$strip.border$col <- "transparent"
    theme$add.line$lwd = 0.4
    theme
}

# @knitr plot_eco
levelplot(ecomask, main="AK level 2 ecoregions", par.settings=revRasterTheme, scales=list(draw=FALSE), contour=F, margin=F, at=at.vals, colorkey=colkey)

# @knitr plot_marginal_tdd_01a
# TDD marginal spatial distributions (across years) by threshold and ecoregion
(p01a <- ggplot(data=d, aes(x=TDD, colour=Threshold)) + theme(legend.position="bottom") +
    geom_line(stat="density", size=1) + facet_wrap(~ Region))

#### SPACE ####
# @knitr tables_marginal_tdd1
# tdd marginal spatial distribution quantiles (across years) by threshold and ecoregion
setorder(d.stats, Threshold, Mean)
subset(d.stats, Threshold=="05pct", -1) %>% knitr::kable(digits=0, caption="5 % Thaw Degree Days.")
subset(d.stats, Threshold=="10pct", -1) %>% knitr::kable(digits=0, caption="10 % Thaw Degree Days.")
subset(d.stats, Threshold=="15pct", -1) %>% knitr::kable(digits=0, caption="15 % Thaw Degree Days.")
subset(d.stats, Threshold=="20pct", -1) %>% knitr::kable(digits=0, caption="20 % Thaw Degree Days.")
#d.stats %>% filter(Threshold=="05pct") %>% arrange(Mean)
#d.stats %>% filter(Threshold=="10pct") %>% arrange(Mean)
#d.stats %>% filter(Threshold=="15pct") %>% arrange(Mean)
#d.stats %>% filter(Threshold=="20pct") %>% arrange(Mean)

#### TIME ####
# @knitr tables_marginal_tdd2
# tdd marginal spatial distribution quantiles (across ecoregions) by threshold and year
subset(d.stats2, Threshold=="05pct", -1) %>% knitr::kable(digits=0, caption="5 % Thaw Degree Days.")
subset(d.stats2, Threshold=="10pct", -1) %>% knitr::kable(digits=0, caption="10 % Thaw Degree Days.")
subset(d.stats2, Threshold=="15pct", -1) %>% knitr::kable(digits=0, caption="15 % Thaw Degree Days.")
subset(d.stats2, Threshold=="20pct", -1) %>% knitr::kable(digits=0, caption="20 % Thaw Degree Days.")
#d.stats %>% filter(Threshold=="05pct") %>% arrange(Mean)
#d.stats %>% filter(Threshold=="10pct") %>% arrange(Mean)
#d.stats %>% filter(Threshold=="15pct") %>% arrange(Mean)
#d.stats %>% filter(Threshold=="20pct") %>% arrange(Mean)

#### SPACE-TIME ####
# @knitr plot_tdd_02a
# Spatial distributions | time
# TDD spatial distributions by ecoregion and year | threshold = X%
dsub <- d[Threshold=="10pct",]
(p02a <- ggplot(data=dsub, aes(x=TDD, group=Year)) + theme(legend.position="bottom") +
    xlab("10% TDD") +
    geom_line(stat="density", alpha=0.5) + facet_wrap(~ Region))

# @knitr plot_tdd_02b
# TDD spatial distributions by ecoregion and year | threshold = X%
(p02b <- ggplot(data=dsub, aes(x=factor(Year), y=TDD)) +
    theme(legend.position="bottom", axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
    xlab("1982 - 2010") + ylab("10% TDD") +
    geom_boxplot(fill="white", outlier.size=NA) + geom_jitter(alpha=0.01) + facet_wrap(~ Region))

# @knitr plot_tdd_02c
# TDD spatial distributions by year | threshold = X% and region = Y
dsub <- d[Region %in% c("Pacific Mountains Transition", "Coastal Rainforests"),]
(p02c <- ggplot(data=dsub, aes(x=factor(Year), y=TDD, fill=Region, colour=Region)) +
    theme(legend.position="bottom", axis.ticks.x = element_blank(), axis.text.x = element_blank()) + guides(colour=guide_legend(override.aes=list(alpha=1))) +
    xlab("1982 - 2010") + ylab("10% TDD") +
    geom_boxplot(fill="white", position=position_dodge(width=0.9), outlier.size=NA) + geom_point(alpha=0.01, position=position_jitterdodge(dodge.width=0.9)) + facet_wrap(~ Threshold))

# @knitr plot_tdd_02d
# Time series by sampled locations | threshold = X% and region = Y
(p02d <- ggplot(data=dsub[Obs <=1000,], aes(x=Year, y=TDD, group=interaction(Region, Threshold, Obs), fill=Region, colour=Region)) +
    theme(legend.position="bottom", axis.ticks.x = element_blank(), axis.text.x = element_blank()) + guides(colour=guide_legend(override.aes=list(alpha=1))) +
    xlab("1982 - 2010") + ylab("10% TDD") +
    geom_line(stat="smooth", alpha=0.1) + facet_wrap(~ Threshold))

# @knitr plot_tdd_02e
# Bivariate space-time heatmap, spatial SD by ecoregion and year
#d %>% filter(Threshold=="10pct") %>% group_by(Region, Year) %>% summarise(SD=sd(TDD, na.rm=T)) %>% group_by() %>% mutate(RelVar=SD/max(SD)) -> dsub
dsub <- subset(d.hm, Threshold=="10pct")
(p02e <- ggplot(dsub, aes(x=factor(Year), y=Region)) +
    theme(legend.position="bottom", axis.ticks.x = element_blank()) +
    xlab("Year") + ylab("Ecoregion") +
    geom_tile(aes(fill=SD), colour = "white") +
    scale_fill_gradient(low = "white", high = "steelblue", limit = c(-0, max(dsub$SD)), name="10% TDD SD") +
    geom_text(aes(factor(Year), Region, label=round(SD, 1)), color="black", size=4))
