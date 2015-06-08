#############################################
#### Growing season thaw degree days EDA ####
#############################################

#### Script author:  Matthew Leonawicz ####
#### Maintainted by: Matthew Leonawicz ####
#### Last updated:   06/03/2015        ####

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
eco_shp <- shapefile(file.path(shpDir, "AK_ecoregions/akecoregions"))
eco_shp <- spTransform(eco_shp, CRS(projection(sos)))
eco_shp <- unionSpatialPolygons(eco_shp, eco_shp@data$LEVEL_2)
eco_IDs <- sapply(slot(eco_shp, "polygons"), function(x) slot(x, "ID"))
ecomask <- rasterize(eco_shp, sos)

# @knitr data_prep
make_TDD_df <- function(d, extractBy, projectTo=NULL, resampleTo=NULL, maskTo=resampleTo){
    d[d <= 1] <- NA
    if(!is.null(projectTo)) d <- projectRaster(d, sos)
    if(!is.null(resampleTo)) d <- resample(d, ecomask, method="bilinear")
    if(!is.null(maskTo)) d <- mask(d, ecomask)
    e <- extract(d, eco_shp)
    d <- rbindlist(lapply(1:length(e),
        function(i, x, years, eco){
            data.table(Region=eco[i], Year=rep(years, each=nrow(x[[i]])), TDD=as.numeric(x[[i]]))
        }, x=e, years=yrs, eco=eco_IDs))
    d
}

d <- lapply(list(tdd05, tdd10, tdd15, tdd20), make_TDD_df, extractBy=eco_shp, projectTo=sos, resampleTo=ecomask)
lapply(1:length(d), function(i, x, pct) x[[i]][, Threshold := pct[i]], x=d, pct=paste0(c("05",10,15,20), "pct"))
d <- rbindlist(d)
setcolorder(d, c("Region", "Year", "Threshold", "TDD"))
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
save(d, d.stats, sos, ecomask, yrs, cbpal, file="data.RData")

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

# ggplot setup
g <- ggplot(data=d, aes(x=TDD)) + theme(legend.position="bottom")

# @knitr plot_eco
levelplot(ecomask, main="AK level 2 ecoregions", par.settings=revRasterTheme, scales=list(draw=FALSE), contour=F, margin=F, at=at.vals, colorkey=colkey)

# @knitr plot_marginal_tdd_01a
# tdd marginal distributions (across years) by threshold and ecoregion
(p01a <- g + geom_line(aes(colour=Threshold), stat="density", size=1) + facet_wrap(~ Region))

# @knitr plot_marginal_tdd_01b
# tdd distributions by ecoregion and year | threshold = 10%
dsub <- d[Threshold=="10pct",]
(p02a <- g + geom_line(data=dsub, aes(x=TDD, group=Year), stat="density", alpha=0.5) + facet_wrap(~ Region))

# @knitr tables_tdd_mean_sd
setorder(d.stats, Threshold, Mean)
subset(d.stats, Threshold=="05pct", -1) %>% knitr::kable(digits=0, caption="5 % Thaw Degree Days.")
subset(d.stats, Threshold=="10pct", -1) %>% knitr::kable(digits=0, caption="10 % Thaw Degree Days.")
subset(d.stats, Threshold=="15pct", -1) %>% knitr::kable(digits=0, caption="15 % Thaw Degree Days.")
subset(d.stats, Threshold=="20pct", -1) %>% knitr::kable(digits=0, caption="20 % Thaw Degree Days.")
#d.stats %>% filter(Threshold=="05pct") %>% arrange(Mean)
#d.stats %>% filter(Threshold=="10pct") %>% arrange(Mean)
#d.stats %>% filter(Threshold=="15pct") %>% arrange(Mean)
#d.stats %>% filter(Threshold=="20pct") %>% arrange(Mean)
