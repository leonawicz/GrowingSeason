setwd("C:/github/GrowingSeason/data")
library(rasterVis)

rTheme <- function(region=colorRampPalette(c("darkred", "firebrick1", "white", "royalblue", "darkblue"))(19), ...){
    theme <- custom.theme.2(region=region, ...)
    theme$strip.background$col <- theme$strip.shingle$col <- theme$strip.border$col <- theme$panel.background$col <- theme$axis.line$col <- "transparent"
    theme
}

r.sos <- calc(brick("sos_1982_2010.tif"), mean)
r.template <- crop(raster("../../DataExtraction/data/tas_mean_C_AR5_GFDL-CM3_rcp60_01_2062.tif"), r.sos)
r <- resample(r.sos, r.template)
shp <- shapefile("../../DataExtraction/data/shapefiles/Political/Alaska.shp")

p <- levelplot(r, maxpixels=ncell(r), par.settings=rTheme, xlab=NULL, ylab=NULL, scales=list(draw=F), margin=F, colorkey=F) +
    layer(sp.polygons(shp, col='gray40'), data=list(shp=shp))
png("../plots/_GrowingSeason.png", res=300, width=3200, height=round(p$aspect.ratio*3200))
p
dev.off()
