setwd("C:/github/GrowingSeason/data")

library(rasterVis)
library(grid)
library(gridExtra)
library(ggplot2)
library(data.table)
library(dplyr)
library(gtable)

# Data prep
b <- brick("sos_1982_2010.tif")
r <- subset(b, 1)
b1 <- brick("pct05_tdd_spring_1979_2010.tif") %>% projectRaster(r) %>% resample(r, method="bilinear") %>% mask(r)
b2 <- brick("pct10_tdd_spring_1979_2010.tif") %>% projectRaster(r) %>% resample(r, method="bilinear") %>% mask(r)
b3 <- brick("pct15_tdd_spring_1979_2010.tif") %>% projectRaster(r) %>% resample(r, method="bilinear") %>% mask(r)
b4 <- brick("pct20_tdd_spring_1979_2010.tif") %>% projectRaster(r) %>% resample(r, method="bilinear") %>% mask(r)

ex_data <- function(b, lab, lab.levels=lab){
    d <- rasterToPoints(subset(b,1)) %>% rbind(rasterToPoints(subset(b,29)))
    d <- data.table(d, Season=rep(as.character(c(1982, 2010)), each=nrow(d)/2), `DOY TDD %`=factor(lab, levels=lab.levels))
    setnames(d, c("x", "y", "DOY TDD", "Year", "Percentile"))
    d
}

lev <- c("5%", "10%", "15%", "20%")
d <- ex_data(b1, "5%", lev) %>% rbind(ex_data(b2, "10%", lev)) %>% rbind(ex_data(b3, "15%", lev)) %>% rbind(ex_data(b4, "20%", lev)) %>% data.table

shp <- shapefile("C:/github/DataExtraction/data/shapefiles/Political/Alaska.shp")
shp <- fortify(shp)

# Prep ggplot theme settings
e <- element_blank()
tp_theme <- theme(axis.title.x=e, axis.title.y=e, axis.text.x=e, axis.text.y=e, axis.ticks=e, panel.grid=e, axis.line=e, strip.background=e,
    panel.background=element_rect(fill="transparent", colour=NA), plot.background=element_rect(fill="transparent", colour=NA))

# Store ggplot objects
g1a <- ggplot(d, aes(x, y)) + geom_raster(aes(fill=`DOY TDD`)) +
    geom_polygon(data=shp, aes(long, lat, group=group), color="gray40", fill="transparent") +
    theme(panel.margin=unit(0, "cm"), plot.margin=unit(c(0,0,0,0), "mm"),
        legend.position="bottom", legend.box="horizontal",
        strip.text.x=element_text(size=12), strip.text.y=element_text(size=12, angle=180)) +
    facet_grid(Percentile~Year, switch="both") + tp_theme + scale_fill_distiller(palette = "Spectral") + labs(title="DOY TDD Percentile Distributions")

val.sq <- seq(0, 200, by=50)
ggbp <- geom_boxplot(fill="gray", width=0.5)
ylm <- ylim(range(val.sq) + c(-1, 1)*0.1*diff(range(val.sq)))
ggann <- annotate("text", label=val.sq, x=1.5, y=val.sq, size=3)
g2a <- ggplot(d %>% filter(Percentile=="5%"), aes(x=Year, y=`DOY TDD`)) + ggbp + ylm + ggann + tp_theme
g2b <- ggplot(d %>% filter(Percentile=="10%"), aes(x=Year, y=`DOY TDD`)) + ggbp + ylm + ggann + tp_theme
g2c <- ggplot(d %>% filter(Percentile=="15%"), aes(x=Year, y=`DOY TDD`)) + ggbp + ylm + ggann + tp_theme
g2d <- ggplot(d %>% filter(Percentile=="20%"), aes(x=Year, y=`DOY TDD`)) + ggbp + ylm + ggann + tp_theme

# Make plot grob and extract component grobs
gt <- ggplotGrob(g1a)
dev.off()
panel <- gtable_filter(gt, "panel")
stripleft <- gtable_filter(gt, "strip-left")
stripbottom <- gtable_filter(gt, "strip-bottom")
guidebox <- gtable_filter(gt, "guide-box")
top <- gtable_filter(gt, "title")

#### Make viewports
# boxplots
just <- c("center", "center")
sq <- c(0.25, 0.5, 0.75, 1)-0.1
w <- 0.3
h <- 0.15
x <- 0.5
vp1 <- viewport(width=w, height=h, x=x, y=sq[1], just=just)
vp2 <- viewport(width=w, height=h, x=x, y=sq[2], just=just)
vp3 <- viewport(width=w, height=h, x=x, y=sq[3], just=just)
vp4 <- viewport(width=w, height=h, x=x, y=sq[4], just=just)
# outer box
outerBox <- viewport(width = unit(1, "npc"), height = unit(1, "npc"))
# inner box
innerBox <- viewport(x=unit(0.525, "npc"), y=unit(0.5, "npc"), width=unit(0.95, "npc"), height=unit(0.9, "npc"))
# left panel labels
Vleft <- viewport(x=unit(0.05, "npc"), y=unit(0.5, "npc"), width=unit(0.05, "npc"), height=unit(0.9, "npc"))
# bottom panel labels and color key
Vbottom <- viewport(x=unit(0.5, "npc"), y=unit(0.025, "npc"), width=unit(1, "npc"), height=unit(1, "npc"))
# top title
Vtop <- viewport(x=unit(0.5, "npc"), y=unit(0.975, "npc"), width=unit(1, "npc"), height=unit(1, "npc"))

# draw plot
png("../plots/test/boxplots_on_maps_test.png", height=3600, width=3000, res=300)
grid.newpage()
pushViewport(outerBox)
#grid.rect(gp = gpar(col = "red", fill = NA))
pushViewport(innerBox)
grid.draw(panel)
#grid.rect(gp = gpar(col = "red", fill = NA, lwd = 2))
print(g2a, vp=vp1)
print(g2b, vp=vp2)
print(g2c, vp=vp3)
print(g2d, vp=vp4)
upViewport()
pushViewport(Vbottom)
grid.draw(stripbottom)
grid.draw(guidebox)
upViewport()
pushViewport(Vleft)
grid.draw(stripleft)
upViewport()
pushViewport(Vtop)
grid.draw(top)
popViewport()
popViewport()
dev.off()
