# @knitr setup
setwd("C:/github/GrowingSeason/workspaces")
pkgs <- list("rasterVis", "maptools", "ggplot2", "data.table", "dplyr", "tidyr")
dummy <- capture.output(lapply(pkgs, function(x) library(x, character.only=T)))

load("data.RData") # d, d.stats, d.stats2, d.hm, sos, ecomask, yrs, cbpal
load("singlePred_outputs/final_gbm_summary_tables_withAK.RData") # ri.out, cv.out, pd.out, d.out, s1, sMean, sMean2002
shpDir <- "C:/github/DataExtraction/data/shapefiles"
eco_shp <- shapefile(file.path(shpDir, "AK_ecoregions/akecoregions.shp")) %>% spTransform(CRS(projection(sos)))
eco_shp <- unionSpatialPolygons(eco_shp, eco_shp@data$LEVEL_2)
sw_shp <- shapefile(file.path(shpDir, "Political/Alaska.shp")) %>% spTransform(CRS(projection(sos)))
dir.create(plotDir <- file.path("../plots/singlePred_outputs_withAK"), recursive=T, showWarnings=F)

# @knitr not_run

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
dx <- filter(d.out, is.na(Run)) %>% dcast(Region + Year + Run ~ Source, value.var="SOS")
rng <- range(c(dx$Observed, dx$Predicted))
png(file.path(plotDir, paste0("gbm_ObsVsFitted_byRegion.png")), width=2000, height=2000, res=200)
ggplot(dx, aes(x=Predicted, y=Observed, colour=Region)) +
  geom_point(position=position_jitter(), size=2) + geom_smooth(method="lm", se=FALSE) + scale_color_manual(values=c("black", brewer.pal(9, "Set1"))) +
  theme_bw() + theme(legend.position="bottom") + coord_cartesian(xlim=rng, ylim=rng) +
  labs(title="Observed vs. fitted values", x="Predicted", y="Observed") #+
#facet_wrap(~Region, ncol=3)
dev.off()

# partial dependence
pd <- pd.out %>% mutate(Var=as.character(Var), Var2=substr(Var, 1, 16)) %>%
  group_by(Region, Var2, Run) %>% mutate(Ymin=min(Prob)) %>% group_by(Region, Var2)
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

f <- function(x, y, df){
  x0 <- x
  x <- x[!is.na(y)]
  y <- y[!is.na(y)]
  s <- smooth.spline(x, y, df=4)
  approx(s$x, s$y, xout=x0)$y
}

pd.interp <- pd %>% split(paste(.$Region, .$Var)) %>% purrr::map(~pd_merge_data(.x)) %>% bind_rows
pd.ci <- pd.interp %>% group_by(Region, Var, x) %>% summarise(
  n_na=length(which(is.na(y))),
  Prob_mean=ifelse(n_na > 0, NA, mean(Prob, na.rm=TRUE)),
  Prob_025=quantile(Prob, 0.025, na.rm=TRUE),
  Prob_975=quantile(Prob, 0.975, na.rm=TRUE),
  Val=ifelse(n_na > 0, NA, Val),
  y_mean=mean(y, na.rm=TRUE),
  y_min=ifelse(n_na > 0, NA, quantile(y, 0.025, na.rm=TRUE)),
  y_max=ifelse(n_na > 0, NA, quantile(y, 0.975, na.rm=TRUE)))



pd.ci <- pd.ci %>% group_by(Region, Var) %>%
  mutate(x=ifelse(n_na > 0, NA, x), y_mean=ifelse(n_na > 0, NA, y_mean),
         Ymin=min(Prob_mean, na.rm=TRUE), y_mean=f(x,y_mean), y_min=f(x,y_min), y_max=f(x,y_max)) %>% ungroup

# bands
save_pdplot <- function(x, outDir){
  region <- gsub(" ", "", unique(x$Region))
  x <- filter(x, Var!="NANA")
  xlm <- range(x$x, na.rm=TRUE)
  ylm <- range(c(x$y_min, x$y_max), na.rm=TRUE)
  g <- ggplot(x, aes(x=Val)) + facet_wrap(~Var, ncol=2, scales="free") + labs(x="x", y="y")
  g <- g + geom_ribbon(aes(ymin=Ymin, ymax=Prob_mean), fill="orange") +
    scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) +
    geom_ribbon(aes(x=x, ymin=y_min, ymax=y_max), fill="black", alpha=0.2) +
    geom_line(aes(x=x, y=y_mean), size=1) +
    geom_line(aes(x=x, y=y_min)) +
    geom_line(aes(x=x, y=y_max))
  png(file.path(outDir, paste0("gbm_PD_", region, ".png")), width=2000, height=2000, res=200)
  print(g)
  dev.off()
}

pd.ci %>% split(.$Region) %>% purrr::walk(~save_pdplot(.x, outDir=plotDir))

x <- tbl_df(data.frame(pd.ci)) %>%
  select(Region, Var, density_x = Val, density_ymin = Ymin, density_ymax = Prob_mean,
         fitted_x = x, fitted_ymean = y_mean, fitted_ymin = y_min, fitted_ymax = y_max)
saveRDS(x, "singlePred_outputs/pd_data.rds")

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

d.proj <- readRDS("singlePred_outputs/sos_projections_withAK.rds")
d.proj.mean <- d.proj %>% group_by(Region, Year, Source) %>% summarise(SOS=mean(SOS), Run=1)
d.smooth <- filter(d.out, is.na(Run) & Source=="Predicted") %>% bind_rows(filter(d.proj, Year > 2010)) %>%
  group_by(Region, Year) %>% summarise(SOS=mean(SOS), Source="Trend", Run=1)

d.proj.table <- d.proj %>% filter(Year >= 1960) %>% mutate(Decade=Year - Year %% 10) %>% group_by(Region, RCP, Model, Decade) %>%
  summarise(SOS_mean=round(mean(SOS)), SOS_sd=round(sd(SOS), 1))

knitr::kable(d.proj.table, format="latex", digits=1, caption='GCM start of season projections')

d.proj.table2 <- d.proj.table %>% dcast(Region + RCP + Model ~ Decade, value.var="SOS_mean") %>%
  mutate(Historical=`1960`) %>% group_by(Region, RCP, Model) %>% summarise(SOS_delta=`2090`-Historical) %>%
  group_by(Region) %>% summarise(SOS_delta=round(mean(SOS_delta)))

knitr::kable(d.proj.table2, format="latex", digits=0, caption='Start of season change in days between historical and 2090s')

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
  geom_line(data=d.proj, aes(group=interaction(Source, Run, RCP, Model)), colour="#00000010", size=1) +
  geom_smooth(data=d.smooth) +
  theme_bw() + theme(legend.position="bottom") + ggtitle("Observed and modeled start of growing season") +
  #scale_x_continuous(breaks=c(1982,1990,2000,2010)) +
  guides(fill=guide_legend(override.aes=list(alpha=1)), colour=guide_legend(override.aes=list(alpha=1))) +
  facet_wrap(~Region, ncol=3, scales="free")
dev.off()

d.sd <- filter(d.out, Source=="Observed" | (!is.na(Run) & Source!="Bias corrected")) %>% group_by(Region, Source) %>% summarise(SD=sd(SOS))

# predictions and observations historical inter-annual variability
png(file.path(plotDir, paste0("gbm_TSpreds_byRegionSD.png")), width=3200, height=1600, res=200)
ggplot(d.sd, aes(x=Region, y=SD, fill=Source)) +
  scale_fill_manual(values=clrs) + geom_bar(stat="identity", position="dodge") +
  theme_bw() + theme(legend.position="bottom") + ggtitle("Observed and modeled start of growing season historical inter-annual SD")
dev.off()
