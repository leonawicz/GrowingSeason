# @knitr setup
setwd("C:/github/GrowingSeason/workspaces")
pkgs <- list("ggplot2", "data.table", "dplyr", "tidyr")
dummy <- capture.output(lapply(pkgs, function(x) library(x, character.only=T)))

load("final_outputs/final_gbm_summary_tables.RData") # ri.out, cv.out, pd.out, d.out, s1, sMean, sMean2002
d.out.eco <- d.out
load("final_outputs/final_gbm_summary_tables_withAK.RData") # ri.out, cv.out, pd.out, d.out, s1, sMean, sMean2002
d.out.ak <- d.out
dir.create(plotDir <- file.path("../plots/final_outputs_regional_vs_ak"), recursive=T, showWarnings=F)

extract_to_dt <- function(x, y, fun, rcp, gcm, run, ...){
  raster::extract(x, y, fun, ...) %>% t %>% data.table %>% setnames(names(y)) %>% mutate(Run=run) %>% melt(id.vars="Run", variable.name="Region", value.name="SOS") %>%
    mutate(Year=as.integer(substr(names(x), 5, 8)), RCP=factor(rcp, levels=c("RCP 6.0", "RCP 8.5")), Model=gcm, Source="Projected") %>% select(RCP, Model, Region, Year, Source, SOS, Run)
}

d.proj.eco <- readRDS("final_outputs/sos_projections.rds")
d.proj.eco.mean <- d.proj.eco %>% group_by(Region, Year, Source) %>% summarise(SOS=mean(SOS), Run=1)
d.eco.smooth <- filter(d.out.eco, is.na(Run) & Source=="Predicted") %>% bind_rows(filter(d.proj.eco, Year > 2010)) %>% group_by(Region, Year) %>% summarise(SOS=mean(SOS), Source="Trend", Run=1)
d.proj.ak <- readRDS("final_outputs/sos_projections_withAK.rds")
d.proj.ak.mean <- d.proj.ak %>% group_by(Region, Year, Source) %>% summarise(SOS=mean(SOS), Run=1)
d.ak.smooth <- filter(d.out.ak, is.na(Run) & Source=="Predicted") %>% bind_rows(filter(d.proj.ak, Year > 2010)) %>% group_by(Region, Year) %>% summarise(SOS=mean(SOS), Source="Trend", Run=1)

deco <- d.proj.eco %>% mutate(GBM="Regional", Year=ifelse(Year<2010, "1958-2009", ifelse(Year>=2055, "2055-2099", "2010-2054")))
dak <- d.proj.ak %>% mutate(GBM="Alaska", Year=ifelse(Year<2010, "1958-2009", ifelse(Year>=2055, "2055-2099", "2010-2054")))

# historical predictions over observations time series
clrs <- c("peru", "royalblue", "black", "red")
png(file.path(plotDir, paste0("gbm_boxplotProjections_byRegion.png")), width=3200, height=1600, res=200)
ggplot(bind_rows(deco, dak), aes(Year, SOS, fill=GBM)) + scale_color_manual(values=clrs) + geom_boxplot(position=position_dodge(width=0.9)) +
  theme_bw() + theme(legend.position="bottom") + ggtitle("Projected start of growing season") +
  facet_wrap(~Region, ncol=3, scales="free")
dev.off()

dx <- filter(d.out, !is.na(Run) & Source!="Bias corrected") %>% mutate(Year = Year - Year %% 10)
dx.proj <- d.proj %>% mutate(Year = Year - Year %% 10)
png(file.path(plotDir, paste0("gbm_TSpreds_byRegion_plusGCMs_boxplots.png")), width=3200, height=1600, res=200)
ggplot(dx, aes(x=Year, y=SOS, colour=Source, group=interaction(Source, Run))) +
  scale_color_manual(values=clrs) +
  geom_line(size=1, alpha=0.25) +
  #geom_line(data=filter(d.out, Source=="Bias corrected" & is.na(Run)), colour="#B8860B", size=1, linetype=1) +
  geom_boxplot(data=dx.proj, colour="#00000030", size=1) +
  geom_line(data=filter(d.out, Source=="Observed" & is.na(Run)), colour="black", size=2) +
  geom_line(data=filter(d.out, Source=="Predicted" & is.na(Run)), colour="black", size=2) +
  geom_line(data=filter(d.out, Source=="Predicted" & is.na(Run)), colour=clrs[2], size=1) +
  geom_line(data=filter(d.out, Source=="Observed" & is.na(Run)), colour=clrs[1], size=1) +
  geom_line(data=filter(d.out, Source=="Predicted" & is.na(Run)), colour=clrs[2], size=1, linetype=2) +
  geom_smooth(data=d.smooth) +
  theme_bw() + theme(legend.position="bottom") + ggtitle("Observed and modeled start of growing season") +
  #scale_x_continuous(breaks=c(1982,1990,2000,2010)) +
  guides(fill=guide_legend(override.aes=list(alpha=1)), colour=guide_legend(override.aes=list(alpha=1))) +
  facet_wrap(~Region, ncol=3, scales="free")
dev.off()
