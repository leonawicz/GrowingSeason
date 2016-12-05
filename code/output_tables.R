# @knitr setup
setwd("C:/github/GrowingSeason/workspaces")
pkgs <- list("rasterVis", "maptools", "ggplot2", "data.table", "dplyr", "tidyr")
dummy <- capture.output(lapply(pkgs, function(x) library(x, character.only=T)))

load("data.RData") # d, d.stats, d.stats2, d.hm, sos, ecomask, yrs, cbpal
load("final_outputs/final_gbm_summary_tables_withAK.RData") # ri.out, cv.out, pd.out, d.out, s1, sMean, sMean2002
dir.create(plotDir <- file.path("../plots/final_outputs_withAK"), recursive=T, showWarnings=F)

# @knitr tables
library(printr)
lab <- c("Ecoregion", "DOY TDD 5%", "DOY TDD 10%", "DOY TDD 15%", "DOY TDD 20%")
ri.table <- ri.out %>% group_by(Region, Predictor) %>% summarise(Mean=mean(RI), SD=sd(RI)) %>%
  mutate(`Relative Influence`=paste0(round(Mean, 1), " (",round(SD, 1), ")")) %>%
  dcast(Region ~ Predictor, value.var="Relative Influence") %>% setnames(lab)
#knitr::kable(ri.table, format="latex", digits=1, caption='GBM predictor relative influence')

rsq.table <- filter(d.out, is.na(Run) & Source!="Bias corrected" & Region!="Alaska") %>% dcast(Region + Year ~ Source, value.var="SOS") %>%
  group_by(Region) %>% summarise(`R^2`=cor(Observed, Predicted)^2) %>% setnames(c("Ecoregion", "R^2")) %>% mutate(`Area (km^2)`=table(ecomask[]))

gbm.table <- left_join(ri.table, rsq.table)
knitr::kable(gbm.table, format="latex", digits=2, caption='GBM relative influence and R2 by ecoregion')

d.proj <- readRDS("C:/github/GrowingSeason/workspaces/final_outputs/sos_projections_withAK.rds")

d.proj.table1a <- d.proj %>% filter(Year >= 1960) %>% mutate(Decade=Year - Year %% 10) %>% group_by(Region, RCP, Model, Decade) %>%
  summarise(SOS_mean=round(mean(SOS)), SOS_sd=round(sd(SOS), 1))

d.proj.table1b <- d.proj %>% filter(Year >= 1960) %>% mutate(Decade=Year - Year %% 10) %>% group_by(Region, Decade) %>%
  summarise(SOS_mean=round(mean(SOS)), SOS_sd=round(sd(SOS), 1))

d.proj.table1c <- d.proj %>% filter(Year >= 1960) %>% mutate(Decade=Year - Year %% 10) %>% group_by(Region, Decade) %>%
  summarise(SOS_mean=round(mean(SOS)), SOS_sd=round(sd(SOS), 1)) %>% dcast(Decade ~ Region, value.var="SOS_mean")
names(d.proj.table1c) <- c("Decade", "AK Range", "Aleut Mdws", "Arc Tun", "Bering Tai", "Bering Tun", "Coast Mt", "Coast Rain", "Boreal", "Pacific Mtn")

knitr::kable(d.proj.table1c, format="latex", digits=1, caption='GCM start of season projections')

d.proj.table2a <- d.proj.table1a %>% dcast(Region + RCP + Model ~ Decade, value.var="SOS_mean") %>%
  mutate(Historical=`1960`) %>% group_by(Region, RCP, Model) %>% summarise(SOS_delta=`2090`-Historical) %>%
  group_by(Region) %>% summarise(SOS_delta=round(mean(SOS_delta))) %>% arrange(desc(SOS_delta))

d.proj.table2b <- d.proj.table1b %>% dcast(Region ~ Decade, value.var="SOS_mean") %>%
  mutate(Historical=`1960`) %>% group_by(Region) %>% summarise(SOS_delta=`2090`-Historical) %>%
  group_by(Region) %>% summarise(SOS_delta=round(mean(SOS_delta))) %>% arrange(desc(SOS_delta))

knitr::kable(d.proj.table2b, format="latex", digits=0, caption='Start of season change in days between historical and 2090s')
