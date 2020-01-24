## Code for: Fugère V, Hébert MP, Costa NB, Xu CCY, Barrett RDH, Beisner BE, Bell G, 
## Fussmann GF, Shapiro BJ, Yargeau V, Gonzalez A. (2020) Community rescue in experimental
## phytoplankton communities facing severe herbicide pollution.

## Code by Vincent Fugère, 2016-2020
## Questions: vincent.fugere@mail.mcgill.ca

rm(list=ls())

#### libraries and functions ####

library(tidyverse)
library(scales)
library(shape)
library(readxl)
library(RColorBrewer)
library(vegan)
library(magrittr)
library(mgcv)
library(itsadug)
library(viridis)
library(codyn)
library(party)
library(plotrix)
library(boot)

sessionInfo()

# R version 3.6.1 (2019-07-05)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS Sierra 10.12.6
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_CA.UTF-8/en_CA.UTF-8/en_CA.UTF-8/C/en_CA.UTF-8/en_CA.UTF-8
# 
# attached base packages:
#   [1] stats4    grid      stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] boot_1.3-23        plotrix_3.7-7      party_1.3-3        strucchange_1.5-2  sandwich_2.5-1     zoo_1.8-6         
# [7] modeltools_0.2-22  mvtnorm_1.0-11     codyn_2.0.3        viridis_0.5.1      viridisLite_0.3.0  itsadug_2.3       
# [13] plotfunctions_1.3  mgcv_1.8-31        nlme_3.1-142       magrittr_1.5       vegan_2.5-6        lattice_0.20-38   
# [19] permute_0.9-5      RColorBrewer_1.1-2 readxl_1.3.1       shape_1.4.4        scales_1.1.0       forcats_0.4.0     
# [25] stringr_1.4.0      dplyr_0.8.3        purrr_0.3.3        readr_1.3.1        tidyr_1.0.0        tibble_2.1.3      
# [31] ggplot2_3.2.1      tidyverse_1.3.0   
# 
# loaded via a namespace (and not attached):
#   [1] httr_1.4.1         jsonlite_1.6       splines_3.6.1      modelr_0.1.5       assertthat_0.2.1   coin_1.3-1        
# [7] cellranger_1.1.0   pillar_1.4.2       backports_1.1.5    glue_1.3.1         rvest_0.3.5        colorspace_1.4-1  
# [13] Matrix_1.2-18      pkgconfig_2.0.3    broom_0.5.2        haven_2.2.0        generics_0.0.2     TH.data_1.0-10    
# [19] withr_2.1.2        lazyeval_0.2.2     cli_2.0.0          survival_3.1-8     crayon_1.3.4       fs_1.3.1          
# [25] fansi_0.4.0        MASS_7.3-51.4      xml2_1.2.2         tools_3.6.1        hms_0.5.2          multcomp_1.4-11   
# [31] matrixStats_0.55.0 lifecycle_0.1.0    munsell_0.5.0      reprex_0.3.0       cluster_2.1.0      compiler_3.6.1    
# [37] rlang_0.4.2        rstudioapi_0.10    codetools_0.2-16   gtable_0.3.0       DBI_1.0.0          R6_2.4.1          
# [43] gridExtra_2.3      lubridate_1.7.4    zeallot_0.1.0      libcoin_1.0-5      stringi_1.4.3      parallel_3.6.1    
# [49] Rcpp_1.0.3         vctrs_0.2.0        dbplyr_1.4.2       tidyselect_0.2.5  

#custom functions
se <- function(x){sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))}
make.italic <- function(x) as.expression(lapply(x, function(y) bquote(italic(.(y)))))
'%!in%' <- function(x,y)!('%in%'(x,y))
log10p <- function(x) log10(1+x)
plot_chull<-function(xcoord, ycoord, lcolor, llty){
  hpts <- chull(x = xcoord, y = ycoord)
  hpts <- c(hpts, hpts[1])
  lines(xcoord[hpts], ycoord[hpts], col = lcolor, lty = llty)
}  
poly <- function(x, upper, lower, fill){
  polygon(x=c(x, rev(x), x[1]), y=c(upper, rev(lower), upper[1]),border=NA,col=fill)
}

#cols <- c(brewer.pal(11, 'RdYlBu')[c(10:7,5:2)],'#000000') #alternative red-blue color palette
cols <- c('#F0D500', '#95E552', '#00E4A5', '#00D2DA','#41B0EF','#A482DD','#BA5AA7','#A3475E','#000000') #viridis palette

#### design data ####

Sampling.dates <-c('17/08/16','23/08/16','31/08/16','15/09/16','20/09/16','23/09/16','26/09/16','28/09/16','30/09/16','04/10/16','12/10/16')
Sampling.dates <- as.Date(Sampling.dates, format = '%d/%m/%y') %>%
  format('%j') %>% as.numeric
Sampling.dates <- Sampling.dates - 229 #so that 1 = beginning of experiment

#glyphosate pulses (2 in Phase 1, 1 in Phase 2)
pulse.dates <- as.Date(c('22/08/2016','19/09/2016','30/09/2016'), format = '%d/%m/%Y') %>% format('%j')
pulse.dates <- as.numeric(pulse.dates) - 229

#glyphosate concentrations. Concentration of active ingredient (glyphosate isopropylamine salt) converted to acid equivalent
gly.conc <- c(0,49.74455746,135.2197266,367.5653257,999.1461456,2715.960812,7382.746921,20068.3868)*0.7365676

#treatment key
treat <- read_excel('Fugere et al 2020 Nature Ecology and Evolution.xlsx', sheet = 'treatments') %>% 
  mutate(nut.f = as.factor(nut.f)) %>% as.data.frame
treat$nut.f <- relevel(treat$nut.f, 'low')

#### Figure S1: glyphosate and nutrient time series ####

#glyphosate data
gly <- read_xlsx('Fugere et al 2020 Nature Ecology and Evolution.xlsx', sheet = 'glyphosate') %>%
  mutate_at(vars(site, time.point), funs(as.factor)) %>% droplevels

#plotting parameters
gly$date <- as.Date(gly$date, format = '%d.%m.%Y') %>% format('%j') %>% as.numeric
gly$date <- gly$date - 229 #day 1 of exp is Julian day 230
gly <- gly[order(gly$date),]
gly$nut.f <- treat$nut.f[match(gly$site,treat$pond)]
gly$gly.trt <- as.factor(treat$gly.lvl[match(gly$site,treat$pond)])
gly$pch <- 16
gly$pch[gly$nut.f == 'high'] <- 15
gly$lty <- 1
gly$col <- as.numeric(gly$gly.trt)
gly$col[gly$site == 'E1' | gly$site == 'H1'] <- 9
gly$date.idx <- gly$date
gly$date.idx[gly$nut.f == 'low'] <- gly$date[gly$nut.f == 'low'] - 0.4
gly$date.idx[gly$nut.f == 'high'] <- gly$date[gly$nut.f == 'high'] + 0.4
gly$gly.measured.ppm <- gly$gly.measured.ppb/1000
gly$log.gly.m <- log10(1+(gly$gly.measured.ppb))

#nutrient data
nutdat <- read_xlsx('Fugere et al 2020 Nature Ecology and Evolution.xlsx', sheet = 'nutrients') 

#time series data frame for plot
nutTS <- filter(nutdat, time.point != 4) %>%
  spread(var, conc) %>% ungroup %>% mutate_at(c('site', 'time.point'), as.factor)
nutTS <- rename(nutTS, 'date' = day)
nutTS$nut.f <- treat$nut.f[match(nutTS$site,treat$pond)]
nutTS$gly.trt <- as.factor(treat$gly.lvl[match(nutTS$site,treat$pond)])
nutTS$pch <- 16
nutTS$pch[nutTS$nut.f == 'high'] <- 15
nutTS$lty <- 1
nutTS$col <- as.numeric(nutTS$gly.trt)
nutTS$col[nutTS$site == 'E1' | nutTS$site == 'H1'] <- 9
nutTS$date.idx <- nutTS$date
nutTS$date.idx[nutTS$nut.f == 'low'] <- nutTS$date[nutTS$nut.f == 'low'] - 0.4
nutTS$date.idx[nutTS$nut.f == 'high'] <- nutTS$date[nutTS$nut.f == 'high'] + 0.4

pdf('FigS1.pdf',width=4.5,height=6,pointsize=8)

par(mfrow=c(3,1),cex=1,mar=c(4,4,1,1))

plot(log.gly.m~date,data = gly,type='n',yaxt='n',xaxt='n',cex.axis=1,ann=F,bty='l',xlim=c(1,52),ylim=range(gly$log.gly.m)*1.1)
title(ylab=log[10](1+glyphosate)~(mu*g~L^-1),line=2.5)
axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
abline(v=pulse.dates,lty=3)
for (i in 1:34){
  tmp <- subset(gly, site == levels(gly$site)[i])
  if(nrow(tmp) == 5){
    points(log.gly.m~date.idx,tmp,type='l',col=alpha(cols[tmp$col],0.8),lwd=0.5,lty=tmp$lty)
  }
  points(log.gly.m~date.idx,tmp,type='p',pch=tmp$pch[1],col=alpha(cols[tmp$col],0.8),cex=1)
}
text(x=pulse.dates[1],y=6000,label=expression(PHASE~I~-~italic(dose~1)),pos=4)
text(x=pulse.dates[2],y=6000,label=expression(italic(dose~2)),pos=4)
text(x=pulse.dates[3],y=6000,'PHASE II',pos=4)

plot(TN~date,data = nutTS,type='n',yaxt='n',xaxt='n',cex.axis=1,ann=F,bty='l',xlim=c(1,52),ylim=range(nutTS$TN)*c(0.9,1.1),log='y')
title(ylab=expression(paste('TN (',mu,g~L^-1,')')),line=2.5)
axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
abline(v=pulse.dates,lty=3)
for (i in 1:34){
  tmp <- subset(nutTS, site == levels(nutTS$site)[i])
  points(TN~date.idx,tmp,type='o',col=alpha(cols[tmp$col],0.8),pch=tmp$pch[1],lwd=0.5,lty=tmp$lty,cex=1)
}

plot(TP~date,data = nutTS,type='n',yaxt='n',xaxt='n',cex.axis=1,ann=F,bty='l',xlim=c(1,52),ylim=range(nutTS$TP)*c(0.9,1.1),log='y')
title(ylab=expression(paste('TP (',mu,g~L^-1,')')),line=2.5)
axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
title(xlab="day", cex.lab=1,line=2.5)
axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
abline(v=pulse.dates,lty=3)
for (i in 1:34){
  tmp <- subset(nutTS, site == levels(nutTS$site)[i])
  points(TP~date.idx,tmp,type='o',col=alpha(cols[tmp$col],0.8),pch=tmp$pch[1],lwd=0.5,lty=tmp$lty,cex=1)
}

dev.off()

rm(nutTS)

#### Figure S2: other environmental parameters ####

pdf('FigS2.pdf',width=6.2,height=5.5,pointsize=8)

layout(cbind(c(1,1,6,2,2,2),c(3,3,4,4,5,5)))
par(cex=1,mar=c(4,4,1,1))

# depth

depth <- read_excel('Fugere et al 2020 Nature Ecology and Evolution.xlsx', sheet = 'depth') %>% mutate(site = as.factor(site))

plot(depth.cm~day,data = depth,type='n',yaxt='n',xaxt='n',cex.axis=1,ann=F,bty='l',ylim=range(depth$depth.cm)*c(0.95,1.04))
title(ylab=depth~(cm),line=2.5)
axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
title(xlab="day", cex.lab=1,line=2.5)
axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
abline(v=pulse.dates,lty=3)
ln.cols <- c(1,1,9,9,1,1)
for (i in 1:6){
  tmp <- subset(depth, site == levels(depth$site)[i])
  points(depth.cm~day,tmp,type='l',col=alpha(cols[ln.cols[i]],0.8))
}
text(x=pulse.dates[1],y=60,label=expression(PHASE~I~-~italic(dose~1)),pos=4)
text(x=pulse.dates[2],y=60,label=expression(italic(dose~2)),pos=4)
text(x=pulse.dates[3],y=60,'PHASE II',pos=4)

## temperature

temp <- read_excel('Fugere et al 2020 Nature Ecology and Evolution.xlsx', sheet = 'temperature') %>% as.data.frame

xmax <- 2750
plot(x=0,type='n',yaxt='n',xaxt='n',xlim=c(1,xmax),ylim=c(7,32),cex.axis=1,ann=F,bty='l')
title(ylab='temperature (ºC)', cex.lab=1)
title(xlab="day", cex.lab=1)
axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
axlab <- c(0,10,20,30,40,50)
axis(1, at = axlab*48, labels = as.character(axlab), cex.axis=1,lwd=0,lwd.ticks=1) #48 measurements per day
abline(v=pulse.dates*48,lty=3)
for(i in 1:n_distinct(temp$mesocosm)){
  site <- unique(temp$mesocosm)[i]
  sample <- filter(temp, mesocosm == site)
  col.idx <- treat[treat$pond  == site, 'gly.lvl']
  if(site %in% c('E1','H1')){col.idx<-9}
  points(y=sample$temp,x=1:nrow(sample),type='l',col=alpha(cols[col.idx],0.3))
}

rm(sample)

## other physico-chem parameters (YSI)

YSI <- read_excel('Fugere et al 2020 Nature Ecology and Evolution.xlsx', sheet = 'YSI') %>% as.data.frame

YSI$nut.f <- treat$nut.f[match(YSI$site,treat$pond)]
YSI$gly.trt <- as.factor(treat$gly.lvl[match(YSI$site,treat$pond)])
YSI$pch <- 16
YSI$pch[YSI$nut.f == 'high'] <- 15
YSI$lty <- 1
YSI$col <- as.numeric(YSI$gly.trt)
YSI$col[YSI$site == 'E1' | YSI$site == 'H1'] <- 9
#adding a small x axis offset to better distinguish high and low nut ponds
YSI$day.idx <- YSI$day
YSI$day.idx[YSI$nut.f == 'low'] <- YSI$day[YSI$nut.f == 'low'] - 0.6
YSI$day.idx[YSI$nut.f == 'high'] <- YSI$day[YSI$nut.f == 'high'] + 0.6

plot(cond~day,data = YSI,type='n',yaxt='n',xaxt='n',cex.axis=1,ann=F,bty='l',xlim=c(1,52),ylim=range(YSI$cond))
title(ylab=expression(paste('SPC (',mu,S~cm^-1,')')),line=2.5)
axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
abline(v=pulse.dates,lty=3)
for (i in 1:34){
  tmp <- subset(YSI, site == unique(YSI$site)[i])
  points(cond~day.idx,tmp,type='o',col=alpha(cols[tmp$col],0.8),pch=tmp$pch[1],lwd=0.5,lty=tmp$lty,cex=1)
}

plot(do~day,data = YSI,type='n',yaxt='n',xaxt='n',cex.axis=1,ann=F,bty='l',xlim=c(1,52),ylim=range(YSI$do))
title(ylab=DO~(mg~L^-1),line=2.5)
axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
abline(v=pulse.dates,lty=3)
for (i in 1:34){
  tmp <- subset(YSI, site == unique(YSI$site)[i])
  points(do~day.idx,tmp,type='o',col=alpha(cols[tmp$col],0.8),pch=tmp$pch[1],lwd=0.5,lty=tmp$lty,cex=1)
}

plot(pH~day,data = YSI,type='n',yaxt='n',xaxt='n',cex.axis=1,ann=F,bty='l',xlim=c(1,52),ylim=range(YSI$pH))
title(ylab='pH',line=2.5)
axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
title(xlab="day", cex.lab=1,line=2.5)
axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
abline(v=pulse.dates,lty=3)
for (i in 1:34){
  tmp <- subset(YSI, site == unique(YSI$site)[i])
  points(pH~day.idx,tmp,type='o',col=alpha(cols[tmp$col],0.8),pch=tmp$pch[1],lwd=0.5,lty=tmp$lty,cex=1)
}

rm(tmp)

dev.off()

#### Figure S4: TP vs SRP ####

nut4 <- filter(nutdat, time.point == 4) %>% spread(var, conc) %>%
  ungroup %>% mutate_at(c('site', 'time.point'), as.factor)

pdf('FigS4.pdf',width=2.5,height=2.2,pointsize=8)
par(mfrow=c(1,1),mar=c(4,4,1,1))
plot(SRP~TP,nut4,pch=c(16,15)[c(rep(1,8),rep(2,8))],col=cols[1:8],xlab=total~phosphorus~(ppb),ylab=soluble~reactive~phosphorus~(ppb),bty='l')
title(main=day~35, cex.main = 1,line=0)
dev.off()

rm(nut4)

#### GAMMs with fluoroprobe data ####

FP <- read_excel('Fugere et al 2020 Nature Ecology and Evolution.xlsx', sheet = 'Fluoroprobe') %>%
  mutate(site = as.factor(site)) %>% as.data.frame

FP$logChla <- log10(FP$total)
FP$nut.f <- treat$nut.f[match(FP$site,treat$pond)]
FP$gly.conc <- log10(1+(treat$gly.conc[match(FP$site,treat$pond)]))
FP$gly.trt <- as.factor(treat$gly.lvl[match(FP$site,treat$pond)])
FP$gly.conc.ppb <- treat$gly.conc[match(FP$site,treat$pond)]

#adding measured glyphosate concentrations
gly.tp1 <- gly %>% filter(time.point == 'T1')
gly.tp4 <- gly %>% filter(time.point == 'T4')
gly.tp5 <- gly %>% filter(time.point == 'T5')
FP$gly.measured.ppb <- 0
FP$gly.measured.ppb[FP$date %in% c(8,15,30)] <- gly.tp1$gly.measured.ppb[match(FP$site[FP$date %in% c(8,15,30)],gly.tp1$site)]
FP$gly.measured.ppb[FP$date %in% c(35,38,41,43)] <- gly.tp4$gly.measured.ppb[match(FP$site[FP$date %in% c(35,38,41,43)],gly.tp4$site)]
FP$gly.measured.ppb[FP$date %in% c(49,57)] <- gly.tp5$gly.measured.ppb[match(FP$site[FP$date %in% c(49,57)],gly.tp5$site)]
FP$log.gly.measured <- log10(1+FP$gly.measured.ppb)
FP$imi <- log10(1+(treat$imi.conc[match(FP$site,treat$pond)]))

#plotting parameters for Figure 2
FP$pch <- 16
FP$pch[FP$nut.f == 'high'] <- 15
FP$lty <- 1
FP$col <- as.numeric(FP$gly.trt)
FP$col[FP$site == 'E1' | FP$site == 'H1'] <- 9
FP$date.idx <- FP$date
FP$date.idx[FP$nut.f == 'low'] <- FP$date[FP$nut.f == 'low'] - 0.4
FP$date.idx[FP$nut.f == 'high'] <- FP$date[FP$nut.f == 'high'] + 0.4

#phase 0 model: before treatments

tp1dat <- FP[FP$date < pulse.dates[1],]

gam0 <- gam(logChla ~ nut.f + s(gly.conc, k = 8), data=tp1dat)
summary(gam0)

#imidacloprid
gam0.i <- gam(logChla ~ nut.f + s(gly.conc, k = 6) + s(imi, k = 6), data=tp1dat)
summary(gam0.i)

rm(tp1dat)

#phase 1 model: during treatment period, including random effects and temporal dimension

P1dat <- FP[FP$date < pulse.dates[3],]
P1dat <- P1dat[P1dat$date > pulse.dates[1],]
P1dat <- droplevels(P1dat)
P1dat$o.nut <- as.ordered(P1dat$nut.f)

gamp1 <- gam(logChla ~ nut.f + s(date,k=6) + ti(date,log.gly.measured, k=6) + ti(date,log.gly.measured, by = o.nut, k=6) + s(date, site, bs='fs',k=5, m=2), data=P1dat, method = 'REML')
summary(gamp1)

#imidacloprid
gamp1.i <- gam(logChla ~ nut.f + s(date,k=5) + ti(date,imi, k=5) + ti(date,log.gly.measured, k=5) + ti(date,log.gly.measured, by = o.nut, k=5) + s(date, site, bs='fs',k=5, m=1), data=P1dat, method = 'REML')
summary(gamp1.i)

#phase 2 model

p3dat <- filter(FP, date > 55) %>% droplevels
p3dat$bstart<- FP$total[FP$date == 43]
p3dat$gly.end <- log10(1+(FP$gly.measured.ppb[FP$date == 43]))

p3dat.mod <- p3dat %>% filter(site %!in% c('E1','H1'))
p3gam <- gam(logChla ~ nut.f + s(log10(bstart)) + s(gly.end), data=p3dat.mod)
summary(p3gam)

#imidacloprid
p3gam.i <- gam(logChla ~ nut.f + s(log10(bstart), k = 8) + s(gly.end, k = 8) + s(imi, k = 8), data=p3dat.mod) #forcing k = 8 to fit
summary(p3gam.i)

#### Figure S3 ####

P1tree <- P1dat %>% select(logChla, date, nut.f, gly.measured.ppb, imi) %>%
  mutate(gly.ppm = gly.measured.ppb/1000) %>% select(-gly.measured.ppb) %>%
  rename('day' = date, 'nutrient' = nut.f, 'glyphosate (mg/L)' = gly.ppm)
fit <- ctree(logChla ~ ., data = P1tree, controls = ctree_control(minsplit = 1, testtype = 'MonteCarlo', maxdepth = 3))

pdf('FigS3.pdf',width=6, height = 4, pointsize = 8)
plot(fit, inner_panel=node_inner(fit,pval = T), terminal_panel=node_boxplot(fit, width=0.4,fill='white',ylines=3,id=F))
dev.off()

#### Figure 2 ####

newD <- expand.grid(nut.f = c('low','high'), o.nut = c('low','high'), date = unique(P1dat$date), log.gly.measured = seq(from=min(P1dat$log.gly.measured),to=max(P1dat$log.gly.measured),length.out = 50), site = 'C1') %>% as.list
fit <- predict(gamp1,newD,type='response',se.fit = T,exclude='s(date,site)')
newD <- cbind(newD,as.data.frame(fit))
newD$lwr <- newD$fit - 1.96*newD$se.fit
newD$upr <- newD$fit + 1.96*newD$se.fit

pdf('Fig2.pdf',width=4.6,height=4.5,pointsize = 6)
layout(rbind(c(1,1,1,1),c(2,3,4,5),c(6,6,7,7)),heights = c(1,0.7,0.8))

par(mar = c(4,4,1,1), cex = 1)

plot(total~date,data = FP,type='n',yaxt='n',xaxt='n',cex.axis=1,ann=F,bty='l',log='y',ylim=range(FP$total)*1.1)
title(ylab=chl.~italic(a)~(mu*g~L^-1),line=2.5)
axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=c(0.1,1,10,100,1000),labels = as.character(c(0.1,1,10,100,1000)))
title(xlab="day", cex.lab=1,line=2.5)
axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
for (i in 1:34){
  tmp <- subset(FP, site == levels(FP$site)[i])
  points(total~date.idx,tmp,type='l',col=alpha(cols[tmp$col],0.8),lwd=0.5,lty=tmp$lty)
  points(total~date.idx,tmp,type='p',pch=tmp$pch[1],col=alpha(cols[tmp$col],0.6),cex=1.2)
}
abline(v=pulse.dates[2],lty=3,lwd=1)
abline(v=pulse.dates[c(1,3)],lty=1,lwd=2)
text(x=pulse.dates[1],y=range(FP$total)[2],label=expression(PHASE~I~-~italic(dose~1)),pos=4)
text(x=pulse.dates[2],y=range(FP$total)[2],label=expression(italic(dose~2)),pos=4)
text(x=pulse.dates[3],y=range(FP$total)[2],'PHASE II',pos=4)

rect(xleft=8:15,xright=9:16,ybottom = rep(130,8),ytop = rep(160,8),col=cols[1:8],border=NULL,lwd=0.2)
text(x=12,y=122,labels = expression(glyphosate~(mg~L^-1)),pos=3,cex=1)
segments(x0=9:15,x1=9:15,y0=rep(130,7),y1=c(110,120,120,110,120,120,110),lwd=c(0.5,0.3,0.3,0.5,0.3,0.3,0.5))
text(x=c(9,12,15),y=rep(120,3),labels = c('0.04','0.74','14.8'), pos=1,cex=0.84)

points(x=c(8,8),y=c(min(FP$total)+0.05,min(FP$total)*1.7+0.1),pch=c(21,22),col=1)
text(x=c(8,8),y=c(min(FP$total)+0.05,min(FP$total)*1.7+0.1),labels = c('low nutrient','high nutrient'),pos=4)

points(x=c(19.5,20.5),y=rep(min(FP$total)+0.05,2),pch=c(16,15),col=cols[9])
text(x=20,y=min(FP$total)*1.7+0.1,labels = 'controls')

par(mar = c(4,4,1,0), cex = 1)

plot_smooth(gam0, view="gly.conc", plot_all="nut.f",cex.main=0.8,rug=F,col=c('black','deeppink1'),print.summary = F,hide.label = T,yaxt='n',xaxt='n',cex.axis=1,ann=F,bty='l',legend_plot_all = F, h0=NA)
title(main=day~2, cex.main = 1,line=0)
title(ylab=log[10]~chl.~italic(a)~(mu*g~L^-1), cex.lab=1,line=2.5)
axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
title(xlab='future gly. dose',cex.lab=1,line = 2.5)
axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
legend('bottomleft',legend=rev(c('low nutrient','high nutrient')), text.col=rev(c('black','deeppink1')),bty='n',inset=c(-0.15,-0.04))

sub<-newD %>% filter(date == 8)
sub.l <- filter(sub, nut.f == 'low', o.nut == 'low')
sub.h <- filter(sub, nut.f == 'high', o.nut == 'high')
sub <- bind_rows(sub.l,sub.h)
emptyPlot(xlim=range(sub$log.gly.measured),ylim=range(c(sub$lwr,sub$upr)),bty='l')
title(main=day~8, cex.main = 1,line=0)
title(xlab=log[10]~1+gly.~(mu*g~L^-1),cex.lab=1,line = 2.5)
poly(sub.l$log.gly.measured,sub.l$upr,sub.l$lwr,alpha(1,0.25))
lines(fit~log.gly.measured,sub.l,lwd=1,col=alpha(1,0.8))
poly(sub.h$log.gly.measured,sub.h$upr,sub.h$lwr,alpha('deeppink1',0.25))
lines(fit~log.gly.measured,sub.h,lwd=1,col=alpha('deeppink1',0.8))

sub<-newD %>% filter(date == 35)
sub.l <- filter(sub, nut.f == 'low', o.nut == 'low')
sub.h <- filter(sub, nut.f == 'high', o.nut == 'high')
sub <- bind_rows(sub.l,sub.h)
emptyPlot(xlim=range(sub$log.gly.measured),ylim=range(c(sub$lwr,sub$upr)),bty='l')
title(main=day~30, cex.main = 1,line=0)
title(xlab=log[10]~1+gly.~(mu*g~L^-1),cex.lab=1,line = 2.5)
poly(sub.l$log.gly.measured,sub.l$upr,sub.l$lwr,alpha(1,0.25))
lines(fit~log.gly.measured,sub.l,lwd=1,col=alpha(1,0.8))
poly(sub.h$log.gly.measured,sub.h$upr,sub.h$lwr,alpha('deeppink1',0.25))
lines(fit~log.gly.measured,sub.h,lwd=1,col=alpha('deeppink1',0.8))

sub<-newD %>% filter(date == 43)
sub.l <- filter(sub, nut.f == 'low', o.nut == 'low')
sub.h <- filter(sub, nut.f == 'high', o.nut == 'high')
sub <- bind_rows(sub.l,sub.h)
emptyPlot(xlim=range(sub$log.gly.measured),ylim=range(c(sub$lwr,sub$upr)),bty='l')
title(main=day~43, cex.main = 1,line=0)
title(xlab=log[10]~1+gly.~(mu*g~L^-1),cex.lab=1,line = 2.5)
poly(sub.l$log.gly.measured,sub.l$upr,sub.l$lwr,alpha(1,0.25))
lines(fit~log.gly.measured,sub.l,lwd=1,col=alpha(1,0.8))
poly(sub.h$log.gly.measured,sub.h$upr,sub.h$lwr,alpha('deeppink1',0.25))
lines(fit~log.gly.measured,sub.h,lwd=1,col=alpha('deeppink1',0.8))

par(mar = c(4,4,1,1), cex = 1)

#plotting parameters
p3dat$gly.end.plot <- p3dat$gly.end
p3dat$gly.end.plot[p3dat$gly.end.plot > 0] <- p3dat$gly.end.plot[p3dat$gly.end.plot > 0] - 1
gly.tick.loc <- c(0,(log10(1+c(100,1000,10000,50000)) - 1))
gly.labs <- c('0','0.1','1','10','50')
p3dat$gly.conc[c(1,9,19,27)] <- 3
p3dat$gly.conc[17:18] <- 11
p3dat$col <- alpha(cols[p3dat$gly.trt],0.5)
p3dat$col[17:18] <- 1

plot(total~bstart,data = p3dat,type='n',yaxt='n',xaxt='n',cex.axis=1,ann=F,bty='l',xlim=c(1,1000),ylim=c(0.1,150),log='xy')
title(ylab=chl.~italic(a)~after~phase~2~(mu*g~L^-1), cex.lab=1,line=2.5)
axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=c(0.1,1,10,100),labels=c(0.1,1,10,100))
title(xlab=chl.~italic(a)~before~phase~2~(mu*g~L^-1),cex.lab=1)
axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=c(1,10,100,1000))
points(total~bstart,data = p3dat,pch=37-p3dat$pch,bg=col,col=alpha(1,0.5),cex=1.5)

plot(total~gly.end.plot,data = p3dat,type='n',yaxt='n',xaxt='n',cex.axis=1,ann=F,bty='l',ylim=c(0.1,150),xlim=c(-0.1,3.8),log='y')
title(ylab=chl.~italic(a)~after~phase~2~(mu*g~L^-1), cex.lab=1,line=2.5)
axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=c(0.1,1,10,100),labels=c(0.1,1,10,100))
title(xlab=max.~phase~1~glyphosate~(mg~L^-1),cex.lab=1)
axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=gly.tick.loc,labels = gly.labs)
axis.break(axis=1,breakpos=0.4,bgcol="white",breakcol="black",style="slash",brw=0.02)
points(total~gly.end.plot,data = p3dat,pch=37-p3dat$pch,bg=col,col=alpha(1,0.5),cex=1.5)

par(mfrow=c(1,1))

dev.off()

#### Figure S4 ####

pdf('FigS4.pdf',width=3,height=3,pointsize=8)
par(mfrow=c(1,1))
gam.gly <- gam(logChla ~ s(log.gly.measured), data=p3dat.mod)
plot_smooth(gam.gly, view="log.gly.measured", xlim = range(p3dat.mod$log.gly.measured), ylim = range(p3dat.mod$logChla), cex.main=0.8,rug=F,print.summary = F,hide.label = T,yaxt='n',xaxt='n',cex.axis=1,ann=F,bty='l',legend_plot_all = F, h0=NA)
summary(gam.gly)
title(ylab=log[10]~1+chl.~italic(a)~after~phase~2~(mu*g~L^-1), cex.lab=1,line=2.5)
axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
title(xlab=log[10]~1+glyphosate~(mu*g~L^-1),cex.lab=1,line = 2.5)
axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
points(logChla~log.gly.measured,p3dat.mod,pch=p3dat.mod$pch,col=alpha(cols[p3dat.mod$col],0.8))
dev.off()

#### Community composition analysis ####

com <- read_excel('Fugere et al 2020 Nature Ecology and Evolution.xlsx', sheet = 'Biovolume') %>% 
  mutate_at(vars(site,nut.f), as.factor) %>% as.data.frame
com$nut.f <- relevel(com$nut.f, 'low')

#indices
firstcol <- which(colnames(com) == 'Ankistrodesmus')
lastcol <- which(colnames(com) == 'Desmodesmus')

#adding a bunch of community variables
com$abund <- rowSums(com[,firstcol:lastcol]) #total biovolume
com$log.abund <- log10(1+com$abund)
com$richness <- specnumber(com[,firstcol:lastcol])
ctemp <- com[,firstcol:lastcol]
com$r.richness <- as.numeric(rarefy(ctemp, round(min(rowSums(ctemp)),0)))
rm(ctemp)
com$sw <- exp(diversity(com[,firstcol:lastcol], index='shannon')) #effective number of genera

#adding plotting parameters
com$pch <- 16
com$pch[com$nut.f == 'high'] <- 15
com$lty <- 1
com$col <- as.numeric(com$gly.trt)
com$col[com$site == 'E1' | com$site == 'H1'] <- 9
com$date.idx <- com$date
com$date.idx[com$nut.f == 'low'] <- com$date[com$nut.f == 'low'] - 0.5
com$date.idx[com$nut.f == 'high'] <- com$date[com$nut.f == 'high'] + 0.5

#long format dataframe
com.long <- com %>% gather(key = taxon, value = Abundance, Ankistrodesmus:Desmodesmus)

#relative.abundance dataframes
com.rel <- com
com.rel[,firstcol:lastcol] <- com.rel[,firstcol:lastcol] / com.rel$abund
com.rel.long <- com.rel %>% gather(key = taxon, value = Abundance, Ankistrodesmus:Desmodesmus)

#prepping data frames for end of phase 1 and end of phase 2, using relative abundances
phase1end <- subset(com, date == 43)
phase1end$gly.end.plot <- log10(1+phase1end$gly.measured.ppb)
phase1end$gly.end.plot[phase1end$gly.end.plot > 0] <- phase1end$gly.end.plot[phase1end$gly.end.plot > 0] - 1.6
phase1end$gly.measured.ppm <- phase1end$gly.measured.ppb/1000
phase1end$p2chla <- p3dat$total[p3dat$gly.trt %in% c(1,4,7,8)]
endp1chla <- P1dat %>% filter(date == 43) %>% select(gly.trt,total)
phase1end$p1chla <- endp1chla$total[endp1chla$gly.trt %in% c(1,4,7,8)]
phase1end$prctloss <- (1-(phase1end$p2chla/phase1end$p1chla))*100

phase2end <- subset(com.rel, date == 57)
phase2end$gly.end.plot <- log10(1+phase2end$gly.measured.ppb)
phase2end$gly.end.plot[phase2end$gly.end.plot > 0] <- phase2end$gly.end.plot[phase2end$gly.end.plot > 0] - 3.6
phase2end$gly.measured.ppm <- phase2end$gly.measured.ppb/1000
phase2end$p2chla <- p3dat$total[p3dat$gly.trt %in% c(1,4,7,8)]
phase2end$p1chla <- endp1chla$total[endp1chla$gly.trt %in% c(1,4,7,8)]
phase2end$prctloss <- (1-(phase2end$p2chla/phase2end$p1chla))*100
phase2end$p1gly <- phase1end$log.gly.measured

#### Diversity gamms and Figure 3 ####

div.gam.df <- select(phase1end, site, nut.f, log.gly.measured, imi, richness, sw)
div.gam.df$o.nut <- as.ordered(div.gam.df$nut.f)

gam.rich <- gam(richness ~ nut.f + s(log.gly.measured, k = 6) + s(log.gly.measured, k = 6, by = o.nut), data=div.gam.df, family='poisson')
summary(gam.rich)

#refit for plotting purposes only: get the same thing as when using predict as above (Fig 2), but this is much simpler
gam.rich <- gam(richness ~ nut.f + s(log.gly.measured, k = 6, by = nut.f), data=div.gam.df, family='poisson')

gam.alpha <- gam(sw ~ nut.f + s(log.gly.measured, k = 6) + s(log.gly.measured, k = 6, by = o.nut), data=div.gam.df)
summary(gam.alpha)

gam.alpha <- gam(sw ~ nut.f + s(log.gly.measured, k = 6,  by = nut.f), data=div.gam.df)

pdf('Fig3.pdf',pointsize=6,width=6.25,height=3)

layout(rbind(c(1,2,3),c(4,5,6)),widths = c(0.9,0.35,0.45))
par(mar = c(4,4,1,1), cex = 1)

#a) rarefied richness

plot(richness~date,data = com,type='n',yaxt='n',xaxt='n',cex.axis=1,ann=F,bty='l',ylim=c(2,15))
title(ylab='genus number',line=2.5)
axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=c(2,6,10,14),labels=c('2','6','10','14'))
title(xlab="day", cex.lab=1,line=2.5)
axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
for (i in 1:nlevels(com$site)){
  tmp <- subset(com, site == levels(com$site)[i])
  points(richness~date.idx,tmp,type='l',col=alpha(cols[tmp$col],0.8),lwd=0.7,lty=tmp$lty)
  points(richness~date.idx,tmp,type='p',pch=tmp$pch[1],col=alpha(cols[tmp$col],0.8),cex=1.2)
}
abline(v=pulse.dates[2],lty=3,lwd=1)
abline(v=pulse.dates[c(1,3)],lty=1,lwd=2)
text(x=pulse.dates[1],y=range(com$richness)[2],label=expression(PHASE~I~-~italic(dose~1)),pos=4)
text(x=pulse.dates[2],y=range(com$richness)[2],label=expression(italic(dose~2)),pos=4)
text(x=pulse.dates[3],y=range(com$richness)[2],'PHASE II',pos=4)

plot_smooth(gam.rich, view="log.gly.measured", plot_all="nut.f",cex.main=0.8,rug=F,col=c('black','deeppink1'),print.summary = F,hide.label = T,yaxt='n',xaxt='n',cex.axis=1,ann=F,bty='l',legend_plot_all = F, h0=NA, transform = exp)
title(main=end~of~phase~1, cex.main = 1,line=0)
title(ylab=genus~number, cex.lab=1,line=2.5)
axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
title(xlab=log[10]~1+gly.~(mu*g~L^-1),cex.lab=1,line = 2.5)
axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
text(x=0,y=c(3,4), labels=c('low nutrient','high nutrient'), col=c('black','deeppink1'),adj=0)

p2sub <- select(phase2end, site, nut.f, gly.trt, richness, sw, pch:col)
p2sub$x.fac <- p2sub$col
p2sub$x.fac[p2sub$x.fac == 1] <- 2
p2sub$x.fac[p2sub$x.fac == 4] <- 3
p2sub$x.fac[p2sub$x.fac == 9] <- 1
p2sub$x.fac[p2sub$x.fac == 7] <- 4
p2sub$x.fac[p2sub$x.fac == 8] <- 5

plot(richness~x.fac,data = p2sub,type='n',yaxt='n',xaxt='n',cex.axis=1,ann=F,bty='l',xlim=c(0.5,5.5),ylim=c(2,10))
axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
xlabs <- c('0','0','0.3','5.4','14.8')
title(xlab=Phase~I~gly.~(mg~L^-1),cex.lab=1,line=2.5)
axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=1:5,labels = xlabs)
axis.break(axis=1,breakpos=1.5,bgcol="white",breakcol="black",style="slash",brw=0.04)
points(richness~jitter(x.fac),data = p2sub,pch=37-p2sub$pch,bg=alpha(cols[col],0.8),col=alpha(1,0.8),cex=1.2)
title(main=end~of~phase~2, cex.main = 1,line=0)
means <- aggregate(richness~x.fac, p2sub, FUN = 'mean')
segments(x0=seq(0.75,4.75,1),x1=seq(1.25,5.25,1),y0=means$richness,lwd=2.5)

#b) alpha diversity

plot(sw~date,data = com,type='n',yaxt='n',xaxt='n',cex.axis=1,ann=F,bty='l',ylim=c(0,8))
title(ylab=expression(alpha~diversity),line=2.5)
axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
title(xlab="day", cex.lab=1,line=2.5)
axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
for (i in 1:nlevels(com$site)){
  tmp <- subset(com, site == levels(com$site)[i])
  points(sw~date.idx,tmp,type='l',col=alpha(cols[tmp$col],0.8),lwd=0.7,lty=tmp$lty)
  points(sw~date.idx,tmp,type='p',pch=tmp$pch[1],col=alpha(cols[tmp$col],0.8),cex=1.2)
}
abline(v=pulse.dates[2],lty=3,lwd=1)
abline(v=pulse.dates[c(1,3)],lty=1,lwd=2)

plot_smooth(gam.alpha, view="log.gly.measured", plot_all="nut.f",cex.main=0.8,rug=F,col=c('black','deeppink1'),print.summary = F,hide.label = T,yaxt='n',xaxt='n',cex.axis=1,ann=F,bty='l',legend_plot_all = F, h0=NA)
title(main=end~of~phase~1, cex.main = 1,line=0)
title(ylab=expression(alpha~diversity), cex.lab=1,line=2.5)
axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
title(xlab=log[10]~1+gly.~(mu*g~L^-1),cex.lab=1,line = 2.5)
axis(1,cex.axis=1,lwd=0,lwd.ticks=1)

plot(sw~x.fac,data = p2sub,type='n',yaxt='n',xaxt='n',cex.axis=1,ann=F,bty='l',xlim=c(0.5,5.5),ylim=c(0.5,5.5))
axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
title(xlab=Phase~I~gly.~(mg~L^-1),cex.lab=1,line=2.5)
axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=1:5,labels = xlabs)
axis.break(axis=1,breakpos=1.5,bgcol="white",breakcol="black",style="slash",brw=0.04)
points(sw~jitter(x.fac),data = p2sub,pch=37-p2sub$pch,bg=alpha(cols[col],0.8),col=alpha(1,0.8),cex=1.2)
title(main=end~of~phase~2, cex.main = 1,line=0)
means <- aggregate(sw~x.fac, p2sub, FUN = 'mean')
segments(x0=seq(0.75,4.75,1),x1=seq(1.25,5.25,1),y0=means$sw,lwd=2.5)

dev.off()

#### Shifts in community composition as a function of treatments (GAMs B-C distance, synchrony & pMANOVA) ####

#adding Bray-Curtis dissimilarity

#to com.rel data frame, for plotting purposes only (Fig.4)
bcvec <- numeric(0)
for(i in 1:nlevels(com.rel$site)){
  com.rel %>% filter(site == levels(com.rel$site)[i]) %>%
    select(Ankistrodesmus:Desmodesmus) %>%
    vegdist(method='bray') -> distmat
  bcvec <- c(bcvec,as.matrix(distmat)[,1])
}
com.rel$BCdis <- bcvec

#for com.rel.long data frame, more convenient to use custom function used in Collins (2018) Ecology
#provide exactly the same B-C distances as vegdist()

source('https://portal.edirepository.org/nis/dataviewer?packageid=edi.16.2&entityid=2fe289cee84c532df32ff8b0141f3721')

rc <- com.rel.long %>% filter(date < pulse.dates[3]) %>%
  rate_change_interval_BC(time.var = "date",species.var = "taxon",abundance.var = "Abundance", replicate.var = "site") %>%
  filter(interval == 4) %>% select(-interval,-date) %>%
  mutate(distance = as.numeric(as.character(distance))) %>%
  rename('BCdist' = distance)
phase1end %<>% left_join(rc,by='site')

# Gross index of synchrony
synch <- com.long %>% filter(date < pulse.dates[3]) %>%
  select(site,date,taxon,Abundance) %>%
  synchrony(time.var = "date",abundance.var = "Abundance", replicate.var = "site",species.var = 'taxon',metric='Gross')
phase1end %<>% left_join(synch,by='site')

bc.gam.df <- select(phase1end, site, nut.f, log.gly.measured, imi, BCdist, synchrony)
bc.gam.df$o.nut <- as.ordered(bc.gam.df$nut.f)

gam.bc <- gam(BCdist ~ nut.f + s(log.gly.measured, k = 6) + s(log.gly.measured, k = 6, by = o.nut), data=bc.gam.df, family=betar(link = 'logit'))
summary(gam.bc)

gam.bc <- gam(BCdist ~ nut.f + s(log.gly.measured, k = 6, by = nut.f), data=bc.gam.df, family=betar(link = 'logit'))

gam.bc.i <- gam(BCdist ~ nut.f + s(log.gly.measured, k = 4) + s(imi, k = 4), data=phase1end)
summary(gam.bc.i)

gam.synch <- gam(synchrony ~ nut.f + s(log.gly.measured, k = 6) + s(log.gly.measured, k = 6, by = o.nut), data=bc.gam.df)
summary(gam.synch)

gam.synch <- gam(synchrony ~ nut.f + s(log.gly.measured, k = 6, by = nut.f), data=bc.gam.df)

gam.synch.i <- gam(synchrony ~ s(log.gly.measured, by = nut.f, k = 4) + s(imi, k = 4), data=phase1end)
summary(gam.synch.i)

# NMDS ordination
spe <- com %>% filter(date < 57) %>% select(Ankistrodesmus:Desmodesmus) %>% as.matrix
row.names(spe) <- com %>% filter(date < 57) %>% pull(site)
ordi <- metaMDS(spe, distance = 'bray', k = 2, autotransform = FALSE, trymax = 100)

#permutational MANOVA
spe2 <- com %>% filter(date == 43) %>% select(Ankistrodesmus:Desmodesmus) %>% as.matrix
spe.w <- vegdist(log1p(spe2),method = 'bray')
treat.a <- com %>% filter(date == 43) %>% select(gly.trt,nut.f) %>% mutate_all(as.factor)
adonis(spe.w~treat.a$gly.trt)
adonis(spe.w~treat.a$nut.f)

#### Figure 4 ####

pdf('Fig4.pdf',width=3,height=5,pointsize = 6)

layout(rbind(c(1,1,1),c(2,5,3),c(4,4,4)), heights = c(0.7,0.7,1), widths = c(0.45,0.1,0.45))
par(mar = c(4,4,1,1), cex = 1, oma=c(0,0,0,0))

plot(BCdis~date,data = com.rel,type='n',yaxt='n',xaxt='n',cex.axis=1,ann=F,bty='l',ylim=c(0,1.05))
title(ylab=Bray-Curtis~dissimilarity,line=2.5)
axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=c(0,0.5,1),labels = c('0','0.5','1'))
title(xlab="day", cex.lab=1,line=2.5)
axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
for (i in 1:nlevels(com.rel$site)){
  tmp <- subset(com.rel, site == levels(com.rel$site)[i])
  points(BCdis~date.idx,tmp,type='l',col=alpha(cols[tmp$col],0.8),lwd=0.5,lty=tmp$lty)
  points(BCdis~date.idx,tmp,type='p',pch=tmp$pch[1],col=alpha(cols[tmp$col],0.6),cex=1.2)
}

abline(v=pulse.dates[2],lty=3,lwd=1)
abline(v=pulse.dates[c(1,3)],lty=1,lwd=2)
text(x=pulse.dates[1],y=1.04,label=expression(PHASE~I~-~italic(dose~1)),pos=4)
text(x=pulse.dates[2],y=1.04,label=expression(italic(dose~2)),pos=4)
text(x=pulse.dates[3],y=1.04,'PHASE II',pos=4)

plot_smooth(gam.bc, view="log.gly.measured", plot_all="nut.f",cex.main=0.8,rug=F,col=c('black','deeppink1'),print.summary = F,hide.label = T,yaxt='n',xaxt='n',cex.axis=1,ann=F,bty='l',legend_plot_all = F, h0=NA, transform = inv.logit)
title(main=end~of~phase~1, cex.main = 1,line=0)
title(ylab=Bray-Curtis~dissimilarity, cex.lab=1,line=2.5)
axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
title(xlab=log[10]~1+gly.~(mu*g~L^-1),cex.lab=1,line = 2.5)
axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
text(x=2,y=c(0.6,0.65), labels=c('low nutrient','high nutrient'), col=c('black','deeppink1'),adj=0)

plot_smooth(gam.synch, view="log.gly.measured", plot_all="nut.f",cex.main=0.8,rug=F,col=c('black','deeppink1'),print.summary = F,hide.label = T,yaxt='n',xaxt='n',cex.axis=1,ann=F,bty='l',legend_plot_all = F, h0=NA)
title(main=during~phase~1, cex.main = 1,line=0)
title(ylab=community~synchrony~(italic(eta)), cex.lab=1,line=2.5)
axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
title(xlab=log[10]~1+gly.~(mu*g~L^-1),cex.lab=1,line = 2.5)
axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
abline(h=0,lty=3)

g<-ordi$points[,1:2]
plot(g[,2] ~ g[,1], type = "n",yaxt='n',xaxt='n',ann=F)
title(xlab='NMDS dimension 1',cex.lab=1,line = 2.5)
axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
title(ylab='NMDS dimension 2',cex.lab=1,line = 2.5)
axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
gtp1 <- g[seq(from=1,by=5,length.out = 18),]
plot_chull(gtp1[,1],gtp1[,2],'gray',2)
gtp5 <- g[seq(from=5,by=5,length.out = 18),]
plot_chull(gtp5[,1],gtp5[,2],'gray',3)
rm(gtp1,gtp5)

com2 <- com %>% filter(date < 57)
com2$fill.alpha <- 0.1
com2$fill.alpha[com2$date == 43] <- 0.8
com2$pt.cex <- 1
com2$pt.cex[com2$date > 1 & com2$date < 43] <- 0.3

for(i in 1:nlevels(com2$site)){
  subdat <- g[row.names(g) == levels(com2$site)[i],]
  subdat <- cbind(subdat, com2[com2$site == levels(com2$site)[i],c('pch','col')])
  lncol <- com2$col[com2$site == levels(com2$site)[i]][1]
  points(subdat[1,2]~subdat[1,1], pch = subdat$pch-15, col = alpha(cols[subdat$col],0.8), cex = 1.3)
  points(subdat[5,2]~subdat[5,1], pch = subdat$pch, col = alpha(cols[subdat$col],0.8), cex = 1.3)
}
common.tax <- com.long %>% filter(date %in% c(1,43)) %>% group_by(taxon) %>% summarize(tot = sum(Abundance)) %>% arrange(desc(tot)) %>% as.data.frame
label.subset <- ordi$species[rownames(ordi$species) %in% common.tax[1:10,'taxon'],]
text(label.subset, rownames(label.subset), cex = 1, col = alpha(1,0.5))
legend('topright',bty='n',legend='Stress = 0.14')

dev.off()

#### Zooplankton models and Figure S6 ####

zoo <- read_xlsx('Fugere et al 2020 Nature Ecology and Evolution.xlsx', sheet = 'Zooplankton') %>% filter(phase == 2) #added last bit for revisions
zoo$log.dens <- log10p(zoo$total.ind.per.L)
phase2end$zoo <- zoo$log.dens

#imidacloprid zoo
zoo.sub <- filter(zoo, site %!in% c('E1','H1'))
zoo.sub$imi <- CRpred$imi
gamx <- gam(zoo.sub$log.dens~s(zoo.sub$imi, k=4))
summary(gamx)

pdf('FigS6.pdf',width=3,height=3,pointsize=8)
par(mfrow=c(1,1))
plot(log10(p2chla)~zoo,data = phase2end,type='n',yaxt='n',xaxt='n',cex.axis=1,ann=F,bty='l')
axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
title(xlab=log[10]~1+crustacean~density~(ind.~L^-1), cex.lab=1,line=2.5)
axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
title(ylab=log[10]~chlorophyll~italic(a)~(mu*g~L^-1),cex.lab=1,line = 2.5)
points(log10(p2chla)~zoo,data = phase2end,pch=phase2end$pch,col=alpha(cols[col],0.8),cex=1)
title(main=end~of~phase~2, cex.main = 1,line=0)
dev.off()

#### Predicting final biomass: regression tree, univariate GAMs and Figure 5 ####

#making a data frame that includes all variables used to predict final biomass

phase1end <- cbind(phase1end, g[seq(from = 5, by = 5, length.out = 18),])
CRpred <- phase1end %>% filter(site %!in% c('E1','H1'))

#for revisions, adding zoops end of Phase I & bacterial abundance Phase I and Phase II
zoo1 <- read_xlsx('Fugere et al 2020 Nature Ecology and Evolution.xlsx', sheet = 'Zooplankton') %>% filter(phase == 1)
FC <- read_xlsx('Fugere et al 2020 Nature Ecology and Evolution.xlsx', sheet = 'bacterial abundance')
ba.p1 <- FC %>% filter(phase == 1)
ba.p2 <- FC %>% filter(phase == 2)

#site with biovolume data
com.ponds <- CRpred$site

#target taxa that are abundant at the end of phase II
target.tax <- c('Selenastrum','Ankistrodesmus','Desmodesmus','Chlorella')

CRpred2 <- select(CRpred, site, gly.measured.ppm, target.tax, richness, sw, p2chla, p1chla, MDS1, MDS2)

CRgly <- gly.tp5 %>% filter(site %in% com.ponds) %>% select(site, gly.measured.ppm) %>% rename(final.gly = gly.measured.ppm)
CRpred2 <- left_join(CRpred2, CRgly, by = 'site')

CRzoo <- zoo %>% filter(site %in% com.ponds) %>% select(site, total.ind.per.L) %>% rename(zoo.p2 = total.ind.per.L)
CRpred2 <- left_join(CRpred2, CRzoo, by = 'site')

CRzoo1 <- zoo1 %>% filter(site %in% com.ponds) %>% select(site, total.ind.per.L) %>% rename(zoo.p1 = total.ind.per.L)
CRpred2 <- left_join(CRpred2, CRzoo1, by = 'site')

CRba <- ba.p1 %>% filter(site %in% com.ponds) %>% select(site, BA) %>% rename(ba.p1 = BA)
CRpred2 <- left_join(CRpred2, CRba, by = 'site')

CRba <- ba.p2 %>% filter(site %in% com.ponds) %>% select(site, BA) %>% rename(ba.p2 = BA)
CRpred2 <- left_join(CRpred2, CRba, by = 'site')

#data frame for regression tree
CRtree <- select(CRpred2, gly.measured.ppm:ba.p2)

#data frame for univariate GAMs
CRaic <- CRpred2 %>% select(-site) %>% 
  mutate_at(vars(gly.measured.ppm:Chlorella,), funs(log10p)) %>%
  mutate_at(vars(p2chla, p1chla, final.gly, zoo.p2, zoo.p1, ba.p1, ba.p2), funs(log)) %>%
  mutate_all(funs(scale))

CRtree <- rename(CRtree, 'phase 1 glyphosate (ppm)' = gly.measured.ppm)
fit2 <- ctree(log10(p2chla) ~ ., data = CRtree, controls = ctree_control(minsplit = 1, testtype = 'MonteCarlo'))

#Figure 5a (couldn't merge both panels in one graphics device using plot.ctree function)
pdf('Fig5a.pdf',width=2.5,height = 3, pointsize = 8)
plot(fit2, inner_panel=node_inner(fit2,pval = T),
     terminal_panel=node_boxplot(fit2, width=0.4,fill='white',ylines=3,id=F))
dev.off()

CRaic %<>% select(p2chla, everything())
AICtable <- data.frame('var' = colnames(CRaic[2:ncol(CRaic)]),
                       'AIC' = numeric(ncol(CRaic)-1),
                       'R2' = numeric(ncol(CRaic)-1),
                       'p' = numeric(ncol(CRaic)-1),
                       stringsAsFactors = F)
for(i in 2:ncol(CRaic)){
  y <- as.numeric(CRaic$p2chla) 
  x <- as.numeric(CRaic[,i]) 
  mod <- gam(y~s(x,k=7))
  #plot(mod, xlab=colnames(CRaic)[i])
  AICtable$AIC[i-1] <- AIC(mod)
  AICtable$R2[i-1] <- summary(mod)$r.sq
  AICtable$p[i-1] <- summary(mod)$s.table[4]
}

AICtable$var[c(1,10,8,11,6,7,12,9,13,14,15)] <- c('Phase 1 glyphosate','NMDS dimension 2','Phase 1 biomass','Phase 2 glyphosate','Genus number','Alpha diversity','Phase 2 zooplankton density','NMDS dimension 1','Phase 1 zooplankton density','Phase 1 bacterial abundance','Phase 2 bacterial abundance')
AICtable <- arrange(AICtable, by=AIC)

#Fig 5b (colours added in inkscape during final revisions)
pdf('Fig5b.pdf',width=3.2,height = 3.5, pointsize = 8)
par(mfrow=c(1,1),mar=c(4,12,1,1))
barplot(rev(AICtable$AIC),xlab="AIC (univariate GAM)", horiz=TRUE,names.arg=rev(AICtable$var),las=1,col='light gray')
dev.off()

#### Figure S7 relationship between biovolume and chlorophyll ####

FP2 <- FP %>% select(date:logChla) %>% filter(date %in% c(2,8,15,30,43,57)) #Phase 1 time points with matching BV data
FP2$date[FP2$date < 15] <- FP2$date[FP2$date < 15]-1
bv.chla <- com %>% left_join(FP2)
rm(FP2)
bv.chla$abund <- bv.chla$abund*(10^-9)
bv.chla$logbv <- log10(bv.chla$abund)
bv.chla$logbm <- log10(bv.chla$greens)
bv.chla <- filter(bv.chla, date < 57, site != 'D1') %>% as.data.frame #D1 excluded because  ~ 20% diatoms for 2 weeks
bv.chla <- bv.chla %>% group_by(site) %>% summarise_if(is.numeric, mean)
cor(bv.chla$logbm,bv.chla$logbv)
gam.bv <- gam(logbm ~ s(logbv, k = 7), data=bv.chla)
summary(gam.bv)

pdf('FigS7.pdf',width=3,height = 3,pointsize = 8)
plot_smooth(gam.bv, view="logbv", cex.main=0.8,rug=F,print.summary = F,hide.label = T,yaxt='n',xaxt='n',cex.axis=1,ann=F,bty='l',legend_plot_all = F, h0=NA)
title(xlab=log[10]~biovolume~(mm^3~L^-1), cex.lab=1,line=2.8)
axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
title(ylab=log[10]~chlorophyll~italic(a)~(mu*g~L^-1),cex.lab=1,line = 2.8)
axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
points(logbm ~ logbv, data=bv.chla, pch = 16)
dev.off()

