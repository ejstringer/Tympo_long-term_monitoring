# METADATA ================================================================
## Filename: short-term_surveys.R
## Description:
## 
## R version: 4.3.2 for Windows
## packages: 
 #  unmarked 1.4.3 
 #  flextable 0.9.5 
 #  lubridate 1.9.3 
 #  forcats 1.0.0 
 #  stringr 1.5.1 
 #  dplyr 1.1.4 
 #  purrr 1.0.2 
 #  readr 2.1.5 
 #  tidyr 1.3.1 
 #  tibble 3.2.1 
 #  ggplot2 3.4.4 
 #  tidyverse 2.0.0
## Author: Emily J Stringer, Bernd Gruber, Richard Duncan
#=========================================================================#

# code --------------------------------------------------------------------

# libraries ---------------------------------------------------------------

library(tidyverse)
library(flextable)
library(unmarked)

# data --------------------------------------------------------------------

## load
ged <- read.csv("./data/ged_occupancy_surveys2018.csv",
                stringsAsFactors = F)

str(ged)
## fix factor levels
ged$bare <- paste0(sub('%', '', ged$bare), '%')
ged$bare <- factor(ged$bare, 
                   levels= c("<1%","2-5%","6-15%","16-30%",">30%"))

ged$vegstruct <- factor(ged$vegstruct,
                        levels = c("<3cm","3-5cm","5-10cm","10-15cm",">15cm"))

## summarise plots
gedplot <- ged %>%
  select(objID, grid, site, vegtype, nburrow, ndragons, padragons, bare, vegstruct) %>%
  group_by(objID,site,  grid, vegstruct, bare, vegtype) %>%
  summarise(mburrows = mean(nburrow), 
            padragons = max(padragons),
            ndragons = sum(ndragons))

levels(gedplot$bare)
levels(gedplot$vegstruct)

# survey times -----------

ged %>% 
  group_by(site, Repeat, date) %>%
  summarise(plots_surveyed = n()) %>% 
  ungroup() %>% 
  group_by(site, Repeat) %>% 
  summarise(mindate = min(ymd(date)),
            maxdate = max(ymd(date)),
           # time = maxdate-mindate+1,
            ndays = n(),
            mean_plots = mean(plots_surveyed),
            plots_surveyed = paste(plots_surveyed, collapse = '-'))

# occupancy and detection -------------------------------------------------

md <- ged

# dragon detections
md_trans <- md %>% 
  dplyr::select(objID, grid, Repeat, padragons, vegtype, site, bare, vegstruct) %>% 
  pivot_wider(names_from = Repeat, values_from = padragons, names_prefix = 'rep')

y <- md_trans[,grep('rep', colnames(md_trans))]

# burrow detections
md_burrows <- md %>% 
  dplyr::select(objID, grid, Repeat, nburrow, vegtype, site) %>% 
  pivot_wider(names_from = Repeat, values_from = nburrow, names_prefix = 'rep')

## observation variable
obsCovs <- list(burrow=md_burrows[,grep('rep', colnames(md_burrows))])


# occupancy variables

siteCovs <- md_trans[,c('site', 'grid', "vegtype", 'bare', 'vegstruct')]

siteCovs$mburrows <- rowMeans(as.matrix(obsCovs$burrow)) # mean burrows per plot

## unmarked data format ----------
mdocc <- unmarkedFrameOccu(y = y, siteCovs = siteCovs, obsCovs = obsCovs)
summary(mdocc)

# structure models

mB <- occu(~1 ~bare, mdocc)
mB

mV <- occu(~1 ~vegstruct, mdocc)
mV

# Null model

fm1 <- occu(~1 ~1, mdocc)

backTransform(fm1, 'state') # occupancy
backTransform(fm1, 'det')   # detection

# full model

fm <- occu(~grid + burrow ~site + grid + vegtype +mburrows, mdocc)
fm


# reduced model
fmd <- occu(~1~grid+mburrows, mdocc)
fmd

backTransform(fmd, 'det')   # detection

#
fmlist <- fitList(null = fm1, reduced = fmd, full = fm)
modSel(fmlist)

aictable <- modSel(fmlist)@Full[, c('model', 'n', 'AIC',
                                    'delta', 'AICwt', 
                                    'cumltvWt', 'formula')]

aictable %>% 
  mutate_if(is.numeric, round, 2) %>% 
flextable() %>% 
  width(j = 7, width = 2.1, unit = "in") %>% 
  bold(part = 'header') %>% 
  border_outer() %>% 
  hline(part = 'header', 
        border = fp_border_default(color = 'grey', width = 2)) %>% 
  save_as_docx(path = './figures/occupany_models_aic.docx')


  ## model table --------

sum_fmd <- summary(fmd)

rbind(cbind(Type = c('Occupancy', NA,NA), 
                Parameters = rownames(sum_fmd$state), 
                round(sum_fmd$state, 3)),
      cbind(Type = c('Detection'), 
                Parameters = rownames(sum_fmd$det), 
                round(sum_fmd$det, 3))) %>% 
  flextable() %>% 
  flextable::hline(i = nrow(sum_fmd$state)) %>% 
  save_as_docx(path = './figures/occupany_model_final.docx')

# predictions -------------------------------------------------------------
# https://cran.r-project.org/web/packages/unmarked/vignettes/unmarked.html

chisq <- function(fm) {
  umf <- fm@data
  y <- umf@y
  y[y>1] <- 1
  fv <- fitted(fm)
  sum((y-fv)^2/(fv*(1-fv)), na.rm=TRUE)
}

set.seed(1)
(pb <- parboot(fmd, statistic=chisq, nsim=100, parallel=FALSE))


## predict detecton ---------
backTransform(fmd, 'det')   # detection


## predict occupancy --------
newData <- data.frame(grid = factor(rep(c('outside','inside'),
                                        each = 26)),mburrows = 0:25)
fit <- predict(fmd, type = 'state', newdata = newData, appendData=TRUE) %>% 
  mutate_if(is.numeric, round, 3) 

### plot occupancy -------------
myCol <- viridisLite::viridis(8)[c(6,4)]

p<- ggplot(fit, aes(mburrows, Predicted, colour = grid, fill = grid))+
  geom_jitter(data = gedplot, aes(y = padragons), height = 0.05, 
              size = 2, alpha = 0.7)+
  geom_ribbon(aes(ymin = lower, ymax = upper), colour= NA,alpha = 0.2)+
  geom_line(linewidth = 1)+
  scale_color_manual(values =myCol, name = 'Monitoring grids')+
  scale_fill_manual(values =myCol, name = 'Monitoring grids')+
  theme_classic()+
  theme(legend.position = 'inside',
        legend.position.inside = c(0.15, 0.7),
        legend.background = element_rect(colour = 'grey'))+
  xlab(expression('Number of burrows / 20 m'^2,')'))+ # mean number of burrows/20m^2
  ylab('Occupancy');p

ggsave('./figures/occupancy_burrows.png',
       p, width = 6, height = 4)

# view plot data for plots overlaying monitoring grids
gedplot %>% filter(grid == 'inside')

# vegetation structure and bare ground ------------------------------------
pp <- ggplot(gedplot, aes(x = (bare), y=(vegstruct)))+
  geom_jitter(width = 0.2, height = 0.2, 
              aes(colour = factor(padragons),
                  size = padragons +1),alpha = 0.5)+
  xlab("Bare ground") +
  ylab("Vegetation structure")+
  scale_colour_manual(values = c("#0000FF","#FF0000"),
                      name = NULL,
                      labels = c('Unoccupied', 'Occupied'))+
  scale_size(range = c(2, 4))+
  guides(size = 'none')+
  theme_bw()+
  theme(legend.position = c(0.88,0.9),
        axis.ticks.length = unit(0.3, units = 'cm'),
        legend.background = element_rect(colour = 'black'),
        legend.key.size = unit(0.1, units = 'cm')) ;pp


ggsave('./figures/veg_bare_pa.png',
       pp, width = 6, height = 4)


# randomisation -----------------------------------------------------------
## function --------
rand.fun <- function(mat, pa, nsim) {
  
  # number of plots
  nplot <- nrow(mat)
  # number of plots with dragons
  npa <- sum(pa)
  
  # compute the mean distance between plots with dragons
  act.dist <- mean(dist(mat[pa == 1, ]))
  
  # sample npa plots at random from mat and get the mean distance between plots. Replicate this nsim times
  sim.mat <- replicate(nsim, mean(dist(mat[sample(1:nplot, npa), ])))
  
  # plot the results
  hist(sim.mat, breaks = 20)
  abline(v = act.dist, col = "red", lwd = 2)
  
  # return proportion of randomised values less than the actual value
  prop <- length(sim.mat[sim.mat <= act.dist]) / nsim
  
  ranlist <- list(sim = data.frame(sim = sim.mat, above = sim.mat <= act.dist),
                  actual = act.dist, prob = prop)
  return(ranlist)
}

## run randomisation -----------

# matrix of coordinates
matc <- cbind(as.numeric(gedplot$bare), as.numeric(gedplot$vegstruct))
# plots with dragons
pad <- gedplot$padragons

prob_dragon <- rand.fun(matc, pad, 100000)

prob_dragon

### plot ---------------
ppp<- ggplot(prob_dragon$sim, aes(x = sim))+
  geom_histogram(fill = 'grey60', colour = 'black', alpha = 0.3)+
  theme_classic()+
  xlab('Mean distance')+
  ylab('Frequency')+
  scale_y_continuous(breaks = seq(0,12000, 2000))+
  scale_x_continuous(breaks = seq(0,3.5, 0.5))+
  geom_vline(xintercept = prob_dragon$actual,
             colour = 'red', linewidth = 1); ppp

ggsave('./figures/veg_bare_randomisation.png',
       ppp, width = 5, height = 3)

prob_dragon$prob

