# landsat
# METADATA ===============================================================
## Filename: landsat_data.R
## Description: This file runs the landsat analysis on NDVI, LST, and dragon 
##              captures.
## 
## R version: 4.3.2 for Windows
## packages: 
 #  brms 2.22.0 
 #  Rcpp 1.0.12 
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
#=======================================================================##

# code -------------------------------------------------------------------
# Google Earth Engine code
# https://code.earthengine.google.com/4e81453e81d8a08c65e3e6166d23400f
# libraries ---------------------------------------------------------------

library(tidyverse)
library(flextable)
library(brms)

# loading -------------------
## landsat ------------
# read in the Landsat data for all grids

landsat <- read.csv("./data/ged_grids_satellite2000-2020_De-identified.csv", 
                    stringsAsFactors = FALSE, strip.white = TRUE) |>
  dplyr::select(year, FVC, LST, NDVI, grid_id)|>
  filter(!is.na(NDVI) & !is.na(LST)) |> 
  glimpse()

# scale satellite variables
landsat$scale_NDVI = (landsat$NDVI - mean(landsat$NDVI)) / sd(landsat$NDVI)
landsat$scale_LST = (landsat$LST - mean(landsat$LST)) / sd(landsat$LST)
landsat$scale_FVC = (landsat$FVC - mean(landsat$FVC)) / sd(landsat$FVC)

## ged captures --------------------------------------------------------------------------
gedcap <- read.csv("./data/ged_monitoring_grids2002-2023_De-identified.csv")
  
capgridyear <- plyr::ddply(gedcap, .variables = c("grid_id",'site','grid', "year"),
                           summarise, ncap=length(grid_id))|>
  filter(complete.cases(grid_id)) %>% 
  mutate(grid_join = paste(gsub(' ', '', site), grid, sep = '_')) %>% 
  select(grid_id, year, grid_join, ncap)

## ged surveys -------------
surveys <- read.csv("./data/ged_survey_numbers2002-2023_De-identified.csv") |>
  filter(checks> 0) %>% 
  mutate(grid_join = paste(site, grid, sep = '_'))
surveys$checks %>% table

# join datasets -----------------------------------------------------------
## captures with surveys --------
N_est <- read.csv('./output/N_estimates.csv')  %>% 
  mutate(grid_join = sub(' ', '_', grid)) %>% 
  select(grid_join, year, N, se, lcl, ucl)
N_est %>% head()
capgridyearsurveys <- full_join(capgridyear[-1], surveys) %>% 
  left_join(unique(capgridyear[,c('grid_join', 'grid_id')])) %>%
  left_join(N_est) %>% 
  select(grid_id, year, checks, ncap, N, se, lcl, ucl) %>% 
  mutate(ncap_std = ncap/checks*18)

table(capgridyearsurveys$ncap, useNA = 'always') # NAs are for surveys 
                                                 # with zero captures
## landsat with geds -----------
landsat_cap <- inner_join(capgridyearsurveys, landsat, by = c('year', 'grid_id')) %>% 
  mutate(ncap = ifelse(is.na(ncap),0, ncap),
         ncap_std = ifelse(is.na(ncap_std),0, ceiling(ncap_std)),
         ncap_pa = ifelse(ncap > 0, 1, 0),
         site = substr(grid_id,1,2),
         N_round = round(N)) %>% 
  ungroup() %>% 
  rename(grid = grid_id) %>% 
  glimpse


# plot captures NDVI and LST ----------------------------------------------

ggplot(landsat_cap, aes(y = scale_LST, x = scale_NDVI, color= N_round)) +
  geom_point(aes(size = (N_round/10 + 1)), show.legend = F) +
  theme_bw()+
  scale_color_gradient(low = 'blue', high = 'red')

# overview ------------------------------------------------------------------------
head(landsat_cap)
table(landsat_cap$grid, landsat_cap$year)


table(landsat_cap$grid)
table(landsat_cap$year)
hist(landsat_cap$N_round, breaks = 20)
# LST and NDVI -------------------------------
# relationship between LST and NDVI
mc <- cor.test(landsat_cap$scale_LST, landsat_cap$scale_NDVI)
mc

#m_predictors <- brm(bf(scale_LST ~ scale_NDVI + (1|year) + (1|grid)), data = landsat_cap)
#summary(m_predictors)

# hurdle model ------------------------------------------------------------------------
# hurdle model with random effects for grid and year

m_hurdle <- brm(bf(N_round ~ s(scale_LST) +s(scale_NDVI) + (1|year) + (1|grid),
             hu ~ s(scale_LST) + s(scale_NDVI) + (1|year) + (1|grid)),
          data = landsat_cap, iter = 4000, 
          family = hurdle_negbinomial(),  # needs integer response
          control = list(adapt_delta = 0.999))

summary(m_hurdle)
m_hurdle

R2_mhurdle <- bayes_R2(m_hurdle)

sum_model <- summary(m_hurdle)
mhurdle_tbl <- rbind(sum_model$spec_pars,
                     sum_model$splines, 
                     sum_model$spec_pars, 
                     sum_model$fixed) %>% 
  mutate_all(round, 2) %>% 
  mutate(Parameter = rownames(.),
         Parameter = sub('_1', '', Parameter)) %>% 
  relocate(Parameter)

mhurdle_tbl[grepl('shape', mhurdle_tbl$Parameter),] <- NA
mhurdle_tbl$Parameter[is.na(mhurdle_tbl$Parameter)] <- c('SPLINES', 'FIXED')
mhurdle_tbl$R2 <- NA
mhurdle_tbl$R2[1] <- round(R2_mhurdle[1],2)
flextable(mhurdle_tbl) %>% 
  autofit() %>% 
  bold(i = c(1,6), j = 1) %>% 
  bold(j = 9, part = 'all') %>% 
  vline(j = 8,
        border = fp_border_default(style = 'solid', color = 'grey',
                                   width = 1)) %>% 
  align(align = 'center', j = 9, part = 'all') %>% 
  hline(i = 5, j = 1:8, border = fp_border_default(style = 'solid', color = 'grey',
                                          width = 2)) %>% 
  border_outer() %>%
  colformat_double(,j = 7:8, big.mark = "", digits = 0) %>% 
  save_as_docx(path = './figures/hurdle_model.docx')

## conditional effects ---------------------
nd1 <- data.frame(scale_NDVI = seq(-2, 3, 0.01),
                  scale_LST = -1)

p1 <-fitted(m_hurdle, newdata = nd1, re_formula = NA)#[, 1]


nd2 <- data.frame(scale_NDVI = -1,
                  scale_LST = seq(-2, 2, 0.01))

p2 <- fitted(m_hurdle, newdata = nd2, re_formula = NA)#[, 1]


landsat_cap$spred <- predict(m_hurdle, re_formula = NA)[, 1]


 nd1$NDVI <-  (nd1$scale_NDVI* sd(landsat$NDVI)  + mean(landsat$NDVI)) 
 nd2$NDVI <-  (nd2$scale_NDVI* sd(landsat$NDVI)  + mean(landsat$NDVI)) 
# 
 nd1$LST <-  (nd1$scale_LST* sd(landsat$LST) + mean(landsat$LST)) 
 nd2$LST <-  (nd2$scale_LST * sd(landsat$LST) + mean(landsat$LST))


## plot --------------------------------------------------------------------------

# set opacity
op <- 0.5
# colour gradient
rbPal <- colorRampPalette(c(rgb(0, 0, 1, op), rgb(1, 0, 0, op)), alpha = TRUE)


png("./figures/landsat_dragons.png", width = 8, height = 8, res = 600, units = 'in')
par(mfrow = c(2, 2))

# plot presence (red) and absence (blue) of dragons
# circle size proportional to number of captures
occ <- ifelse(landsat_cap$ncap > 0, 1, 0)
cl <- ifelse(occ == 1, rgb(1, 0, 0, op), rgb(0, 0, 1, op))
cx <- ifelse(occ == 0, 0.7, landsat_cap$N_round/30 + 0.7)
plot(scale_LST ~ scale_NDVI, data = landsat_cap, pch = 19, cex = cx, col = cl,
     ylab = "scaled LST", xlab = "scaled NDVI",
     main = "Occupancy and abundance")
legend(1, 2.2, legend = c("Occupied", "Unoccupied"), fill = c("#FF0000", "#0000FF"))
text(1.5, 2.4, substitute(italic(r)["LST, NDVI"]~" = "~x~", "~italic(P)~"<0.001", list(x = round(mc$estimate, 2))), 
     cex = 1)

cx <- landsat_cap$spred/5 + 0.7
cl <- rbPal(10)[as.numeric(cut(landsat_cap$spred, breaks = 10))]
plot(scale_LST ~ scale_NDVI, data = landsat_cap, pch = 19, cex = cx, col = cl,
     ylab = "scaled LST", xlab = "scaled NDVI",
     main = "Predicted abundance")

plot(p1[,1] ~ nd1$scale_NDVI, pch = 19, ylim = c(0, 26),col = 'white',
     ylab = "Abundance", xlab = "scaled NDVI", main = "NDVI")
polygon(c(nd1$scale_NDVI, rev(nd1$scale_NDVI)),
        c(p1[,3], rev(p1[,4])),
        col = adjustcolor("grey", 0.3),
        border = adjustcolor("grey", 0.3))  
lines(nd1$scale_NDVI, p1[,1], lwd = 3)

plot(p2[,1] ~ nd2$scale_LST, pch = 19, ylim = c(0, 26), col = 'white',
     ylab = "Abundance", xlab = "scaled LST", main = "LST")

polygon(c(nd2$scale_LST, rev(nd2$scale_LST)),
        c(p2[,3], rev(p2[,4])),
        col = adjustcolor("grey", 0.3),
        border = adjustcolor("grey", 0.3))  
lines(nd2$scale_LST, p2[,1], lwd = 3)

dev.off()



## plot v2 -----------------------------------------------------------------

op <- 0.5
# colour gradient
rbPal <- colorRampPalette(c(rgb(0, 0, 1, op), rgb(1, 0, 0, op)), alpha = TRUE)


png("./figures/landsat_dragons_notscaled.png", width = 8, height = 8, res = 600, units = 'in')
par(mfrow = c(2, 2))

# plot presence (red) and absence (blue) of dragons
# circle size proportional to number of captures
occ <- ifelse(landsat_cap$ncap > 0, 1, 0)
cl <- ifelse(occ == 1, rgb(1, 0, 0, op), rgb(0, 0, 1, op))
cx <- ifelse(occ == 0, 0.7, landsat_cap$N_round/30 + 0.7)
plot(scale_LST ~ scale_NDVI, data = landsat_cap, pch = 19, cex = cx, col = cl,
     ylab = "scaled LST", xlab = "scaled NDVI",
     main = "Occupancy and abundance")
legend(1, 2.2, legend = c("Occupied", "Unoccupied"), fill = c("#FF0000", "#0000FF"))
text(1.5, 2.4, substitute(italic(r)["LST, NDVI"]~" = "~x~", "~italic(P)~"<0.001", list(x = round(mc$estimate, 2))), 
     cex = 1)

cx <- landsat_cap$spred/5 + 0.7
cl <- rbPal(10)[as.numeric(cut(landsat_cap$spred, breaks = 10))]
plot(scale_LST ~ scale_NDVI, data = landsat_cap, pch = 19, cex = cx, col = cl,
     ylab = "scaled LST", xlab = "scaled NDVI",
     main = "Predicted abundance")


plot(p1[,1] ~ nd1$NDVI, pch = 19, ylim = c(0, 26),col = 'white',
     ylab = "Abundance", xlab = "NDVI", main = "NDVI")
polygon(c(nd1$NDVI, rev(nd1$NDVI)),
        c(p1[,3], rev(p1[,4])),
        col = adjustcolor("grey", 0.3),
        border = adjustcolor("grey", 0.3))  
lines(nd1$NDVI, p1[,1], lwd = 3)

plot(p2[,1] ~ nd2$LST, pch = 19, ylim = c(0, 26), col = 'white',
     ylab = "Abundance", xlab = "LST", main = "LST")

polygon(c(nd2$LST, rev(nd2$LST)),
        c(p2[,3], rev(p2[,4])),
        col = adjustcolor("grey", 0.3),
        border = adjustcolor("grey", 0.3))  
lines(nd2$LST, p2[,1], lwd = 3)

dev.off()




