# METADATA ================================================================
## Filename: long-term_monitoring
## Description: This file estimates dragon abundance based on the long-term
##              monitoring grid data. 
## 
## R version: 4.3.2 for Windows
## packages: 
 #  jagsUI 1.6.2 
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

# code -------------------------------------------------------------------

# libraries ---------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(jagsUI)

# functions ---------------------------------------------------------------
rm(list = ls())

unlogit <- function(x) exp(x) / (1 + exp(x))

# load ged data ---------------------------------------------------------------

ged <- read.csv("./data/ged_monitoring_grids2002-2023_De-identified.csv") |>
  rename(grid_code = grid_id)|> 
  glimpse() 

table(ged$site, exclude = NULL)



# remove blanks from site names so these match the surveys names (below)
ged$site <- gsub(" ", "", ged$site)

table(ged$site, exclude = NULL)
# rename grid 
# remove juveniles
# create a variable to distinguish each site*grid
# number grids within each site

ged$svl <- as.numeric(ged$svl)
ged$svl[is.na(ged$svl)] <- 0

ged <- ged |>
  rename(grid.id = grid) |>
  group_by(site, grid.id, year, animal.id) %>% 
  mutate(svl.max = max(as.numeric(svl), na.rm = T),
         svl.max = ifelse(svl.max==0, NA, svl.max)) %>%  
  mutate(age.calc = case_when(
    svl.max <38 ~ 'J',
    svl.max <=48 ~ 'SA',
    svl.max > 48 ~ 'A',
    .default = NA
  )) %>%
   ungroup() %>% 
   filter(age.calc != 'J') %>% 
  mutate(grid = paste(site, grid.id)) |>
  group_by(site) |>
  mutate(grid.n = as.numeric(factor(grid.id))) |>
  glimpse()

# number of years and grids
n.year <- length(levels(factor(ged$year)))
n.grid <- length(levels(factor(ged$grid)))
n.year
n.grid

#load survey data -------------------------------------------------------------
# read in list of surveys
surveys <- read.csv("./data/ged_survey_numbers2002-2023_De-identified.csv") |>
  rename(grid.id = grid) |>
  mutate(grid = paste(site, grid.id),
         grid.n = as.numeric(factor(grid.id))) |>
  group_by(site) |>
  mutate(grid.n = as.numeric(factor(grid.id))) |>
  glimpse()

# table of number of checks by year (rows) and grid (columns)
nn <- surveys |>
  ungroup() |>
  select(year, grid, checks) |>
  pivot_wider(names_from = "year", values_from = "checks") |>
  arrange(grid) |>
  glimpse()

n.nights <- t(as.matrix(nn[, -1]))
n.nights

# ind captures/grid/year ----------------------------------------------------------------------
# number of times individuals were captured on each grid in each year
ind <- ged |>
  group_by(site, grid.n, grid, year, animal.id) |>
  summarise(n = n()) |>
  glimpse()

table(ind$n)
# evidence of overdispersion in captures

## minimum number captures -------------------------------------------------------------------
# minimum number of individuals per grid*year
min.ind <- ind |>
  group_by(site, grid.n, grid, year) |>
  summarise(mna = n()) |>
  glimpse()

ggplot(min.ind, aes(y = mna, x = year, colour = site)) +
  geom_point(size = 3) +
  geom_line() +
  facet_grid(site ~ grid.n) +
  theme_bw() +
  theme(legend.position = "none")

## max caught ---------------------------------------------------------------------
# maximum number caught on any grid
max.n <- max(min.ind$mna)
max.n

# set up array Y with 
# rows i = individuals 
# columns j = year 
# k = grids

# Set up for model -----------------
# include pseudo-individuals up to 4 * max.n
Y <- array(0, dim = c(max.n * 4, n.year, n.grid))

# array with the number of checks/surveys undertaken for each year and grid
J <- Y

# matrix to record the number of individuals captured in each year on each grid
ind.cap <- matrix(0, nrow = n.year, ncol = n.grid)

# make a numeric variable for each grid and year
ind$n.grid <- as.numeric(factor(ind$grid))
ind$n.year <- as.numeric(factor(ind$year))


for (j in 1:n.year) {
  for (k in 1:n.grid) {
    # dragons caught on grid k during year j, with n = number of times caught
    cap <- ind$n[ind$n.year == j & ind$n.grid == k]
    n.cap <- length(cap)
    # write values to the arrays
    if (n.cap > 0) {
      Y[1:n.cap, j, k] <- cap
      ind.cap[j, k] <- n.cap
    }
    # write J from number of checks/surveys calculated earlier
    J[, j, k] <- n.nights[j, k]
  }
}

Y[1:20, , 1]
# check that Y does not exceed J
table(Y > J)

# variation in recapture rates
table(Y)

# number of animals captured each year x grid
ind.cap

# J is number of checks for each animal on each grid in each year
# all years on grid 1
J[1:20, , 1]

## ind animal ####################################################
# row i is individual animal
# column j is year
# sheet k is grid

# variables we need for the model

# Y = number of captures for each individual (i) in each year (j) on each grid (k)
# z = latent variable for whether individuals were present or not
# J = number of checks for each individual (i) in each year (j) on each grid (k)

# latent indicator variable
z <- ifelse(Y > 0, 1, NA)

N_grid <- dim(Y)[3]
N_year <- dim(Y)[2]

# set total number of animals to mna * 4 with a minimum = 20
my <- ifelse(Y > 0, 1, 0)
N_animals <- matrix(nrow = N_year, ncol = N_grid)
for (i in 1:N_grid) {
  N_animals[, i] <- colSums(my[, , i])
}

N_animals <- N_animals * 4
N_animals <- ifelse(N_animals < 20, 20, N_animals)

# set z to 0 for animals outside the bounds above
for (k in 1:N_grid) {
  for (j in 1:N_year) {
    if(N_animals[j, k] < dim(z)[1]) z[((N_animals[j, k] + 1):dim(z)[1]), j, k] <- 0
  }
}

## years ---------------------
# number of years that each grid was sampled
my <- ifelse(n.nights > 0, 1, 0)
ny <- colSums(my)
ny

# years that were sampled for each grid
samp <- matrix(nrow = N_year, ncol = N_grid)
for (k in 1:N_grid) {
  samp[1:ny[k], k] <- which(n.nights[, k] > 0)
}
samp

#--------------------------------- CMR model ----------------------------------
CR_dragon <- 
"model {
  for(k in 1:N_grid){
    for(j in samp[1:ny[k], k]){    # only include those years that were sampled for each grid
      for(i in 1:N_animals[j, k])  {
        Y[i, j, k] ~ dbin(mu[i, j, k], J[i, j, k])
        mu[i, j, k] <- z[i, j, k] * p[i, j, k]
        z[i, j, k] ~ dbern(psi[j, k])
        logit(p[i, j, k]) <- lp[i, j, k]
        # allow for heterogeneity among individuals in capture probability
        lp[i, j, k] ~ dnorm(mu.p, tau.p)
      } #i
    
    # prior for presence probability differs betwen trips and grids
    # use the Link prior
    psi[j, k] ~ dbeta(0.001, 1)
    
    # Derived parameters: the estimated number of individuals on each grid at each trip
    N[j, k] <- sum(z[, j, k])
    
    } #j
  } #k
    
  # priors
  mu.p ~ dnorm(0, 0.001)
  tau.p <- 1 / (sigma.p * sigma.p)
  sigma.p ~ dunif(0, 10)

}"  #model finish

model_temp_file <- tempfile(pattern = 'CR_dragon_', fileext = '_model.txt')
write(CR_dragon, model_temp_file)

mod <- jags(
  model.file = model_temp_file,
  data = list(
    Y = Y,
    J = J,
    z = z,
    N_grid = N_grid,
    N_animals = N_animals,
    samp = samp,
    ny = ny
  ),
  parameters.to.save = c("mu.p", "sigma.p", "N"),
  n.chains = 4,
  n.iter = 10000,
  n.burnin = 4000,
  n.thin = 1,
  parallel = T
)

mod
mod.sum <- mod$summary

mod.sum[1:30, ]

## capture probability -------------------------------------------------------------------------------

unlogit(mod.sum[1, 1])

# plot N
out <-  data.frame(mod.sum[substr(rownames(mod.sum), 1, 2) == "N[", c(1, 2, 3, 7)])
names(out) <- c("N", "se", "lcl", "ucl")
out$name <- rownames(out)
glimpse(out)

# extract years
ny <- unlist(lapply(strsplit(out$name, split = ","), function(x) x[1]))
out$n.year <- as.numeric(substr(ny, 3, nchar(ny)))

# extract grids
ng <- unlist(lapply(strsplit(out$name, split = ","), function(x) x[2]))
out$n.grid <- as.numeric(substr(ng, 1, nchar(ng) - 1))

out$year <- as.numeric(levels(factor(ged$year))[out$n.year])
out$grid <- levels(factor(ged$grid))[out$n.grid]

# merge with site and grid names
sg <- surveys |>
  select(grid, site, grid.id, grid.n) |>
  unique()

out <- left_join(out, sg) |>
  glimpse()

head(out)

ggplot(out, aes(y = N, x = year, colour = site)) +
  geom_line(col = "grey20", lty = 1) +
  geom_errorbar(aes(ymin = lcl, ymax = ucl), colour = "grey20", width = 0) +
  geom_point(aes(fill = site),
             shape = 21,
             col = "black",
             size = 2.5) +
  facet_grid(site ~ grid.n, scales = "free_y") + 
  xlab("Year") +
  ylab("Abundance (N)") +
  theme_bw() +
  theme(legend.position = "none")


## figure abundance --------------------------------------------------------

grid.labs <- c("Grid one", "Grid two","Grid three", "Grid four")
names(grid.labs) <- c(1,2,3,4)

de_identified_sites <- levels(factor(out$site))
site.labs <-  paste('Site', LETTERS[1:length(de_identified_sites)])
names(site.labs) <- de_identified_sites

abundanceFig <- ggplot(out, aes(x = year, y = N, group = site))+
  geom_line()+
  geom_errorbar(aes(ymin = lcl, ymax = ucl), colour = "grey", width = 0) +
  geom_point(shape =16, size = 1.5, colour = 'white')+
  geom_point(shape = 1, size = 1.5, colour = "black")+
  facet_grid(site ~ grid.n, scales = "free",labeller=labeller(grid.n = grid.labs, site = site.labs))+
  scale_x_continuous(breaks = seq(2002, 2022, 2), 
                     labels = c(2002, '', 2006, '',2010, '', 2014, '', 2018, '', 2022 )) +
  labs(x="Year", y="Abundance (N)") +
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.x = element_text(size=8, face = "bold" ),
        strip.text.y = element_text(size=9, face = "bold"),
        strip.background = element_blank(),
        axis.line.y = element_line(colour = 'black', size = 1),
        legend.position = "none");abundanceFig

ggsave('./figures/abundance_estimates.png', abundanceFig, width = 6, height = 7.5)

# correlation ------------------------------------------------------------------
# correlation in abundance in each year among grids
# matrix of abundances for year by grid

ma <- out |>
  select(grid, year, N) |>
  pivot_wider(names_from = grid, values_from = N) |>
  glimpse()

ma <- as.matrix(ma[, -1])
ma[ma<1]<- NA
ma

# number of non NA observations for correlations
no <- crossprod(!is.na(ma))
no

# correlations exluding grids with 3 or fewer years data
cm <- cor(ma, use = "pairwise.complete.obs")

# only pairs with >= x observations
cm[no < 6] <- NA

# remove sites with no obs
cm <- cm[colSums(!is.na(cm)) > 0, colSums(!is.na(cm)) > 0]

# upper triangle
cm[lower.tri(cm, diag = T)] <- NA
cm

par(mfrow = c(3, 1))
hist(cm, breaks = seq(-1, 1, 0.1))
abline(v = 0, col = "red", lwd = 2)

site <- unlist(lapply(strsplit(colnames(cm), split = " "), function(x) x[1]))
site.names <- unique(site)

# within site correlations
ws <- numeric()
ws <- as.vector(cm[site == site.names[1], site == site.names[1]])
for(i in 2:length(site.names)) {
  ws <- c(ws, as.vector(cm[site == site.names[i], site == site.names[i]]))
}

ws <- ws[!is.na(ws)]
hist(ws, breaks = seq(-1, 1, 0.1))
abline(v = 0, col = "red", lwd = 2)

# among site correlations
as <- cm
as[as %in% ws] <- NA
as

hist(as, breaks = seq(-1, 1, 0.1))
abline(v = 0, col = "red", lwd = 2)


## figure correlation ------------------------------------------------------
as <- as.vector(as)
as <- as[!is.na(as)]
cor_data <- data.frame(type = 'Within sites', correlation = ws,
                       mean = mean(ws), se = sd(ws)/sqrt(length(ws)),
                       x = -0.75, y = 2.7, t = qt(0.975, df = length(ws)-1)) %>% 
  rbind(data.frame(type = 'Among sites', correlation = as,
                   mean = mean(as), se = sd(as)/sqrt(length(as)),
                   x = -0.75, y = 2.7, t = qt(0.975, df = length(as)-1))) %>% 
  filter(complete.cases(correlation)) %>% 
  mutate(lwr = mean - se*t,
         upr = mean + se*t,
         est = paste0('mean = ', round(mean,2),
                     ' (CI: ', round(lwr,2), ' - ',round(upr,2), ')'),
         est = ifelse(duplicated(est), NA, est))
cor_data$est[!is.na(cor_data$est)]
table(cor_data$type)


tapply(cor_data$correlation, cor_data$type,
       t.test, alternative = 'greater')


corFig <- ggplot(cor_data, aes(x = correlation))+
  geom_histogram(colour = 'black')+
  facet_wrap(~type, ncol = 1, scale = 'free')+
  geom_vline(xintercept = 0, colour = 'red', size = 1)+
  xlim(-1.1,1.1)+
  ylim(0,3)+
  ylab('Frequency')+
  geom_text(aes(label = est,  x = x, y = y), size = 3)+
  xlab(expression(Correlation~~italic(group('(', r, ')'))))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size=8, face = "bold" ),
        strip.text.y = element_text(size=9, face = "bold"),
        strip.background = element_blank(),
        axis.line.y = element_line(colour = 'black', size = 1),
        legend.position = "none");corFig

ggsave('./figures/correlation_N.png', corFig,width = 6, height = 5)


# save N -------------------------------------------------------------------------
write.csv(out[, c('N', 'se', 'lcl', 'ucl', 'year', 'grid')], "./output/N_estimates.csv", row.names = F)

