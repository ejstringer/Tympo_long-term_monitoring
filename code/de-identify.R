library(tidyverse)

# load --------------------------------------------------------------------


deidentify <- read.csv("./data/ged_monitoring_grids2002-2023.csv")
head(deidentify)

sat <- read.csv('./data/ged_grids_satellite2000-2020.csv')
head(sat)

surveys <- read.csv('./data/ged_survey_numbers2002-2023.csv')
head(surveys)



# codes -------------------------------------------------------------------

sitenames <- levels(factor(deidentify$site))
sitecoded <- paste('Site', LETTERS[1:length(sitenames)])

animalnames <- levels(factor(deidentify$animal.id))
animalcoded <- paste('id', 1:length(animalnames))

gridnames <- levels(factor(c(deidentify$grid_id, sat$grid_id)))
gridcoded <- paste('grid', 1:length(gridnames), sep = '-')



# new data ----------------------------------------------------------------


deidentify$site <-factor(deidentify$site, levels = sitenames, labels = sitecoded)

deidentify$animal.id <-factor(deidentify$animal.id, 
                                levels = animalnames, labels = animalcoded)

deidentify$grid_id <- factor(deidentify$grid_id, 
                             levels = gridnames, labels = gridcoded)

sat$new_site <- factor(sat$new_site, levels = sitenames, labels = sitecoded)

sat$grid_id <- factor(sat$grid_id, 
                      levels = gridnames, labels = gridcoded)

surveys$site <- factor(surveys$site, levels = sub(' ', '', sitenames),
                       labels = sub(' ', '', sitecoded))



# filter ------------------------------------------------------------------

deidentify <- deidentify[deidentify$site!="Site H" & deidentify$site!="Site I",]


# save --------------------------------------------------------------------

write.csv(deidentify, "./data/ged_monitoring_grids2002-2023_De-identified.csv",
          row.names = F)

write.csv(sat, './data/ged_grids_satellite2000-2020_De-identified.csv',
          row.names = F)

write.csv(surveys, './data/ged_survey_numbers2002-2023_De-identified.csv',
          row.names = F)


