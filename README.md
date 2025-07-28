# Tympo_long-term_monitoring

manuscript - Adding value to imperfect long-term animal monitoring data for conservation insight

There are three analysis script files associated with this manuscript:

1.  long-term_monitoring.R
2.  short-term_surveys
3.  landsat_data.R

There are four data files associated with this analysis.

1.  ged_grids_satellite2000-2020.csv
2.  ged_monitoring_grids2002-2023.csv <sup>1</sup>
3.  ged_occupancy_surveys2018.csv\
4.  ged_survey_numbers2002-2023.csv <sup>1</sup>

<sup>1</sup> **Data supplied by “Australian Capital Territory (City and Environment Directorate, 2023)”**

and one output file included.

1.  N_estimates.csv

### Details for: ged_grids_satellite2000-2020.csv

Description: a comma-delimited file containing satellite imagery data.

Format(s): .csv

Size(s): 54.7 KB

Dimensions: 672 rows x 9 columns

-   Variables:

    -   objectid: index for grid_id
    -   oldgrid_id: old grid code
    -   year: year of dragon survey
    -   FVC: Fractional Vegetation Cover
    -   LST: Land Surface Temperature
    -   NDVI: Normalised Difference Vegetation Index
    -   area: metres squared
    -   grid_id: unique identifier for site and grid
    -   new_site: site name

-   Missing data codes: NA

### Details for: ged_monitoring_grids2002-2023.csv

Description: a comma-delimited file containing capture data from grassland earless dragon (ged) surveys

Format(s): .csv

Size(s): 284 KB

Dimensions: 2199 rows x 15 columns

-   Variables:

    -   date: dd/mm/yyyy
    -   site: site name
    -   year: year of survey
    -   time: 24hr time
    -   animal.id: unique animal identifier of captures within years and sites
    -   grid: grid code
    -   trap.id: pitfall trap code
    -   shelter: shelter covering pitfall trap (Y = yes, N = no)
    -   recap: whether capture is a recapture or new capture
    -   sex: dragon sex (female/F, male/M, blank = unknown/unsure)
    -   svl: snout-vent length (mm)
    -   weight: grams
    -   check: which survey check day the dragon was captured
    -   grid_new: new grid code
    -   grid_id: site code and grid code

-   Missing data codes: blank cell and NA

Data supplied by “Australian Capital Territory (City and Environment Directorate, 2023)”.

### Details for: ged_occupancy_surveys2018.csv

Description: a comma-delimited file containing occupancy data from grassland earless dragon (ged) natural burrow surveys

Format(s): .csv

Size(s): 50.9 KB

Dimensions: 500 rows x 21 columns

-   Variables:

    -   date: date of survey (yyyy-mm-dd)
    -   objID: 20 by 20 metre survey plot code
    -   Repeat: survey number (1-5)
    -   temp: temperature when plot was surveyed (Celsius - taken from BOM)
    -   site: site name
    -   nburrow: total number of burrows (all types)
    -   wsburrow: burrows with wolf spiders
    -   rcburrow: burrows with raspy crickets
    -   empty: empty burrows
    -   otherburrow: burrows with other insects
    -   aburrow: artificial burrows in the survey plot
    -   natburrow: number of natural burrows in the plot
    -   time: time of survey
    -   grid: whether the survey plot overlayed the long-term monitoring grids
    -   vegtype: vegetation type (native pasture or NTG = natural temperate grassland)
    -   vegstruct: vegetation structure
    -   bare: bare ground percentage
    -   dragons: number of dragons found in natural burrows
    -   outside: number of dragons found inside the plot but outside a natural burrow
    -   adragon: number of dragons found inside an artificial burrow
    -   ndragons: total number of dragons found in the plot
    -   padragons: presence/absence of dragons in a plot (1/0)

### Details for: ged_survey_numbers2002-2023.csv

Description: a comma-delimited file containing number of checks of long-term monitoring grids (survey effort)

Format(s): .csv

Size(s): 9.7 KB

Dimensions: 483 rows x 4 columns

-   Variables:

    -   site: site name
    -   grid: grid code
    -   year: year of dragon surveys
    -   checks: number of times monitoring grids were checked during the survey period

Data supplied by “Australian Capital Territory (City and Environment Directorate, 2023)”.

### Details for: N_estimates.csv

Description: a comma-delimited file containing abundance estimates calculated from the long-term_monitoring.R file. N_estimates.R is located in the output folder.

Format(s): .csv

Size(s): 24.4 KB

Dimensions: 277 rows x 6 columns

-   Variables:

    -   N: abundance estimate
    -   se: standard error
    -   lcl: lower 95% confidence interval
    -   ucl: upper 95% confidence interval
    -   year: year of dragon survey
    -   grid: site name and monitoring grid id
