library(ktools)

setwd('~/GitHub/mx_paper/')
roxygen2::roxygenise('mixtmb'); devtools::load_all('mixtmb')

fbi <- fread('data/bi_report.csv.bz2') %>% 
    filter(direction=='self')

cop <- dodo(fbi)