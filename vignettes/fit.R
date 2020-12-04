library(ktools)

setwd('~/GitHub/mx_paper/')
roxygen2::roxygenise('mixtmb'); devtools::load_all('mixtmb')

dat <- fread('/Volumes/kklot/dhs/PA.csv.bz2')

mn_selfreport <- dodo(sx=1, dat, backward=FALSE)
mn_backreport <- dodo(sx=1, dat, backward=TRUE)
wm_selfreport <- dodo(sx=2, dat, backward=FALSE)
wm_backreport <- dodo(sx=2, dat, backward=TRUE)

mn_selfreport$obj$env$data %>% str
mn_backreport$obj$env$data %>% str
wm_selfreport$obj$env$data %>% str
wm_backreport$obj$env$data %>% str

dir.create('fit/skewlogis')
save.image('fit/skewlogis/separate.RData')
