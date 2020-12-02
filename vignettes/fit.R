source("dodom.R")

library(ktools)
setwd('~/GitHub/mx_paper/')
devtools::load_all('mixtmb')
cd('../..')

mn_selfreport <- dodo(sx=1, backward=FALSE)
mn_backreport <- dodo(sx=1, backward=TRUE)
wm_selfreport <- dodo(sx=2, backward=FALSE)
wm_backreport <- dodo(sx=2, backward=TRUE)

mn_selfreport$obj$env$data %>% str
mn_backreport$obj$env$data %>% str
wm_selfreport$obj$env$data %>% str
wm_backreport$obj$env$data %>% str

dir.create('fit/skewlogis')
save.image('fit/skewlogis/separate.RData')
