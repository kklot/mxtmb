library(ktools)
library(data.table)
library(mixtmb)

dat <- fread('data/PA.csv.bz2')

n_cores=20

dir.create("fit/sasrw_interaction_location", 0, recursive=1)

mn_selfreport <- mixtmb:::dodo(sx=1, dat, backward=FALSE, n_cores=n_cores)
FreeADFun(mn_selfreport$obj)
saveRDS(mn_selfreport, 'fit/sasrw_interaction_location/mn_selfreport.rds')
mn_backreport <- mixtmb:::dodo(sx=1, dat, backward=TRUE, n_cores=n_cores)
FreeADFun(mn_backreport$obj)
saveRDS(mn_backreport, 'fit/sasrw_interaction_location/mn_backreport.rds')
wm_selfreport <- mixtmb:::dodo(sx=2, dat, backward=FALSE, n_cores=n_cores)
FreeADFun(wm_selfreport$obj)
saveRDS(wm_selfreport, 'fit/sasrw_interaction_location/wm_selfreport.rds')
wm_backreport <- mixtmb:::dodo(sx=2, dat, backward=TRUE, n_cores=n_cores)
FreeADFun(wm_backreport$obj)
saveRDS(wm_backreport, 'fit/sasrw_interaction_location/wm_backreport.rds')