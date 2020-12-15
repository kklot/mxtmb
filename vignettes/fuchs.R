# kklot/mixtmb@nospace
# 
library(ktools)
library(mixtmb)
library(data.table)

dat <- fread('data/PA.csv.bz2')

n_cores=20

dir.create("fit/indiv/", 0, recursive=1)

for (cc in unique(dat$ISO_A3)) {
    cat(cc, '\n')
    mn_selfreport <- mixtmb:::dodo(cc, sx=1, dat, backward=0, n_cores=n_cores)
    FreeADFun(mn_selfreport$obj)
    saveRDS(mn_selfreport, paste0('fit/indiv', cc, '_mn_selfreport.rds'))
    mn_backreport <- mixtmb:::dodo(cc, sx=1, dat, backward=1, n_cores=n_cores)
    FreeADFun(mn_backreport$obj)
    saveRDS(mn_backreport, paste0('fit/indiv', cc, '_mn_backreport.rds'))
    wm_selfreport <- mixtmb:::dodo(cc, sx=2, dat, backward=0, n_cores=n_cores)
    FreeADFun(wm_selfreport$obj)
    saveRDS(wm_selfreport, paste0('fit/indiv', cc, '_wm_selfreport.rds'))
    wm_backreport <- mixtmb:::dodo(cc, sx=2, dat, backward=1, n_cores=n_cores)
    FreeADFun(wm_backreport$obj)
    saveRDS(wm_backreport, paste0('fit/indiv', cc, '_wm_backreport.rds'))
    gc()
}