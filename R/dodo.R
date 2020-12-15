#' @export
dodo = function(cc, sx=1, dhs, backward=FALSE, fixpars = FALSE, n_cores=4,
    test=FALSE, rw_order=1) {

library(tidyverse)
library(magrittr)
set.seed(2020)

meta <- list() # track used metadata for post-processing
meta$seed <- 2020

# Read and prep data

dhs %<>% 
    mutate_if(is.integer, as.double) %>%
    mutate(yob=cmc_to_year(dob), age=cmc_to_year(doi)-yob) %>%
    mutate(sex =  if_else(sex=='male', 1, 2)) %>% 
    filter(partner < 60 & partner > 10)

# who report?
if (!backward) {
    dt = dhs %>% filter(sex==sx) %>% filter(ISO_A3 == cc)
} 
else {
    sy = ifelse(sx==1, 2, 1) # other sex
    dt = dhs %>%
        mutate(age2 = if_else(sex==sy, partner,     age)) %>% 
        mutate(par2 = if_else(sex==sy,     age, partner)) %>% 
        dplyr::select(ISO_A3, sex, age2, par2) %>% 
        rename(age = age2, partner=par2) %>%
        filter(sex == sy) %>% 
        filter(ISO_A3 == cc)
}

if (!nrow(dt)) return(NULL)

# Random walk
age_id   <- dt %>% pull(age) %>% range2seq
meta$min_age <- min(age_id)
n_age_id <- length(age_id)
R_age    <- genR(n_age_id, rw_order, scale=FALSE)

# TMB metadata and data
data = with(dt, 
    list(
        # meta
        age_id     = age_id - min(age),
        # response
        pna        = partner,
        age        = age - min(age),
        # random walk
        mu_beta0   = c( 3,   0,   0,   -1),
        sd_beta0   = c( 1,   1,   1,    1),
        rw_order   = rw_order,
        R_age      = R_age,
        sd_age     = c(.005, 0.01)
    )
)

init = list(
    beta0        = data$mu_beta0,
    mu_sm        = rep(0, n_age_id),
    si_sm        = rep(0, n_age_id),
    nu_sm        = rep(0, n_age_id),
    ta_sm        = rep(0, n_age_id),
    ln_sm        = rep(log(sd2prec(.005)), 4)
)

if (fixpars)
    fixpars = tmb_fixit(init, char(ln_sm))
else
    fixpars = NULL

opts = list(
    data       = data,
    parameters = init,
    random     = char(mu_sm, si_sm, nu_sm, ta_sm),
    silent     = 0,
    DLL        = 'mixtmb', 
    map        = fixpars
)

if (test)
    browser()

# Fit
library(TMB)
openmp(n_cores)
config(tape.parallel=0, optimize.instantly=1, DLL="mixtmb")

obj = do.call('MakeADFun', opts)
runSymbolicAnalysis(obj)
newtonOption(obj, smartsearch=TRUE)
fit = nlminb(obj$par, obj$fn, obj$gr)
rp  = obj$report()

o = list(obj=obj, fit=fit, meta=meta, rp=rp)
class(o) <- 'mixtmb'
o
}