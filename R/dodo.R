#' @export
dodo = function(sx=1, dhs, backward=FALSE, fixpars = FALSE, ICAR=TRUE, n_cores=4,
    test=FALSE, sub_set=0, rw_order=1) {

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

# SSA only
ISO_SSA = name2iso(ktools:::.UN_SSA)

# who report?
if (!backward) {
    dt = dhs %>% filter(sex==sx & ISO_A3 %in% ISO_SSA)
} 
else {
    sy = ifelse(sx==1, 2, 1) # other sex
    dt = dhs %>%
        mutate(age2 = if_else(sex==sy, partner,     age)) %>% 
        mutate(par2 = if_else(sex==sy,     age, partner)) %>% 
        dplyr::select(ISO_A3, sex, age2, par2) %>% 
        rename(age = age2, partner=par2) %>%
        filter(sex == sy)
}

if (sub_set > 0)
    dt %<>% group_by(ISO_A3) %>% sample_n(min(n(), sub_set))

# Country things
isoindata  <- unique(dt$ISO_A3)
n_cc_id    <- length(isoindata)
cc_id      <- data.table(ISO_A3 = isoindata, cc_id = 1:n_cc_id)
meta$cc_id <- cc_id
R_cc       <- diag(n_cc_id)
R_cc_rank  <- qr(R_cc)$rank # not really needed

dt %<>% left_join(cc_id, 'ISO_A3')

age_id <- dt %>% pull(age) %>% unique %>% sort

# Random walk
n_age_id <- length(age_id)
R_age    <- genR(n_age_id, rw_order)

# TMB metadata and data
data = with(dt, 
    list(
        # meta
        age_id     = age_id - min(age),
        # response
        pna        = partner,
        age        = age - min(age),
        # random walk
        mu_beta0   = c( 3,   0,   0,    0),
        sd_beta0   = c( 1,   1,   1,    1),
        rw_order   = rw_order,
        R_age      = R_age,
        # spatial
        sd_age     = c(1, 0.1),
        sd_cc      = c(1, 0.1),
        cc_id      = cc_id - 1,
        R_cc       = R_cc, 
        R_cc_rank  = R_cc_rank
    )
)

init = list(
    beta0        = data$mu_beta0,
    mu_sm        = rep(0, n_age_id),
    si_sm        = rep(0, n_age_id),
    nu_sm        = rep(0, n_age_id),
    ta_sm        = rep(0, n_age_id),
    ln_sm        = rep(log(sd2prec(1)), 4),
    cc_vec       = rep(0, n_cc_id),
    log_cc_e     = log(sd2prec(10))
)

if (fixpars)
    fixpars = tmb_fixit(init, char(ln_sm))
else
    fixpars = NULL

opts = list(
    data       = data,
    parameters = init,
    random     = char(beta0, cc_vec, mu_sm, si_sm, nu_sm, ta_sm, ln_sm),
    silent     = 0,
    DLL        = 'mixtmb', 
    map        = fixpars
)

if (test)
    browser()

# Fit
library(TMB)
openmp(n_cores)
config(tape.parallel=0, optimize.instantly=0, DLL="mixtmb")

obj = do.call('MakeADFun', opts)
obj$env$inner.control$tol10 = 0
obj$env$tracepar = FALSE
fit = nlminb(obj$par, obj$fn, obj$gr)

rp  = obj$report()

est <- to_array(rp, 'mu_vec', isoindata) %>% 
    left_join(to_array(rp, 'si_vec', isoindata)) %>% 
    left_join(to_array(rp, 'nu_vec', isoindata)) %>%
    left_join(to_array(rp, 'ta_vec', isoindata)) %>%
    arrange(ISO_A3)

o = list(obj=obj, fit=fit, meta=meta, rp=rp, est=est)
class(o) <- 'mixtmb'
o
}