to_array <- function(rp, par, iso3) {
    array(rp[[par]], rp$rdims) %>% 
        set_colnames(iso3) %>% 
        set_rownames(rp$age) %>% 
        as.data.table(1) %>% 
        rename(age=rn) %>% mutate(age=as.numeric(age)) %>% 
        pivot_longer(-age, names_to='ISO_A3', values_to=par)
}

dodo = function(sx=1, backward=FALSE, fixpars = FALSE, ICAR=TRUE, n_cores=4,
    test=FALSE, sub_set=0) {

library(ktools)
library(data.table)
library(tidyverse)
library(magrittr)
set.seed(2020)

meta <- list() # track used metadata for post-processing
meta$seed <- 2020

# Read and prep data
dhs <- fread('/Volumes/kklot/dhs/PA.csv.bz2')

dhs %<>% 
    mutate_if(is.integer, as.double) %>%
    mutate(yob=cmc_to_year(dob), age=cmc_to_year(doi)-yob) %>%
    mutate(sex =  if_else(sex=='male', 1, 2)) %>% 
    filter(partner < 60 & partner > 10)

# SSA only
ISO_SSA = name2iso(ktools:::.UN_SSA)

# who report?
if (!backward) {
    dt = dhs[sex==sx & ISO_A3 %in% ISO_SSA]
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
    dt = dt[, .SD[sample(.N, ifelse(.N > sub_set,  sub_set, .N))], ISO_A3]

# Country things
isoindata  <- unique(dt$ISO_A3)
n_cc_id    <- length(isoindata)
cc_id      <- data.table(ISO_A3 = isoindata, cc_id = 1:n_cc_id)
meta$cc_id <- cc_id
R_cc       <- diag(n_cc_id)
R_cc_rank  <- qr(R_cc)$rank # not really needed

dt[cc_id, on="ISO_A3", cc_id := i.cc_id]

# Spline
meta$num_knots <- 10
beta_knots <- dt[, sort(unique(age))] %>% 
    {seq(min(.), max(.), length.out=meta$num_knots)}

# TMB metadata and data
data = list(
    # data
    pna           = dt[, partner],
    log_age       = dt[, log(age)], 
    age_id        = dt[, sort(unique(age))],
    beta_knots    = beta_knots,
    mu_beta0      = c( 0,  5),
    sd_beta0      = c(.1,  1),
    mu_beta1      = c(-1,  1),
    sd_beta1      = c(.1, .1),
    sd_cc         = c(1, .1),
    cc_id         = dt[, cc_id-1],
    R_cc          = R_cc,
    R_cc_rank     = R_cc_rank
)

init = list(
    beta0        = data$mu_beta0,
    beta1        = data$mu_beta1,
    beta_sm      = rep(1, meta$num_knots),
    cc_vec       = rep(0, n_cc_id),
    log_cc_e     = log(sd2prec(1))
)

if (fixpars)
    fixpars = tmb_fixit(init, char(log_cc_e))
else
    fixpars = NULL

opts = list(
    data       = data,
    parameters = init,
    random     = char(beta0, beta1, cc_vec, beta_sm),
    silent     = 0,
    DLL        = 'mixtmb', 
    map        = fixpars
)

if (test) return(0)

# Fit
library(TMB)
openmp(n_cores)
config(tape.parallel=0, DLL="mixtmb")
obj = do.call('MakeADFun', opts)

obj$env$inner.control$tol10 = 0
fit = nlminb(obj$par, obj$fn, obj$gr)
# Report
rp  <- obj$report()
est <- to_array(rp, 'a_vec', isoindata) %>% 
    left_join(to_array(rp, 'b_vec', isoindata)) %>% 
    left_join(to_array(rp, 'g_vec', isoindata)) %>%
    arrange(ISO_A3)

o = list(obj=obj, fit=fit, meta=meta, rp=rp, est=est)
class(o) <- 'mixtmb'
o
}
