#' @export
dodo = function(fbi, fixpars = FALSE, ICAR=TRUE, n_cores=4, test=FALSE, sub_set=0, rw_order=1) {

library(tidyverse)
library(magrittr)
set.seed(2020)

meta <- list() # track used metadata for post-processing
meta$seed <- 2020

# Read and prep data

dt = fbi %>% 
  filter(age <= 80) %>% filter(partner <= 90) %>%
  mutate(w_age = if_else(sex=='male', partner, age)) %>% 
  mutate(m_age   = if_else(sex=='male', age, partner)) %>% 
  dplyr::select(ISO_A3, w_age, m_age, sex, direction)

# SSA only
dt = dt %>% filter(ISO_A3 %in% name2iso(ktools:::.UN_SSA))

if (sub_set > 0)
    dt %<>% group_by(ISO_A3) %>% sample_n(min(n(), sub_set))

# Country things
isoindata  <- unique(dt$ISO_A3)
n_cc_id    <- length(isoindata)
cc_id      <- data.table(ISO_A3 = isoindata, cc_id = 1:n_cc_id)
meta$cc_id <- cc_id
R_cc       <- diag(n_cc_id)
dt %<>% left_join(cc_id, 'ISO_A3')

# TMB metadata and data
data = with(dt, 
    list(
        # meta
        w_age = w_age,
        m_age = m_age,
        mu_beta0   = c( 2,  2,  1, .5),
        sd_beta0   = c(.02, .02, .02, .02),
        # spatial
        sd_cc      = c(.1, 0.1),
        cc_id      = cc_id - 1,
        R_cc       = R_cc
    )
)

init = list(
    beta0_w      = data$mu_beta0,
    beta0_m      = data$mu_beta0,
    alpha        = 8,
    cc_vec       = rep(0, n_cc_id), # alpha
    cc_mu_m      = rep(0, n_cc_id),
    cc_mu_w      = rep(0, n_cc_id),
    cc_si_m      = rep(0, n_cc_id),
    cc_si_w      = rep(0, n_cc_id),
    log_cc_e     = log(sd2prec(.1)),
    ln_sd_iid    = rep(log(.02),  4)
)

if (fixpars)
    fixpars = tmb_fixit(init, char(log_cc_e))
else
    fixpars = NULL

opts = list(
    data       = data,
    parameters = init,
    random     = char(cc_vec, cc_mu_m, cc_mu_w, cc_si_m, cc_si_w),
    silent     = 0,
    DLL        = 'mixtmb', 
    map        = fixpars
)

if (test)
    browser()

# Fit
library(TMB)
TMB::openmp(1)
TMB::config(tape.parallel=0, optimize.instantly=1, DLL="mixtmb")

openmp(n_cores)
config(tape.parallel=0, optimize.instantly=1, DLL="mixtmb")

obj = do.call('MakeADFun', opts)
runSymbolicAnalysis(obj)
fit = nlminb(obj$par, obj$fn, obj$gr)
rp  = obj$report()

o = list(obj=obj, fit=fit, meta=meta, rp=rp)
class(o) <- 'mixtmb'
o
}