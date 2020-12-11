#' @export
dodo = function(sx=1, dhs, backward=FALSE, fixpars = FALSE, ICAR=TRUE, n_cores=4,
    test=FALSE, sub_set=0) {

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

# Spline
gam_S <- mgcv::gam(partner ~ -1 + s(age, bs = 'cs'), data=dt, fit=FALSE)
S     <- gam_S$smooth[[1]]$S[[1]] %>% as('sparseMatrix')
X     <- gam_S$X
P     <- mgcv::PredictMat(gam_S$smooth[[1]], data = data.frame(age=age_id))
num_basis <- nrow(S)

# TMB metadata and data
data = with(dt, 
    list(
        # meta
        age_id     = age_id,
        # response
        pna        = partner,
        # splines and design matrix (gamlss)
        mu_beta0   = c( 2,  -2,   0, -0.5),
        sd_beta0   = c(.1,  .1, 0.1,  0.1),
        S          = as(S, 'sparseMatrix'),
        X          = X,
        P          = P,
        rw_sd      = 0.1, 
        # spatial
        sd_cc      = c(1, 0.1),
        cc_id      = cc_id-1,
        R_cc       = R_cc,
        R_cc_rank  = R_cc_rank
    )
)

init = list(
    beta0        = data$mu_beta0,
    mu_sm        = rep(0, num_basis),
    si_sm        = rep(0, num_basis),
    nu_sm        = rep(0, num_basis),
    ta_sm        = rep(0, num_basis),
    la_sm        = rep(5, 1), # penalize splines/variance of spline coef.
    cc_vec       = rep(0, n_cc_id),
    log_cc_e     = log(sd2prec(1))
)

if (fixpars)
    fixpars = tmb_fixit(init, char(la_sm))
else
    fixpars = NULL

opts = list(
    data       = data,
    parameters = init,
    random     = char(beta0, cc_vec, mu_sm, si_sm, nu_sm, ta_sm, la_sm),
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
est <- to_array(rp, 'mu_vec', isoindata) %>% 
    left_join(to_array(rp, 'si_vec', isoindata)) %>% 
    left_join(to_array(rp, 'nu_vec', isoindata)) %>%
    left_join(to_array(rp, 'ta_vec', isoindata)) %>%
    arrange(ISO_A3)

o = list(obj=obj, fit=fit, meta=meta, rp=rp, est=est)
class(o) <- 'mixtmb'
o
}