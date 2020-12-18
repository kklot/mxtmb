#include <TMB.hpp>
#include "ktools.hpp"

template<class Type>
Type objective_function<Type>::operator() ()
{
  parallel_accumulator<Type> dll(this);
  Type prior = 0.0;

  // data
  DATA_VECTOR(m_age);
  DATA_VECTOR(w_age);

  int N = m_age.size();

  // Intercept
  PARAMETER_VECTOR(beta0_m);
  PARAMETER_VECTOR(beta0_w);
  DATA_VECTOR(mu_beta0); DATA_VECTOR(sd_beta0);
  for (int i = 0; i < beta0_m.size(); ++i) {
    prior -= dnorm(beta0_m[i], mu_beta0[i], sd_beta0[i], true);
    prior -= dnorm(beta0_w[i], mu_beta0[i], sd_beta0[i], true);
  }

  // countries spatial
  DATA_IVECTOR     (cc_id);
  DATA_MATRIX      (R_cc);
  DATA_VECTOR      (sd_cc);
  PARAMETER_VECTOR (cc_vec);
  PARAMETER        (log_cc_e);
  Type cc_e = exp(log_cc_e);
  prior -= ktools::pc_prec(cc_e, sd_cc(0), sd_cc(1));
  prior -= ktools::soft_zero_sum(cc_vec);
  prior += density::GMRF(ktools::prepare_Q(R_cc, cc_e))(cc_vec); 

  // iid mu
  PARAMETER_VECTOR(cc_mu_m);
  DATA_VECTOR(sd_iid);
  prior -= ktools::soft_zero_sum(cc_mu_m);
  prior -= dnorm(cc_mu_m, Type(0), sd_iid[0], true).sum();

  PARAMETER_VECTOR(cc_mu_w);
  prior -= ktools::soft_zero_sum(cc_mu_w);
  prior -= dnorm(cc_mu_w, Type(0), sd_iid[1], true).sum();

  PARAMETER_VECTOR(cc_si_m);
  prior -= ktools::soft_zero_sum(cc_si_m);
  prior -= dnorm(cc_si_m, Type(0), sd_iid[2], true).sum();

  PARAMETER_VECTOR(cc_si_w);
  prior -= ktools::soft_zero_sum(cc_si_w);
  prior -= dnorm(cc_si_w, Type(0), sd_iid[3], true).sum();

  // Copula
  PARAMETER(alpha);
  prior -= dnorm(alpha,  Type(8), Type(1), true);
  alpha += 1e-5; // avoid alpha to exactly be zero
  vector<Type> u_m(N), u_w(N);
  for (int i = 0; i < N; i++) {
    int j = cc_id[i];
    Type
      mum = exp(beta0_m[0] + cc_mu_m[j]), 
      muw = exp(beta0_w[0] + cc_mu_w[j]),
      sim = exp(beta0_m[1] + cc_si_m[j]), 
      siw = exp(beta0_w[1] + cc_si_w[j]),
      num =     beta0_m[2], 
      nuw =     beta0_w[2], 
      tam = exp(beta0_m[3]),
      taw = exp(beta0_w[3]);
    // Marginal likelihood
    dll -= dSHASHo(w_age[i], muw, siw, nuw, taw, true);
    dll -= dSHASHo(m_age[i], mum, sim, num, tam, true);
    // Copula
    Type
      u_w = pSHASHo(w_age[i], muw, siw, nuw, taw),
      u_m = pSHASHo(m_age[i], mum, sim, num, tam);
    dll -= ktools::dfrankCopula(u_w, u_m, alpha + cc_vec[j], true);
  }
  
  Type llrp = dll;
  dll += prior;
  vector<Type>
    a_v  = alpha + cc_vec.array(),
    mu_m = exp(beta0_m[0] + cc_mu_m.array()), 
    mu_w = exp(beta0_w[0] + cc_mu_w.array()), 
    si_m = exp(beta0_m[1] + cc_si_m.array()),
    si_w = exp(beta0_w[1] + cc_si_w.array());
  Type
    nu_m =     beta0_m[2],
    nu_w =     beta0_w[2],
    ta_m = exp(beta0_m[3]),
    ta_w = exp(beta0_w[3]);

  REPORT(prior); REPORT(llrp);
  REPORT(cc_e); 
  REPORT(alpha); 
  REPORT(beta0_w); REPORT(beta0_m); 
  REPORT(a_v);
  REPORT(mu_m); REPORT(mu_w);
  REPORT(si_m); REPORT(si_w);
  REPORT(nu_m); REPORT(nu_w);
  REPORT(ta_m); REPORT(ta_w);
  return dll;
}