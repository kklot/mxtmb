#include <TMB.hpp>
#include "ktools.hpp"

template<class Type>
Type objective_function<Type>::operator() ()
{
  parallel_accumulator<Type> dll(this);
  Type prior = 0.0;

  // data
  DATA_VECTOR(pna);
  DATA_IVECTOR(age);
  DATA_IVECTOR(age_id);

  // Intercept
  PARAMETER_VECTOR(beta0);
  DATA_VECTOR(mu_beta0); DATA_VECTOR(sd_beta0);
  for (int i = 0; i < beta0.size(); ++i)
    prior -= dnorm(beta0[i], mu_beta0[i], sd_beta0[i], true);

  // Random walk
  PARAMETER_VECTOR(mu_sm);
  PARAMETER_VECTOR(si_sm);
  PARAMETER_VECTOR(nu_sm);
  PARAMETER_VECTOR(ta_sm);

  PARAMETER_VECTOR(ln_sm);
  DATA_VECTOR(sd_age);
  vector<Type> prec_sm = exp(ln_sm);
  for (int i = 0; i < prec_sm.size(); ++i)
    prior -= ktools::pc_prec(prec_sm[i], sd_age(0), sd_age(1));
  
  DATA_MATRIX(R_age);
  DATA_SCALAR(rw_order);
  prior += ktools::rw(mu_sm, R_age, prec_sm[0], rw_order);
  prior += ktools::rw(si_sm, R_age, prec_sm[1], rw_order);
  prior += ktools::rw(nu_sm, R_age, prec_sm[2], rw_order);
  prior += ktools::rw(ta_sm, R_age, prec_sm[3], rw_order);

  // Data likelihood
  for (int i = 0; i < pna.size(); i++) {
    Type 
      mu = exp(beta0[0] + mu_sm[age[i]]),
      si = exp(beta0[1] + si_sm[age[i]]), 
      nu =     beta0[2] + nu_sm[age[i]], 
      ta = exp(beta0[3] + ta_sm[age[i]]);    
    dll -= dSHASHo(pna[i], mu, si, nu, ta, true);
  }
  Type llrp = dll;
  dll += prior;

  // Reporting
  int nA = age_id.size();
  vector<Type> si_vec(nA), nu_vec(nA), ta_vec(nA), mu_vec(nA);

  for (int a = 0; a < nA; ++a) {
      mu_vec[a] = exp(beta0[0] + mu_sm[age_id[a]]);
      si_vec[a] = exp(beta0[1] + si_sm[age_id[a]]);
      nu_vec[a] =     beta0[2] + nu_sm[age_id[a]];
      ta_vec[a] = exp(beta0[3] + ta_sm[age_id[a]]);
  }

  REPORT(llrp);
  REPORT(prior);
  REPORT(prec_sm);
  REPORT(mu_sm); REPORT(si_sm); REPORT(nu_sm); REPORT(ta_sm);
  REPORT(beta0);
  REPORT(age_id);
  REPORT(mu_vec); REPORT(si_vec); REPORT(nu_vec); REPORT(ta_vec);
  return dll;
}