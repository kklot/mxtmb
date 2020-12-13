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

  // countries spatial
  DATA_IVECTOR     (cc_id);
  DATA_MATRIX      (R_cc);
  DATA_VECTOR      (sd_cc);
  PARAMETER_VECTOR (cc_vec);
  PARAMETER        (log_cc_e);
  Type cc_e = exp(log_cc_e);
  prior -= ktools::pc_prec(cc_e, sd_cc(0), sd_cc(1));
  prior -= ktools::soft_zero_sum(cc_vec);
  prior -= density::GMRF(ktools::prepare_Q(R_cc, cc_e))(cc_vec); 

  // ccxage interaction
  DATA_IVECTOR     (ccxage_id);
  DATA_MATRIX      (R_ccxage);
  PARAMETER_VECTOR (ccxage_vec);
  PARAMETER        (log_ccxage_e);
  Type ccxage_e = exp(log_ccxage_e);
  prior -= ktools::pc_prec(ccxage_e, sd_cc(0), sd_cc(1));
  prior -= ktools::constraint2D(ccxage_vec.data(), age_id.size(), cc_vec.size());
  prior -= density::GMRF(ktools::prepare_Q(R_ccxage, ccxage_e))(ccxage_vec); 

  // Data likelihood
  for (int i = 0; i < pna.size(); i++) {
    Type 
      mu = exp(beta0[0] + mu_sm[age[i]] + cc_vec[cc_id[i]] + ccxage_vec[ccxage_id[i]]),
      si = exp(beta0[1] + si_sm[age[i]]), 
      nu =     beta0[2] + nu_sm[age[i]], 
      ta = exp(beta0[3] + ta_sm[age[i]]);    
    dll -= dSHASHo(pna[i], mu, si, nu, ta, true);
  }
  Type llrp = dll;
  dll += prior;

  // Reporting
  int nC = cc_vec.size(), nA = age_id.size();
  vector<Type> si_vec(nA), nu_vec(nA), ta_vec(nA), mu_vec(nC * nA), rdims(2);
  rdims << nA, nC;
  for (int c = 0; c < nC; ++c)
    for (int a = 0; a < nA; ++a)
      mu_vec[c * nA + a] = exp( beta0[0] + mu_sm[age_id[a]] + cc_vec[c] + ccxage_vec[c * nA + a] );

  for (int a = 0; a < nA; ++a) {
      si_vec[a] = exp(beta0[1] + si_sm[age_id[a]]);
      nu_vec[a] =     beta0[2] + nu_sm[age_id[a]];
      ta_vec[a] = exp(beta0[3] + ta_sm[age_id[a]]);
  }

  REPORT(llrp);
  REPORT(prior);
  REPORT(prec_sm);
  REPORT(mu_sm); REPORT(si_sm); REPORT(nu_sm); REPORT(ta_sm);
  REPORT(beta0); REPORT(cc_vec); 
  REPORT(rdims);
  REPORT(age_id);
  REPORT(mu_vec); REPORT(si_vec); REPORT(nu_vec); REPORT(ta_vec);
  return dll;
}