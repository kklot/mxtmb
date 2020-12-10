#include <TMB.hpp>
#include "ktools.hpp"

template<class Type>
Type objective_function<Type>::operator() ()
{
  parallel_accumulator<Type> dll(this);

  // data
  DATA_VECTOR(pna);
  DATA_IVECTOR(age_id);

  DATA_MATRIX(X);
  DATA_SPARSE_MATRIX(S);

  DATA_IVECTOR(cc_id);
  DATA_MATRIX(R_cc);
  DATA_INTEGER(R_cc_rank);

  PARAMETER_VECTOR(beta0);
  DATA_VECTOR(mu_beta0); DATA_VECTOR(sd_beta0);

  Type prior = 0.0;
  for (int i = 0; i < beta0.size(); ++i)
    prior -= dnorm(beta0[i], mu_beta0[i], sd_beta0[i], true);

  // Splines
  PARAMETER_VECTOR(mu_sm);
  PARAMETER_VECTOR(si_sm);
  PARAMETER_VECTOR(nu_sm);
  PARAMETER_VECTOR(ta_sm);
  
  // random walk spline
  prior -= dnorm(diff(mu_sm), Type(0), Type(1), true).sum();
  prior -= dnorm(diff(si_sm), Type(0), Type(1), true).sum();
  prior -= dnorm(diff(nu_sm), Type(0), Type(1), true).sum();
  prior -= dnorm(diff(ta_sm), Type(0), Type(1), true).sum();

  // Penalty
  PARAMETER_VECTOR(la_sm);
  vector<Type> a_sm = exp(la_sm);
  prior -= ktools::p_spline(mu_sm, a_sm(0), S);
  prior -= ktools::p_spline(si_sm, a_sm(0), S);
  prior -= ktools::p_spline(nu_sm, a_sm(0), S);
  prior -= ktools::p_spline(ta_sm, a_sm(0), S);

  // countries spatial
  PARAMETER_VECTOR  (cc_vec);
  PARAMETER         (log_cc_e);
  DATA_VECTOR       (sd_cc);
  Type cc_e = exp(log_cc_e);
  prior -= ktools::pc_prec(cc_e, sd_cc(0), sd_cc(1));
  prior -= ktools::soft_zero_sum(cc_vec);
  prior -= density::GMRF(ktools::prepare_Q(R_cc, cc_e))(cc_vec); 
  // prior -= (R_cc_rank - cc_vec.size()) * log(sqrt(2*M_PI));

  // Data likelihood
  vector<Type> 
    mu =     beta0[0] + X * mu_sm, // model log(partner) see below
    si = exp(beta0[1] + X * si_sm), 
    nu =     beta0[2] + X * nu_sm, 
    ta = exp(beta0[3] + X * ta_sm);
  for (int i = 0; i < pna.size(); i++)
    dll  -= dSHASHo(pna[i], exp(mu[i] + cc_vec[cc_id[i]]), si[i], nu[i], ta[i], true);
  Type llrp = dll;
  dll += prior;

  // Reporting

  DATA_MATRIX(P);
  int nC = cc_vec.size(), nA = age_id.size();
  vector<Type>
    mu_v   =     beta0[0] + P * mu_sm, 
    si_vec = exp(beta0[1] + P * si_sm), 
    nu_vec =     beta0[2] + P * nu_sm, 
    ta_vec = exp(beta0[3] + P * ta_sm),
    mu_vec(nC * nA),
    rdims(2);
  rdims << nA, nC;
  for (int cc = 0; cc < nC; ++cc)
    for (int aa = 0; aa < nA; ++aa)
      mu_vec[cc * nA + aa] = exp(mu_v[aa] + cc_vec[cc]);

  REPORT(llrp);
  REPORT(prior);
  REPORT(a_sm);
  REPORT(mu_sm); REPORT(si_sm); REPORT(nu_sm); REPORT(ta_sm);
  REPORT(mu); REPORT(si); REPORT(nu); REPORT(ta);
  REPORT(beta0); REPORT(cc_vec); REPORT(age_id); REPORT(rdims);
  REPORT(mu_vec); REPORT(si_vec); REPORT(nu_vec); REPORT(ta_vec);
  return dll;
}