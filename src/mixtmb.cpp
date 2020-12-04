#include <TMB.hpp>
#include "ktools.hpp"

template<class Type>
Type objective_function<Type>::operator() ()
{
  parallel_accumulator<Type> dll(this);

  // data
  DATA_VECTOR(pna);
  DATA_VECTOR(log_age);
  DATA_IVECTOR(age_id);

  DATA_IVECTOR(cc_id);
  DATA_MATRIX(R_cc);
  DATA_INTEGER(R_cc_rank);

  // Data model - log-logistic parameters
  PARAMETER_VECTOR(beta0);
  PARAMETER_VECTOR(beta1);
  
  DATA_VECTOR(beta_knots);
  PARAMETER_VECTOR(beta_sm);
  tmbutils::splinefun<Type> beta_spline(beta_knots, beta_sm);

  // priors
  DATA_VECTOR(mu_beta0);
  DATA_VECTOR(mu_beta1);
  DATA_VECTOR(sd_beta0);
  DATA_VECTOR(sd_beta1);

  Type prior = 0.0;
  for (int i = 0; i < beta0.size(); ++i) {
    prior -= dnorm(beta0[i], mu_beta0[i], sd_beta0[i], true);
    prior -= dnorm(beta1[i], mu_beta1[i], sd_beta1[i], true);
  }

  // countries spatial
  PARAMETER_VECTOR  (cc_vec);
  PARAMETER         (log_cc_e);
  DATA_VECTOR       (sd_cc);
  Type cc_e = exp(log_cc_e);
  prior -= ktools::pc_prec(cc_e, sd_cc(0), sd_cc(1));
  prior -= ktools::soft_zero_sum(cc_vec);
  prior += density::GMRF(ktools::prepare_Q(R_cc, cc_e))(cc_vec);
  prior += (R_cc_rank - cc_vec.size()) * log(sqrt(2*M_PI)); // ktools::GMRF would be nice
  
  // Data likelihood
  for (int i = 0; i < pna.size(); i++) {
    Type 
      eta   = beta0[0] + beta1[0] * log_age(i) + cc_vec(cc_id(i)),
      alpha = exp(eta), 
      beta  = exp(beta_spline(Type(exp(log_age[i])))),
      gamma = exp(beta0[1] + beta1[1] * log(beta));
      dll  -= log(ktools::ft_llogisI(pna(i), beta, alpha, gamma));
  }
  dll += prior;
  // Reporting
  int nC = cc_vec.size(), nA = age_id.size();
  vector<Type> rdims(2), a_vec(nC * nA), b_vec(nC * nA), g_vec(nC * nA), beta_age(nA);
  rdims << nA, nC;

  for (int i = 0; i < nA; ++i)
    beta_age[i] = exp(beta_spline(Type(age_id[i])));

  for (int cc = 0; cc < nC; ++cc) {
    for (int aa = 0; aa < nA; ++aa) {
      Type ate = beta0[0] + beta1[0] * log(age_id[aa]) + cc_vec[cc];
      a_vec[cc * nA + aa] = exp(ate);
      b_vec[cc * nA + aa] = beta_age[aa];
      g_vec[cc * nA + aa] = exp(beta0[1] + beta1[1] * log(beta_age[aa]));
    }
  }
  REPORT(beta_age);
  REPORT(beta0); REPORT(beta1); REPORT(cc_vec);
  REPORT(age_id); REPORT(rdims);
  REPORT(a_vec); REPORT(b_vec); REPORT(g_vec);
  return dll;
}