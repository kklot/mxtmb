#include <TMB.hpp>
#include "ktools.hpp"

template<class Type>
Type objective_function<Type>::operator() ()
{
  parallel_accumulator<Type> dll(this);

  // data
  DATA_VECTOR(pna);
  DATA_VECTOR(age);
  DATA_IVECTOR(age_id);

  DATA_IVECTOR(cc_id);
  DATA_MATRIX(R_cc);
  DATA_INTEGER(R_cc_rank);

  PARAMETER_VECTOR(beta0);
  PARAMETER_VECTOR(beta1);  
  // priors
  DATA_VECTOR(mu_beta0); DATA_VECTOR(sd_beta0);
  DATA_VECTOR(mu_beta1); DATA_VECTOR(sd_beta1);

  Type prior = 0.0;
  for (int i = 0; i < beta0.size(); ++i)
    prior -= dnorm(beta0[i], mu_beta0[i], sd_beta0[i], true);
  prior -= dnorm(beta1[0], mu_beta1[0], sd_beta1[0], true);

  // Splines
  DATA_VECTOR(age_knots);
  PARAMETER_VECTOR(beta_sm);
  PARAMETER_VECTOR(alpha_sm);
  
  prior -= dnorm(beta_sm, Type(0), Type(1), true).sum() -
           dnorm(alpha_sm, Type(0), Type(1), true).sum(); // estimate sd?

  tmbutils::splinefun<Type>
    beta_spline(age_knots, beta_sm), 
    alpha_spline(age_knots, alpha_sm);

  // countries spatial
  PARAMETER_VECTOR  (cc_vec);
  PARAMETER         (log_cc_e);
  DATA_VECTOR       (sd_cc);
  Type cc_e = exp(log_cc_e);
  prior -= ktools::pc_prec(cc_e, sd_cc(0), sd_cc(1));
  prior -= ktools::soft_zero_sum(cc_vec);
  prior -= density::GMRF(ktools::prepare_Q(R_cc, cc_e))(cc_vec); 
  prior -= (R_cc_rank - cc_vec.size()) * log(sqrt(2*M_PI));
  
  // Data likelihood
  for (int i = 0; i < pna.size(); i++) {
    Type 
      alpha = exp(beta0[0] + alpha_spline(Type(age[i])) + cc_vec[cc_id[i]]),
      beta  = exp(beta0[1] +  beta_spline(Type(age[i]))),
      gamma = exp(beta0[2] + beta1[0] * log(beta) );
    dll  -= log(ktools::ft_llogisI(pna(i), beta, alpha, gamma));
  }
  dll += prior;

  // Reporting
  int nC = cc_vec.size(), nA = age_id.size();
  vector<Type> rdims(2), a_vec(nC * nA), b_vec(nA), g_vec(nA);
  rdims << nA, nC;

  for (int i = 0; i < nA; ++i) {
    b_vec[i] = exp(beta0[1] +  beta_spline(Type(age_id[i])));
    g_vec[i] = exp(beta0[2] + beta1[0] * log(b_vec[i]));
  }

  for (int cc = 0; cc < nC; ++cc)
    for (int aa = 0; aa < nA; ++aa)
      a_vec[cc * nA + aa] = exp(beta0[0] + alpha_spline(Type(age_id[aa])) + cc_vec[cc]);

  REPORT(beta_sm); REPORT(alpha_sm);
  REPORT(beta0); REPORT(beta1); REPORT(cc_vec);
  REPORT(age_id); REPORT(rdims);
  REPORT(a_vec); REPORT(b_vec); REPORT(g_vec);
  return dll;
}