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

  // Marginal likelihood
  vector<Type> u_m(N), u_w(N);
  for (int i = 0; i < N; i++) {
    Type
      mum = exp(beta0_m[0]), 
      muw = exp(beta0_w[0]),
      sim = exp(beta0_m[1]), 
      siw = exp(beta0_w[1]), 
      num =     beta0_m[2], 
      nuw =     beta0_w[2], 
      tam = exp(beta0_m[3]),
      taw = exp(beta0_w[3]);
    dll -= dSHASHo(w_age[i], muw, siw, nuw, taw, true);
    dll -= dSHASHo(m_age[i], mum, sim, num, tam, true);
    u_w[i] = pSHASHo(w_age[i], muw, siw, nuw, taw);
    u_m[i] = pSHASHo(m_age[i], mum, sim, num, tam);
  }
  // Copula
  PARAMETER(alpha);
  for (int i = 0; i < N; i++)
    dll -= ktools::dfrankCopula(u_w[i], u_m[i], alpha + cc_vec[cc_id[i]], true);
  
  Type llrp = dll;
  dll += prior;
  
  REPORT(prior); REPORT(llrp);
  REPORT(cc_e); 
  REPORT(alpha); 
  REPORT(beta0_w); REPORT(beta0_m); REPORT(cc_vec); 
  return dll;
}