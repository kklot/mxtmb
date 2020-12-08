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
  DATA_VECTOR(mu_beta0); DATA_VECTOR(sd_beta0);

  Type prior = 0.0;
  for (int i = 0; i < beta0.size(); ++i)
    prior -= dnorm(beta0[i], mu_beta0[i], sd_beta0[i], true);

  // Splines
  DATA_VECTOR(age_knots);
  PARAMETER_VECTOR(mu_sm);
  PARAMETER_VECTOR(si_sm);
  PARAMETER_VECTOR(nu_sm);
  PARAMETER_VECTOR(ta_sm);

  prior -= dnorm(mu_sm[0], Type(0), Type(1), true) -
    dnorm(si_sm[0], Type(0), Type(1), true) -
    dnorm(nu_sm[0], Type(0), Type(1), true) -
    dnorm(ta_sm[0], Type(0), Type(1), true);

  for (int i = 1; i < age_knots.size(); ++i) {
    mu_sm[i] = mu_sm[i-1] + mu_sm[i] * 0.1;
    si_sm[i] = si_sm[i-1] + si_sm[i] * 0.1;
    nu_sm[i] = nu_sm[i-1] + nu_sm[i] * 0.1;
    ta_sm[i] = ta_sm[i-1] + ta_sm[i] * 0.1;
  }
  
  // prior -= dnorm(mu_sm, Type(0), Type(1), true).sum() -
  //          dnorm(si_sm, Type(0), Type(1), true).sum() - 
  //          dnorm(nu_sm, Type(0), Type(1), true).sum() - 
  //          dnorm(ta_sm, Type(0), Type(1), true).sum();// estimate sd?

  tmbutils::splinefun<Type>
     mu_spline(age_knots, mu_sm, 1), 
     si_spline(age_knots, si_sm, 1),
     nu_spline(age_knots, nu_sm, 1),
     ta_spline(age_knots, ta_sm, 1); // natural spline

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
      age_i = age[i],
      mu    = exp(beta0[0] + mu_spline(age_i) + cc_vec[cc_id[i]]),
      si    = exp(beta0[1] + si_spline(age_i)),
      nu    =     beta0[2] + nu_spline(age_i),
      ta    = exp(beta0[3] + ta_spline(age_i)); 
    dll  -= dSHASHo(pna[i], mu, si, nu, ta, true);
  }
  dll += prior;

  // Reporting
  
  int nC = cc_vec.size(), nA = age_id.size();

  vector<Type> 
    mu_vec(nC * nA),
    si_vec(nA), 
    nu_vec(nA),
    ta_vec(nA),
    rdims(2);

  rdims << nA, nC;

  for (int i = 0; i < nA; ++i) {
    Type aa  = age_id[i];
    si_vec[i]    = exp(beta0[1] + si_spline(aa)); 
    nu_vec[i]    =     beta0[2] + nu_spline(aa); 
    ta_vec[i]    = exp(beta0[3] + ta_spline(aa)); 
  }

  for (int cc = 0; cc < nC; ++cc)
    for (int aa = 0; aa < nA; ++aa)
      mu_vec[cc * nA + aa] = exp(beta0[0] + mu_spline(Type(age_id[aa])) + cc_vec[cc]);

  REPORT(mu_sm); REPORT(si_sm); REPORT(nu_sm); REPORT(ta_sm);
  REPORT(beta0); REPORT(cc_vec); REPORT(age_id); REPORT(rdims);
  REPORT(mu_vec); REPORT(si_vec); REPORT(nu_vec); REPORT(ta_vec);
  return dll;
}