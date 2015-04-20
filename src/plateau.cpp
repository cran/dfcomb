#include<R.h>

#define BOOST_DISABLE_ASSERTS
#define ARMA_32BIT_WORD
#define ARMA_DONT_USE_BLAS
#define ARMA_DONT_USE_LAPACK
#include <RcppArmadillo.h>

#undef dnorm

#include <cmath>
#include <string>
#include <functional>
#include <tuple>

#include <cppbugs/cppbugs.hpp>

#include <boost/random/exponential_distribution.hpp>
#include <cstdlib>

#ifdef _OPENMP
#include<omp.h>
#endif

#include <progress.hpp>

extern "C" {
#include "dfcomb.h"
}

#include<cstdio>

using namespace std;
using namespace arma;
using namespace cppbugs;

namespace dfcomb { namespace plateau {

double TARG_SUP, EFF_MIN;
double TIMEFULL, CYCLE;
int COHORT_START, COHORT;

template<typename T>
typename T::value_type median(const T& x) {
  if(x.size() == 0)
    return 0;
  std::vector<typename T::value_type> v(x.begin(), x.end());
  std::nth_element(v.begin(), v.begin()+v.size()/2, v.end());
  if(v.size() % 2 == 1)
    return v[v.size()/2];
  else
    return (v[v.size()/2] +
            *std::max_element(v.begin(), v.begin()+v.size()/2))/2.;
}

/* True toxicity and response rates from input */
struct true_data {
  vector<vector<double>> piV, respV;
  double incl_per_week;
};

/* A struct representing the parameters of the toxicity model */
struct toxicity_parameters {
  double beta0, beta1, beta2;

  toxicity_parameters(double beta0, double beta1, double beta2) : beta0(beta0), beta1(beta1), beta2(beta2) {}
  toxicity_parameters() : beta0(0), beta1(0), beta2(0) {}

  template<typename U>
  void proba_tox(const U& dose_tox1, const U& dose_tox2, U& out) {
    out = 1-1/(1+exp_approx(beta0 + beta1*dose_tox1 + beta2*dose_tox2));
  }
};

/* A struct representing the parameters of the progression model */
struct progression_parameters {
  double lambda0, gamma1, gamma2;
  int tau;

  progression_parameters(double lambda0, double gamma1, double gamma2, int tau) :
    lambda0(lambda0), gamma1(gamma1), gamma2(gamma2), tau(tau) {}
  progression_parameters() :
    lambda0(0), gamma1(0), gamma2(0), tau(0) {}

  template<typename U>
  void calc_hasard(const U& dose_prog1, const vector<U>& dose_prog2_tau, U& out) const {
    unsigned tau = min((int)dose_prog2_tau.size()-1, max(0, this->tau));
    out = lambda0 * exp_approx(gamma1*dose_prog1 + gamma2*dose_prog2_tau[tau]);
  }

  double survival(double t, double d1, const vector<double>& d2) const {
    double hasard;
    calc_hasard(d1, d2, hasard);
    return exp_approx(-t*hasard);
  }

  double responseRate(double d1, const vector<double>& d2) const {
    return survival(TIMEFULL, d1, d2);
  }
};

struct estimations {
  /* Estimated parameters */
  toxicity_parameters tox_params;
  progression_parameters prog_params;

  /* Estimated information on doses */
  vector<vector<double>> pi, ptox_sup, resp, qeff_min;
  vector<double> proba_tau;

  estimations(int ndose1, int ndose2):
    pi(ndose1, vector<double>(ndose2)),
    ptox_sup(ndose1, vector<double>(ndose2)),
    resp(ndose1, vector<double>(ndose2)),
    qeff_min(ndose1, vector<double>(ndose2)),
    proba_tau(ndose2)
  { }
};

enum trial_stage { FSTHALF, SNDHALF, FINAL };

/* All the data needed dring a trial */
struct trial_data{
  vector<double> dose1T, dose2T;
  vector<double> dose1E, dose2E;
  int cdose1, cdose2;
  int startup_end;
  double time_cur;
  int pat_incl;
  vector<unsigned> dose_adm1, dose_adm2;
  vector<double> time_prog, time_incl;
  vector<int> toxicity;
  trial_stage stage;

  std::mt19937_64 r;

  trial_data(const vector<double>& dose1TV, const vector<double>& dose2TV,
             const vector<double>& dose1EV, const vector<double>& dose2EV,
             uint_fast64_t seed):
    dose1T(dose1TV.size()), dose2T(dose2TV.size()),
    dose1E(dose1EV.size()), dose2E(dose2EV.size()),
    cdose1(0), cdose2(0), startup_end(-1), time_cur(0), pat_incl(0),
    stage(FSTHALF), r(seed)
  {
    for(int i=0; i<dose1TV.size(); i++){
      dose1T[i] = log(dose1TV[i]/(1-dose1TV[i]));
      dose1E[i] = -log(-log(dose1EV[i])/TIMEFULL);
    }
    for(int i=0; i<dose2TV.size(); i++){
      dose2T[i] = log(dose2TV[i]/(1-dose2TV[i]));
      dose2E[i] = -log(-log(dose2EV[i])/TIMEFULL);
    }
  }
};

/* Results of a set of trials */
struct results {
  int inconc;
  vector<vector<int>> nptsdose, nrecdose;
  vector<vector<int>> ntox, neff;
  int n_trials;
  double duration;

  results(int ndose1, int ndose2) :
    nptsdose(ndose1, vector<int>(ndose2, 0)),
    nrecdose(ndose1, vector<int>(ndose2, 0)),
    ntox(ndose1, vector<int>(ndose2, 0)),
    neff(ndose1, vector<int>(ndose2, 0)),
    inconc(0), n_trials(0),
    duration(0)
  { }

  void accumulate(const trial_data& trial_data, const pair<int, int>& recom) {
    if(recom.first == -1)
      inconc++;
    else
      nrecdose[recom.first][recom.second]++;
    for(int pat = 0; pat < trial_data.pat_incl; pat++) {
      int dose1 = trial_data.dose_adm1[pat];
      int dose2 = trial_data.dose_adm2[pat];
      if(trial_data.time_prog[pat] > TIMEFULL)
        neff[dose1][dose2]++;
      if(trial_data.toxicity[pat])
        ntox[dose1][dose2]++;
      nptsdose[dose1][dose2]++;
    }
    n_trials++;
    duration += trial_data.time_cur;
  }
};

// Wait for the next patient to arrive. Rejects it if we did not spend
// CYCLE week since the last inclusion in this cohort. Returns true if the
// trial ends.
bool wait_patient(trial_data& trial_data, const true_data& true_data) {
  boost::exponential_distribution<double> exp_rng;

  // The trial is finished.
  if(trial_data.cdose1 < 0) {
    if(trial_data.pat_incl > 0 &&
       trial_data.time_cur -
       trial_data.time_incl[trial_data.pat_incl-1] <= TIMEFULL)
      trial_data.time_cur = TIMEFULL + 0.01 + trial_data.time_incl[trial_data.pat_incl-1];
    return true;
  } else {
    while(true) {
      /* How much time do we wait for this patient ? */
      double time_incl = exp_rng(trial_data.r, true_data.incl_per_week);

      /* We shift everything by this amount of time. */
      trial_data.time_cur += time_incl;

      /* Should we reject this patient ? */
      /* Because we do not already know whether it is toxic for the latest cohort ? */
      int startup_end = trial_data.startup_end;
      if(!((startup_end == -1 && trial_data.pat_incl % COHORT_START == 0) ||
           (startup_end != -1 && (trial_data.pat_incl - startup_end) % COHORT == 0)))
        break;

      if(trial_data.pat_incl == 0)
        break;

      if(trial_data.time_cur - trial_data.time_incl[trial_data.pat_incl-1] >= CYCLE)
        break;
    }
    return false;
  }
}

void include_patient(trial_data& trial_data, const true_data& true_data) {
  trial_data.dose_adm1.push_back(trial_data.cdose1);
  trial_data.dose_adm2.push_back(trial_data.cdose2);

  /* Is there a toxicity ? */
  double pi_ij = true_data.piV[trial_data.cdose1][trial_data.cdose2];
  std::uniform_real_distribution<double> uni_rng;
  trial_data.toxicity.push_back(uni_rng(trial_data.r) < pi_ij);

  /* When will there be a progression ? */
  double resp_ij = true_data.respV[trial_data.cdose1][trial_data.cdose2];
  boost::exponential_distribution<double> exp_rng;
  double param_exp = -log(resp_ij)/TIMEFULL;
  trial_data.time_prog.push_back(exp_rng(trial_data.r, param_exp));

  trial_data.time_incl.push_back(trial_data.time_cur);

  trial_data.pat_incl++;
}

/* Estimate toxicity and efficacy parameters. */
estimations estimate(const trial_data& trial_data) {
  const int ndose1 = trial_data.dose1T.size();
  const int ndose2 = trial_data.dose2T.size();

  /* Read usefull data from trial_data */
  int NPatients = trial_data.pat_incl;

  const vec toxicity = conv_to<vec>::from(trial_data.toxicity);

  vec dose_tox1 = vec(trial_data.dose1T).elem(conv_to<uvec>::from(trial_data.dose_adm1));
  vec dose_tox2 = vec(trial_data.dose2T).elem(conv_to<uvec>::from(trial_data.dose_adm2));

  vec dose_prog1 = vec(trial_data.dose1E).elem(conv_to<uvec>::from(trial_data.dose_adm1));
  vector<vec> dose_prog2_tau;
  for(int tau = 0; tau < ndose2; tau++)
    dose_prog2_tau.push_back(
      vec(trial_data.dose2E).elem(clamp(conv_to<uvec>::from(trial_data.dose_adm2), 0, tau)));

  vec time_min_prog(NPatients), progression(NPatients);
  for(unsigned pat = 0; pat < NPatients; pat++) {
    double time_follow =
      min(TIMEFULL, trial_data.time_cur - trial_data.time_incl[pat]);
    if(trial_data.time_prog[pat] < time_follow) {
      progression(pat) = 1;
      time_min_prog(pat) = trial_data.time_prog[pat];
    } else {
      progression(pat) = 0;
      time_min_prog(pat) = time_follow;
    }
  }

  /* MCMC initalization and sampling */
  progression_parameters prog_params(1, -1, -1, 0);
  toxicity_parameters tox_params(0, 1, 1);

  double gamma1Opp = 1, gamma2Opp = 1;

  vec hasard;
  vec pi_jk;
  std::function<void ()> model = [&]() {
    prog_params.gamma1 = -gamma1Opp;
    prog_params.gamma2 = -gamma2Opp;
    prog_params.calc_hasard(dose_prog1, dose_prog2_tau, hasard);
    tox_params.proba_tox(dose_tox1, dose_tox2, pi_jk);
  };
  MCModel<std::mt19937> m(model);

  m.track<Exponential>(prog_params.lambda0).dexp(0.01);
  m.track<Exponential>(gamma1Opp).dexp(0.1);
  m.track<Exponential>(gamma2Opp).dexp(0.1);
  m.track<Deterministic>(prog_params.gamma1);
  m.track<Deterministic>(prog_params.gamma2);
  vec distr(ndose2);
  for(int i = 0; i < ndose2; i++)
    distr(i) = pow(1.4,i);
  m.track<Discrete>(prog_params.tau).ddiscr(distr);

  m.track<ObservedExponentialCensored>(make_pair(time_min_prog, progression)).dexpcens(hasard);

  m.track<ObservedBernoulli>(toxicity).dbern(pi_jk);
  m.track<Normal>(tox_params.beta0).dnorm(0, 0.01);
  m.track<Exponential>(tox_params.beta1).dexp(1);
  m.track<Exponential>(tox_params.beta2).dexp(1);

  m.sample(1e6, 1e5, 1e4, 10);

  /* Get list of prog_params and tox_params */
  vector<pair<progression_parameters, toxicity_parameters> > params_draws_list;
  auto itbeta0 = m.getNode(tox_params.beta0).history.begin(), itbeta0end = m.getNode(tox_params.beta0).history.end();
  auto itbeta1 = m.getNode(tox_params.beta1).history.begin();
  auto itbeta2 = m.getNode(tox_params.beta2).history.begin();
  auto itlambda0 = m.getNode(prog_params.lambda0).history.begin();
  auto itgamma1 = m.getNode(prog_params.gamma1).history.begin();
  auto itgamma2 = m.getNode(prog_params.gamma2).history.begin();
  auto ittau = m.getNode(prog_params.tau).history.begin();

  while(itbeta0 != itbeta0end) {
    params_draws_list.emplace_back(
      progression_parameters(*itlambda0, *itgamma1, *itgamma2, *ittau),
      toxicity_parameters(*itbeta0, *itbeta1, *itbeta2));
    itbeta0++; itbeta1++; itbeta2++;
    itlambda0++; itgamma1++; itgamma2++; ittau++;
  }

  int n_samp = params_draws_list.size();

  /* Estimations */
  estimations result(ndose1, ndose2);
  result.tox_params.beta0 = median(m.getNode(tox_params.beta0).history);
  result.tox_params.beta1 = median(m.getNode(tox_params.beta1).history);
  result.tox_params.beta2 = median(m.getNode(tox_params.beta2).history);

  vector<int> counts(ndose2, 0);
  for(int tau : m.getNode(prog_params.tau).history)
    counts[tau]++;
  for(int d2 = 0; d2 < ndose2; d2++)
    result.proba_tau[d2] = (double)counts[d2]/n_samp;
  result.prog_params.tau = max_element(counts.begin(), counts.end())-counts.begin();

  result.prog_params.lambda0 = median(m.getNode(prog_params.lambda0).history);
  result.prog_params.gamma1 = median(m.getNode(prog_params.gamma1).history);
  result.prog_params.gamma2 = median(m.getNode(prog_params.gamma2).history);

  for(int d1 = 0; d1 < ndose1; d1++){
    for(int d2 = 0; d2 < ndose2; d2++) {
      int count_tox = 0, count_eff = 0;
      std::vector <double> pi_median, q_median;
      vector<double> dose2E_tau;
      for(int tau = 0; tau < ndose2; tau++)
        dose2E_tau.push_back(trial_data.dose2E[min(tau, d2)]);

      for(auto draw : params_draws_list){
        double proba_tox;
        draw.second.proba_tox(trial_data.dose1T[d1], trial_data.dose2T[d2],
                              proba_tox);
        count_tox += proba_tox > TARG_SUP;
        pi_median.push_back(proba_tox);

        double resp_rate =
          draw.first.responseRate(trial_data.dose1E[d1], dose2E_tau);
        count_eff += resp_rate > EFF_MIN;
        q_median.push_back(resp_rate);
      }

      result.ptox_sup[d1][d2] = double(count_tox)/n_samp;
      result.qeff_min[d1][d2] = double(count_eff)/n_samp;
      result.pi[d1][d2] = median(pi_median);
      result.resp[d1][d2] = median(q_median);
    }
  }

  for(int d1 = 0; d1 < ndose1; d1++)
    for(int d2 = 0; d2 < ndose2; d2++)
      if(d2 > result.prog_params.tau)
        result.resp[d1][d2] = result.resp[d1][d2-1];

  return result;
}

void take_if_better(const estimations& estim,
                    int& nextdose1, int& nextdose2,
                    int candidate_dose1, int candidate_dose2) {
  const int ndose1 = estim.pi.size(), ndose2 = estim.pi[0].size();

  if(nextdose1 == -1 && nextdose2 == -1) {
    nextdose1 = candidate_dose1;
    nextdose2 = candidate_dose2;
    return;
  }

  if(nextdose1 < 0 || nextdose2 < 0 ||
     candidate_dose1 < 0 || candidate_dose2 < 0 ||
     nextdose1 >= ndose1 || nextdose2 >= ndose2 ||
     candidate_dose1 >= ndose1 || candidate_dose2 >= ndose2)
    throw std::logic_error("Internal error: invalid nextdose or candidate_dose.");

  double candidateResp = estim.resp[candidate_dose1][candidate_dose2];
  double bestResp = estim.resp[nextdose1][nextdose2];

  if(candidateResp > bestResp ||
     (candidateResp == bestResp &&
      estim.pi[nextdose1][nextdose2] > estim.pi[candidate_dose1][candidate_dose2])) {
    nextdose1 = candidate_dose1;
    nextdose2 = candidate_dose2;
  }
}

pair<int, int> find_next_dose(trial_data& trial_data, double c_tox, double c_eff,
                              estimations* estimRet = NULL){
  const int ndose1 = trial_data.dose1E.size(), ndose2 = trial_data.dose2E.size();

  if(trial_data.startup_end == -1 && trial_data.stage == FINAL)
    trial_data.startup_end = trial_data.pat_incl;

  // Warning : this code should detect whether the startup has ended, even
  // if it was a long time ago.
  if(trial_data.startup_end == -1) {
    if(trial_data.pat_incl == 0)
      return make_pair(0, 0);

    vector<vector<int>> npat(ndose1, vector<int>(ndose2, 0));
    vector<vector<int>> ntox(ndose1, vector<int>(ndose2, 0));
    for(int pat = 0; pat < trial_data.pat_incl; pat++) {
      npat[trial_data.dose_adm1[pat]][trial_data.dose_adm2[pat]]++;
      if(trial_data.toxicity[pat])
        ntox[trial_data.dose_adm1[pat]][trial_data.dose_adm2[pat]]++;
    }

    vector<vector<bool>> admissible(ndose1, vector<bool>(ndose2, true));
    for(int d1d2 = 0; d1d2 <= ndose1+ndose2-2; d1d2++)
      for(int d1 = max(0, d1d2-(ndose2-1)); d1 <= min(ndose1-1, d1d2); d1++) {
        int d2 = d1d2 - d1;

        if(3*ntox[d1][d2] >= 2*COHORT_START)
          for(int dd1 = d1; dd1 < ndose1; dd1++)
            for(int dd2 = d2; dd2 < ndose2; dd2++)
              admissible[dd1][dd2] = false;

        if(!admissible[d1][d2])
          continue;

        if(npat[d1][d2] == 0 ||
           (npat[d1][d2] <= COHORT_START && 3*ntox[d1][d2] >= COHORT_START))
          return make_pair(d1, d2);
      }

    trial_data.startup_end = trial_data.pat_incl;
  }

  estimations estim = estimate(trial_data);
  if(estimRet != NULL)
    *estimRet = estim;

  int nextdose1 = -1, nextdose2 = -1;
  const int cdose1 = trial_data.cdose1;
  const int cdose2 = trial_data.cdose2;

  bool final = trial_data.stage == FINAL;
  for(int candidate_dose1 = 0; candidate_dose1 < ndose1; candidate_dose1++) {
    for(int candidate_dose2 = 0; candidate_dose2 < ndose2; candidate_dose2++) {
      bool candidate = false;
      if(!final && candidate_dose1 <= cdose1 && candidate_dose2 <= cdose2)
        candidate = true;
      if(!final && candidate_dose1 == cdose1+1 && candidate_dose2 <= cdose2)
        candidate = true;
      if(!final && candidate_dose1 <= cdose1 && candidate_dose2 == cdose2+1)
        candidate = true;
      for(int pat = 0; !candidate && pat < trial_data.pat_incl; pat++)
        if(trial_data.dose_adm1[pat] == candidate_dose1 &&
           trial_data.dose_adm2[pat] == candidate_dose2)
          candidate = true;
      if(candidate &&
         estim.ptox_sup[candidate_dose1][candidate_dose2] < c_tox &&
         (estim.qeff_min[candidate_dose1][candidate_dose2] >= c_eff ||
          trial_data.stage == FSTHALF)) {
        take_if_better(estim, nextdose1, nextdose2, candidate_dose1, candidate_dose2);
      }
    }
  }

  return make_pair(nextdose1, nextdose2);
}

}}

using namespace dfcomb::plateau;

extern "C" {

R_NativePrimitiveArgType plateau_next_args[] =
  {INTSXP, INTSXP,
   REALSXP, REALSXP,
   REALSXP, REALSXP,
   REALSXP, REALSXP,
   REALSXP, REALSXP,
   INTSXP, INTSXP,
   REALSXP, REALSXP,

   INTSXP, INTSXP,
   LGLSXP,
   REALSXP,
   INTSXP,
   INTSXP, INTSXP,
   REALSXP, REALSXP,
   LGLSXP,
   INTSXP,

   REALSXP, REALSXP,
   REALSXP, REALSXP,
   REALSXP };
const int plateau_next_nargs = 30;

void plateau_next(int* ndose1, int* ndose2,
                  double* targ_sup, double* eff_min,
                  double* dose1TV0, double* dose2TV0,
                  double* dose1EV0, double* dose2EV0,
                  double* time_full, double* cycle,
                  int* cohort_start, int* cohort,
                  double* c_tox, double* c_eff,

                  int* cdose1, int* cdose2,
                  int* in_startup,
                  double* time_cur,
                  int* pat_incl,
                  int* dose_adm1, int* dose_adm2,
                  double* time_prog, double* time_incl,
                  int* toxicity,
                  int* stage, /* 0 = first half of trial */
                              /* 1 = second half of trial */
                              /* 2 = end of the trial */

                  double* pi, double* ptox_sup,
                  double* resp, double* qeff_min,
                  double* proba_tau) {
  try {

  CYCLE = *cycle;
  TARG_SUP = *targ_sup;
  EFF_MIN = *eff_min;
  TIMEFULL = *time_full;

  vector<double> dose1TV(dose1TV0, dose1TV0+*ndose1);
  vector<double> dose2TV(dose2TV0, dose2TV0+*ndose2);
  vector<double> dose1EV(dose1EV0, dose1EV0+*ndose1);
  vector<double> dose2EV(dose2EV0, dose2EV0+*ndose2);

  COHORT_START = *cohort_start;
  COHORT = *cohort;

  struct trial_data trial_data(dose1TV, dose2TV,
                               dose1EV, dose2EV, 0/*seed - not used */);
  trial_data.pat_incl = *pat_incl;
  trial_data.cdose1 = *cdose1;
  trial_data.cdose2 = *cdose2;
  trial_data.startup_end = *in_startup ? -1 : trial_data.pat_incl;
  trial_data.time_cur = *time_cur;
  trial_data.dose_adm1 = vector<unsigned int>(dose_adm1, dose_adm1+*pat_incl);
  trial_data.dose_adm2 = vector<unsigned int>(dose_adm2, dose_adm2+*pat_incl);
  trial_data.time_prog = vector<double>(time_prog, time_prog+*pat_incl);
  trial_data.time_incl = vector<double>(time_incl, time_incl+*pat_incl);
  trial_data.toxicity = vector<int>(toxicity, toxicity+*pat_incl);
  trial_data.stage = (trial_stage)*stage;

  estimations estim(*ndose1, *ndose2);
  std::tie(*cdose1, *cdose2) =
    find_next_dose(trial_data, *c_tox, *c_eff, &estim);
  *in_startup = trial_data.startup_end == -1;

  if(!*in_startup)
    for(int d2 = 0; d2 < *ndose2; d2++) {
      for(int d1 = 0; d1 < *ndose1; d1++) {
        pi[d1 + d2**ndose1] = estim.pi[d1][d2];
        ptox_sup[d1 + d2**ndose1] = estim.ptox_sup[d1][d2];
        resp[d1 + d2**ndose1] = estim.resp[d1][d2];
        qeff_min[d1 + d2**ndose1] = estim.qeff_min[d1][d2];
      }
      proba_tau[d2] = estim.proba_tau[d2];
    }

  }
  catch (std::logic_error &e) {
    error("Internal error in dfcomb (details: %s)", e.what());
  }
  catch (...) {
    error("Internal error in dfcomb");
  }
}

R_NativePrimitiveArgType plateau_sim_args[] =
  {INTSXP, INTSXP,
   REALSXP, REALSXP,
   REALSXP, REALSXP,
   REALSXP, REALSXP,
   REALSXP, REALSXP,
   REALSXP, INTSXP,
   REALSXP, REALSXP,
   INTSXP, INTSXP, INTSXP,
   REALSXP, REALSXP,
   INTSXP, INTSXP,

   INTSXP, INTSXP, INTSXP,
   INTSXP, INTSXP,
   REALSXP};
const int plateau_sim_nargs = 27;

void plateau_sim(int* ndose1, int* ndose2,
                 double* piV, double* respV,
                 double* targ_sup, double* eff_min,
                 double* dose1TV0, double* dose2TV0,
                 double* dose1EV0, double* dose2EV0,
                 double* incl_per_week, int* npatients,
                 double* time_full, double* cycle,
                 int* cohort_start, int* cohort, int* ntrial,
                 double* c_tox, double* c_eff,
                 int* seed, int* nthreads,

                 int* inconc, int* n_pat_dose, int* rec_dose,
                 int* n_tox, int* n_eff,
                 double* duration) {
  string errstr;

  {
    struct true_data true_data;
    CYCLE = *cycle;
    TIMEFULL = *time_full;

    true_data.piV = vector<vector<double>>(*ndose1, vector<double>(*ndose2));
    true_data.respV = vector<vector<double>>(*ndose1, vector<double>(*ndose2));
    for(int d1 = 0; d1 < *ndose1; d1++)
      for(int d2 = 0; d2 < *ndose2; d2++) {
        true_data.piV[d1][d2] = piV[d1 + d2**ndose1];
        true_data.respV[d1][d2] = respV[d1 + d2**ndose1];
      }

    TARG_SUP = *targ_sup;
    EFF_MIN = *eff_min;

    vector<double> dose1TV(dose1TV0, dose1TV0+*ndose1);
    vector<double> dose2TV(dose2TV0, dose2TV0+*ndose2);
    vector<double> dose1EV(dose1EV0, dose1EV0+*ndose1);
    vector<double> dose2EV(dose2EV0, dose2EV0+*ndose2);
    true_data.incl_per_week = *incl_per_week;

    COHORT_START = *cohort_start;
    COHORT = *cohort;

    std::mt19937_64 global_rng(*seed);
    vector<uint_fast64_t> seeds;
    for(int trial=0; trial<*ntrial; trial++)
      seeds.push_back(global_rng());

    struct results results(*ndose1, *ndose2);

    #ifdef _OPENMP
    if(*nthreads > 0)
      omp_set_num_threads(*nthreads);
    else
      omp_set_num_threads(omp_get_num_procs());
#endif

    bool err = false;
    Progress prog(*ntrial);
#pragma omp parallel for schedule(dynamic, 1)
    for(int trial=0; trial<*ntrial; trial++) {
      try {
        bool err2;
#pragma omp critical
        err2 = err;
        if(prog.check_abort() || err2)
          continue;

        struct trial_data trial_data(dose1TV, dose2TV, dose1EV, dose2EV, seeds[trial]);
        while(true) {
#pragma omp critical
          err2 = err;
          if(prog.check_abort() || err2)
            goto aborted;

          if(wait_patient(trial_data, true_data))
            break;

          if(2*trial_data.pat_incl <= *npatients)
            trial_data.stage = FSTHALF;
          else
            trial_data.stage = SNDHALF;

          int startup_end = trial_data.startup_end;
          if( (startup_end == -1 && trial_data.pat_incl % COHORT_START == 0) ||
              (startup_end != -1 && (trial_data.pat_incl - startup_end) % COHORT == 0))
            std::tie(trial_data.cdose1, trial_data.cdose2) =
              find_next_dose(trial_data, *c_tox, *c_eff);

          if(trial_data.cdose1 >= 0) {
            include_patient(trial_data, true_data);
            if(trial_data.pat_incl >= *npatients)
              /* No more patiens for this group */
              trial_data.cdose1 = trial_data.cdose2 = -2;
          }
        }

        trial_data.stage = FINAL;

        pair<int, int> recom;
        if(trial_data.cdose1 == -2)
          recom = find_next_dose(trial_data, *c_tox, *c_eff);
        else
          recom = make_pair(-1, -1);

#pragma omp critical
        results.accumulate(trial_data, recom);

      } catch (std::logic_error &e) {
#pragma omp critical
        { err = true;
          errstr = e.what(); }
      } catch(...) {
#pragma omp critical
        err = true;
      }

#pragma omp critical
      prog.increment();
      aborted: ;
    }

    if(err) goto errlbl;

    *inconc = results.inconc;
    for(int d1 = 0; d1 < *ndose1; d1++)
      for(int d2 = 0; d2 < *ndose2; d2++) {
        n_pat_dose[d1 + d2**ndose1] = results.nptsdose[d1][d2];
        rec_dose[d1 + d2**ndose1] = results.nrecdose[d1][d2];
        n_tox[d1 + d2**ndose1] = results.ntox[d1][d2];
        n_eff[d1 + d2**ndose1] = results.neff[d1][d2];
      }
    *ntrial  = results.n_trials;
    *duration  = results.duration;
  }

  if(false) {
  errlbl:
    error("Internal error in dfcomb (details: %s)", errstr.c_str());
  }
}

}
