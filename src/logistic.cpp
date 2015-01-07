#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <algorithm>
#include <random>

#define BOOST_DISABLE_ASSERTS

#include <boost/random/exponential_distribution.hpp>

extern "C" {
  #include "arms.h"
  #include "dfcomb.h"
}

#include <R.h>
#include <Rcpp.h>

#include <progress.hpp>

using namespace std;

namespace dfcomb { namespace logistic {

struct datastru{
  double d0, a0, b0, c0;
  vector<double> dose_scale_1, dose_scale_2;
  vector<vector<double> > pi, ptox, ptox_inf_targ, ptox_targ, ptox_sup_targ;
  vector<vector<int> > y, n;
  vector<bool> delta;
  vector<int> dose_adm1, dose_adm2;
  vector<double> time_ev, time_follow, time_min;
  int pat_incl;
  int cdose1, cdose2;

  datastru(int NDOSE1, int NDOSE2) :
    dose_scale_1(NDOSE1), dose_scale_2(NDOSE2),
    pi(NDOSE1, vector<double>(NDOSE2)), ptox(NDOSE1, vector<double>(NDOSE2)),
    ptox_inf_targ(NDOSE1, vector<double>(NDOSE2)), ptox_targ(NDOSE1, vector<double>(NDOSE2)), ptox_sup_targ(NDOSE1, vector<double>(NDOSE2)),
    y(NDOSE1, vector<int>(NDOSE2)), n(NDOSE1, vector<int>(NDOSE2))
  { }
};

enum para_ty {
  PARA_d0,
  PARA_a0,
  PARA_b0,
  PARA_c0
} PARA;
int NDOSE1, NDOSE2;
double TARGET, TARGET_MIN, TARGET_MAX, WEEK_incl, TIMEFULL;
int NCOHORT, COHORT;
int TITE;
double ESCP, DESP, ARRET, NMIN;

std::minstd_rand r;

// Cohort inclusion
void genpopn(datastru * datapt, const vector<vector<double> >& piV){
  vector<double> time_incl_cohort(COHORT);
  double pi_ij = piV[datapt->cdose1][datapt->cdose2];

  boost::exponential_distribution<double> exp_rng;
  std::uniform_real_distribution<double> uni_rng;
  for(int i=0; i<COHORT; i++){
    datapt->dose_adm1.push_back(datapt->cdose1);
    datapt->dose_adm2.push_back(datapt->cdose2);
    if(TITE) {
      datapt->delta.push_back(false);
      double time_btw_incl = exp_rng(r, WEEK_incl);
      time_incl_cohort[i] = time_btw_incl;
      double time_tox = exp_rng(r, -log(1-pi_ij)/TIMEFULL);
      datapt->time_ev.push_back(time_tox);
      datapt->time_follow.push_back(0);
      datapt->time_min.push_back(0);
      for(int j=i; j>=0; j--){
        datapt->time_follow[datapt->pat_incl+j] += time_incl_cohort[i];
      }
    } else {
      bool tox = uni_rng(r)<pi_ij;
      datapt->y[datapt->cdose1][datapt->cdose2] += (int)tox;
      datapt->delta.push_back(tox);
    }
  }

  if(TITE) {
    for(int j=0; j<COHORT; j++){
      for(int i=0; i<datapt->pat_incl; i++){
        datapt->time_follow[i] += time_incl_cohort[j];
      }
    }
  }

  datapt->n[datapt->cdose1][datapt->cdose2] += COHORT;
  datapt->pat_incl += COHORT;

  if(TITE) {
    if(datapt->pat_incl == COHORT*NCOHORT){
      for(int i=0; i<datapt->pat_incl; i++){
        datapt->time_follow[i] = INFINITY;
      }
    }

    for(int i=0; i<NDOSE1; i++){
      for(int j=0; j<NDOSE2; j++){
        datapt->y[i][j] = 0;
      }
    }

    for(int i=0; i<datapt->pat_incl; i++){
      datapt->delta[i] = datapt->time_ev[i] <= min(datapt->time_follow[i], TIMEFULL);
      datapt->time_min[i] = min(datapt->time_ev[i], min(datapt->time_follow[i],TIMEFULL));
      datapt->y[datapt->dose_adm1[i]][datapt->dose_adm2[i]] += (int)datapt->delta[i];
    }
  }
}


// Start-up phase
void startup(datastru *datapt, const vector<vector<double> >& piV){
  while(true){
    genpopn(datapt, piV);

    if(datapt->cdose1==NDOSE1-1 && datapt->cdose2==NDOSE2-1){
      break;
    }

    if(datapt->y[datapt->cdose1][datapt->cdose2] != 0){
      break;
    }

    if(datapt->cdose1<NDOSE1-1){
      datapt->cdose1++;
    }
    if(datapt->cdose2<NDOSE2-1){
      datapt->cdose2++;
    }
  }
}

// Toxicity probability model
double proba_tox(double x_p, double y_q, double d0, double a0, double b0, double c0){
  double x = d0+a0*x_p+b0*y_q+c0*x_p*y_q;
  return( exp(x)/(1+exp(x)) );
}

double density(double x,  void *datapt){
  struct datastru *dp;
  dp = (datastru *) datapt;
  double d0, a0, b0, c0;
  d0=dp->d0; a0=dp->a0; b0=dp->b0; c0=dp->c0;

  double logprior;
  switch(PARA) {
    case PARA_d0:
      d0=x;
      logprior=-0.05*d0*d0;
      break;
    case PARA_a0:
      a0=x;
      logprior=-a0;
      break;
    case PARA_b0:
      b0=x;
      logprior=-b0;
      break;
    case PARA_c0:
      c0=x;
      logprior=-0.05*c0*c0;
      break;
    default:
      throw std::logic_error("Internal error: invalid PARA.");
  }

  double pi_ij;
  double loglike=0;

  if(TITE) {
    int card2=0;
    double weight;
    double p_tox;
    for(int j=0; j<dp->pat_incl; j++){
      if(dp->time_follow[j] >= TIMEFULL && dp->delta[j]){
        card2 += 1;
      }
    }
    for(int p=0; p<dp->pat_incl; p++){
      p_tox = proba_tox(dp->dose_scale_1[dp->dose_adm1[p]],
                        dp->dose_scale_2[dp->dose_adm2[p]], d0, a0, b0, c0);
      if(dp->delta[p]){
        weight = 1;
      }
      else{
        int card1=0;
        for(int j=0; j<dp->pat_incl; j++){
          if(dp->time_follow[j] >= TIMEFULL &&
             dp->time_ev[j] <= min(dp->time_follow[p],TIMEFULL)){
            card1 += 1;
          }
        }

        double w_ini = dp->time_min[p] / TIMEFULL;
        weight = (double)(card1+w_ini)/(card2+1);
        if(weight < 0 || weight > 1){
          throw std::logic_error("Internal error: invalid weight.");
        }
      }
      loglike += dp->delta[p] ? log(weight*p_tox) : log(1-weight*p_tox);
    }
  } else {
    for(int i=0; i<NDOSE1; i++){
      for(int j=0; j<NDOSE2; j++){
        if(dp->n[i][j] != 0){
          pi_ij = proba_tox(dp->dose_scale_1[i], dp->dose_scale_2[j], d0, a0, b0, c0);
          loglike += dp->y[i][j]*log(pi_ij) + (dp->n[i][j]-dp->y[i][j])*log(1-pi_ij);
        }
      }
    }
  }

  return(loglike+logprior);
}

std::minstd_rand r_estim;
std::uniform_real_distribution<double> uni_rng_estim;

double u_random() {
  return uni_rng_estim(r_estim);
}

// Estimations
void estimation(datastru * datapt){
  r_estim.seed(53425);
  uni_rng_estim.reset();

  int nburn=2000;
  int niter=5000;
  datapt->d0=1, datapt->a0=1, datapt->b0=1, datapt->c0=0;

  // paramters for ARMS
  int ninit=4, dometrop=1;
  double d0L=-8, d0R=8, d0prev=1, d0prev2=1;
  double a0L=0.01, a0R=8, a0prev=1, a0prev2=1;
  double b0L=0.01, b0R=8, b0prev=1, b0prev2=1;
  double c0L=-8, c0R=8, c0prev=0, c0prev2=0;

  for(int i=0; i<NDOSE1; i++){
    for(int j=0; j<NDOSE2; j++){
      datapt->pi[i][j]=0;
      datapt->ptox[i][j]=0;
      datapt->ptox_inf_targ[i][j]=0;
      datapt->ptox_targ[i][j]=0;
      datapt->ptox_sup_targ[i][j]=0;
    }
  }

  for(int iter=0; iter<nburn+niter; iter++){
    int err=0;
    double d0samp;
    double a0samp;
    double b0samp;
    double c0samp;

    while(true){
      PARA = PARA_d0;
      err = arms_simple(ninit, &d0L, &d0R, density, datapt, dometrop, &d0prev, &d0samp, u_random);
      if(err>0) {
        char buf[100]; sprintf(buf, "%d", err);
        throw std::logic_error(string("arms internal error (d0): ") + string(buf));
      }
      datapt->d0 = d0samp;
      d0prev2 = d0prev;
      d0prev = d0samp;

      PARA = PARA_a0;
      err = arms_simple(ninit, &a0L, &a0R, density, datapt, dometrop, &a0prev, &a0samp, u_random);
      if(err>0) {
        char buf[100]; sprintf(buf, "%d", err);
        throw std::logic_error(string("arms internal error (a0): ") + string(buf));
      }
      datapt->a0 = a0samp;
      a0prev2 = a0prev;
      a0prev =  a0samp;

      PARA = PARA_b0;
      err = arms_simple(ninit, &b0L, &b0R, density, datapt, dometrop, &b0prev, &b0samp, u_random);
      if(err>0) {
        char buf[100]; sprintf(buf, "%d", err);
        throw std::logic_error(string("arms internal error (b0): ") + string(buf));
      }
      datapt->b0 = b0samp;
      b0prev2 = b0prev;
      b0prev = b0samp;

      PARA = PARA_c0;
      err = arms_simple(ninit, &c0L, &c0R, density, datapt, dometrop, &c0prev, &c0samp, u_random);
      if(err>0) {
        char buf[100]; sprintf(buf, "%d", err);
        throw std::logic_error(string("arms internal error (c0): ") + string(buf));
      }
      datapt->c0 = c0samp;
      c0prev2 = c0prev;
      c0prev = c0samp;

      if(a0samp+c0samp*(datapt->dose_scale_2[0]) >= 0 &&
         a0samp+c0samp*(datapt->dose_scale_2[NDOSE2-1]) >= 0 &&
         b0samp+c0samp*(datapt->dose_scale_1[0]) >= 0 &&
         b0samp+c0samp*(datapt->dose_scale_1[NDOSE1-1]) >= 0)
        break;

      datapt->d0 = d0prev2;
      d0prev = d0prev2;
      datapt->a0 = a0prev2;
      a0prev = a0prev2;
      datapt->b0 = b0prev2;
      b0prev = b0prev2;
      datapt->c0 = c0prev2;
      c0prev = c0prev2;
    }

    if(iter>=nburn){
      for(int i=0; i<NDOSE1; i++){
        for(int j=0; j<NDOSE2; j++){
          double pi_ij;
          pi_ij=proba_tox(datapt->dose_scale_1[i], datapt->dose_scale_2[j], d0samp, a0samp, b0samp, c0samp);
          datapt->pi[i][j] += pi_ij;
          datapt->ptox[i][j] += (pi_ij<TARGET)?1:0;
          datapt->ptox_inf_targ[i][j] += (pi_ij<TARGET_MIN)?1:0;
          datapt->ptox_targ[i][j] += (pi_ij>=TARGET_MIN && pi_ij<=TARGET_MAX)?1:0;
          datapt->ptox_sup_targ[i][j] += (pi_ij>TARGET_MAX)?1:0;
        }
      }
    }
  }
  for(int i=0; i<NDOSE1; i++){
    for(int j=0; j<NDOSE2; j++){
      datapt->pi[i][j] /= niter;
      datapt->ptox[i][j] /= niter;
      datapt->ptox_inf_targ[i][j] /= niter;
      datapt->ptox_targ[i][j] /= niter;
      datapt->ptox_sup_targ[i][j] /= niter;
    }
  }
}


void take_if_better(datastru * datapt, int& recomdose1, int& recomdose2, int candidate_dose1, int candidate_dose2) {
  if(recomdose1 == -1 && recomdose2 == -1) {
    recomdose1 = candidate_dose1;
    recomdose2 = candidate_dose2;
    return;
  }

  double candidateTarg = fabs(datapt->pi[candidate_dose1][candidate_dose2]-TARGET);
  double bestTarg = fabs(datapt->pi[recomdose1][recomdose2]-TARGET);

  if(candidateTarg < bestTarg) {
    recomdose1 = candidate_dose1;
    recomdose2 = candidate_dose2;
  }
}


// Determination of the next combination
bool newdfrule(datastru * datapt){
  int recomdose1=-1, recomdose2=-1;
  int cdose1 = datapt->cdose1;
  int cdose2 = datapt->cdose2;

  if(datapt->ptox[cdose1][cdose2] > ESCP){
    for(int candidate_dose1 = max(cdose1-1,0); candidate_dose1 <= min(NDOSE1-1, cdose1+1); candidate_dose1++) {
      for(int candidate_dose2 = max(cdose2-1,0); candidate_dose2 <= min(NDOSE2-1, cdose2+1); candidate_dose2++) {
        if(candidate_dose1 == cdose1+1 && candidate_dose2 == cdose2+1)
          continue;
        if(datapt->pi[candidate_dose1][candidate_dose2] <= datapt->pi[cdose1][cdose2])
          continue;
        if(candidate_dose1 == cdose1+1 || candidate_dose2 == cdose2+1) {
          take_if_better(datapt, recomdose1, recomdose2, candidate_dose1, candidate_dose2);
        }
      }
    }

    if(recomdose1==-1 && recomdose2==-1){
      recomdose1 = cdose1;
      recomdose2 = cdose2;
    }
  }
  else if(datapt->ptox[cdose1][cdose2] < DESP){
    for(int candidate_dose1 = max(cdose1-1,0); candidate_dose1 <= min(NDOSE1-1, cdose1+1); candidate_dose1++) {
      for(int candidate_dose2 = max(cdose2-1,0); candidate_dose2 <= min(NDOSE2-1, cdose2+1); candidate_dose2++) {
        if(candidate_dose1 == cdose1-1 && candidate_dose2 == cdose2-1)
          continue;
        if(datapt->pi[candidate_dose1][candidate_dose2] >= datapt->pi[cdose1][cdose2])
          continue;
        if(candidate_dose1 == cdose1-1 || candidate_dose2 == cdose2-1) {
          take_if_better(datapt, recomdose1, recomdose2, candidate_dose1, candidate_dose2);
        }
      }
    }

    if(recomdose1==-1 && recomdose2==-1){
      if(datapt->n[cdose1][cdose2] > NMIN && 1-datapt->ptox[cdose1][cdose2] >= ARRET){
        return true;
      }
      else{
        recomdose1 = cdose1;
        recomdose2 = cdose2;
      }
    }

  }
  else{
    recomdose1 = cdose1;
    recomdose2 = cdose2;
  }

  datapt->cdose1 = recomdose1;
  datapt->cdose2 = recomdose2;
  return false;
}

}}

using namespace dfcomb::logistic;

R_NativePrimitiveArgType logistic_sim_args[] =
  {LGLSXP, INTSXP, INTSXP,
   REALSXP, REALSXP,
   REALSXP,
   REALSXP, REALSXP, REALSXP,
   REALSXP, REALSXP,
   INTSXP, INTSXP, INTSXP,
   REALSXP, REALSXP, REALSXP,
   INTSXP, INTSXP,

   REALSXP, REALSXP, REALSXP,
   REALSXP};
const int logistic_sim_nargs = 23;

void logistic_sim(int* tite, int* ndose1, int* ndose2,
                  double* timefull, double* week_incl,
                  double* piv,
                  double* target, double* target_max, double* target_min,
                  double* prior_1, double* prior_2,
                  int* ncohort, int* cohort, int* ntrial,
                  double* escp, double* desp, double* arret,
                  int* nmin, int* seed,

                  double* nrecdoserat, double* nptsdoserat, double* ntoxrat,
                  double* inconcrat)
{
  try {

  ESCP = *escp;
  DESP = *desp;
  ARRET = *arret;
  NMIN = *nmin;
  r.seed(*seed);

  TITE = *tite;
  NDOSE1 = *ndose1;
  NDOSE2 = *ndose2;

  struct datastru data(NDOSE1, NDOSE2);

  vector<vector<double> > piV(NDOSE1, vector<double>(NDOSE2));

  for(int i=0; i<NDOSE1; i++){
    for(int j=0; j<NDOSE2; j++){
      piV[i][j] = piv[i + j*NDOSE1];
    }
  }

  if(TITE) {
    TIMEFULL = *timefull;
    WEEK_incl = *week_incl;
  }

  TARGET = *target;
  TARGET_MIN = *target_min;
  TARGET_MAX = *target_max;

  for(int i=0; i<NDOSE1; i++){
    data.dose_scale_1[i] = log(prior_1[i]/(1-prior_1[i]));
  }

  for(int i=0; i<NDOSE2; i++){
    data.dose_scale_2[i] = log(prior_2[i]/(1-prior_2[i]));
  }

  NCOHORT = *ncohort;
  COHORT = *cohort;

  int inconc=0;
  int ntreated=0;
  vector<vector<int> > nptsdose(NDOSE1, vector<int>(NDOSE2, 0)),
    ntox(NDOSE1, vector<int>(NDOSE2, 0)), nrecdose(NDOSE1, vector<int>(NDOSE2, 0));

  // Trials simulations
  Progress prog(*ntrial);
  for(int trial=0; trial<*ntrial; trial++) {
    data.pat_incl = 0;
    for(int i=0; i<NDOSE1; i++){
      for(int j=0; j<NDOSE2; j++){
        data.y[i][j]=0;
        data.n[i][j]=0;
      }
    }

    data.cdose1 = 0;
    data.cdose2 = 0;

    // Start-up phase
    startup(&data, piV);

    // Cohort inclusions and estimations, based on the model
    while(data.pat_incl < COHORT * NCOHORT){
      if(prog.check_abort()) {
        *ntrial = trial;
        goto aborted;
      }

      estimation(&data);

      if(newdfrule(&data)) {
        inconc++;
        goto trial_end;
      }

      genpopn(&data, piV);
    }

    if(TITE) {
      for(int i=0; i<data.pat_incl; i++){
        data.time_follow[i] = INFINITY;
      }

      for(int i=0; i<NDOSE1; i++){
        for(int j=0; j<NDOSE2; j++){
          data.y[i][j] = 0;
        }
      }

      for(int i=0; i<data.pat_incl; i++){
        data.delta[i] = data.time_ev[i] <= TIMEFULL;
        data.time_min[i] = min(data.time_ev[i], TIMEFULL);
        data.y[data.dose_adm1[i]][data.dose_adm2[i]] += (int)data.delta[i];
      }
    }

    estimation(&data);

    if(newdfrule(&data)) {
      inconc++;
      goto trial_end;
    }

    { int recom1 = -1, recom2 = -1;
      double closest=0.0;
      for(int i=0; i<NDOSE1; i++){
        for(int j=0; j<NDOSE2; j++){
          if(data.n[i][j]!=0){
            double dis = data.ptox_targ[i][j];
            if(dis >= closest){
              closest = dis;
              recom1=i;
              recom2=j;
            }
          }
        }
      }
      if(recom1 == -1 || recom2 == -1)
        throw std::logic_error("Internal error: no recommended dose.");
      nrecdose[recom1][recom2]++;
    }

    trial_end:

    ntreated += data.pat_incl;
    for(int i=0; i<NDOSE1; i++){
      for(int j=0; j<NDOSE2; j++){
        ntox[i][j] += data.y[i][j];
        nptsdose[i][j] += data.n[i][j];
      }
    }

    prog.increment();
  }

  aborted:

  for(int i=0; i<NDOSE1; i++){
    for(int j=0; j<NDOSE2; j++){
      nrecdoserat[i + j*NDOSE1] = (double)nrecdose[i][j] / *ntrial;
      nptsdoserat[i + j*NDOSE1] = (double)nptsdose[i][j] / *ntrial;
      ntoxrat[i + j*NDOSE1] = (double)ntox[i][j] / *ntrial;
    }
  }

  *inconcrat = (double)inconc / *ntrial;

  }
  catch (std::logic_error &e) {
    error("Internal error in dfcomb (details: %s)", e.what());
  }
  catch (...) {
    error("Internal error in dfcomb");
  }

  return;
}

R_NativePrimitiveArgType logistic_next_args[] =
  {LGLSXP, INTSXP, INTSXP,
   REALSXP,
   REALSXP, REALSXP, REALSXP,
   REALSXP, REALSXP,
   LGLSXP,
   REALSXP, REALSXP, REALSXP,
   INTSXP,

   INTSXP,
   INTSXP, INTSXP,
   INTSXP, INTSXP,
   REALSXP, REALSXP,
   LGLSXP,

   LGLSXP, LGLSXP,
   REALSXP, REALSXP, REALSXP,
   REALSXP, REALSXP };
const int logistic_next_nargs = 29;

void logistic_next(int* tite, int* ndose1, int* ndose2,
                   double* timefull,
                   double* target, double* target_max, double* target_min,
                   double* prior_1, double* prior_2,
                   int* trial_end,
                   double* escp, double* desp, double* arret,
                   int* nmin,

                   int* pat_incl,
                   int* cdose1, int* cdose2,
                   int* dose_adm1, int* dose_adm2,
                   double* time_ev, double* time_follow /* Only for tite */,
                   int* delta /* Only for non-tite */,

                   int* in_startup, int* inconc,
                   double* pi, double* ptox, double* ptox_inf_targ,
                   double* ptox_targ, double* ptox_sup_targ)
{
  try {

  ESCP = *escp;
  DESP = *desp;
  ARRET = *arret;
  NMIN = *nmin;

  TITE = *tite;
  NDOSE1 = *ndose1;
  NDOSE2 = *ndose2;

  struct datastru data(NDOSE1, NDOSE2);

  if(TITE) {
    TIMEFULL = *timefull;
  }

  TARGET = *target;
  TARGET_MIN = *target_min;
  TARGET_MAX = *target_max;

  for(int i=0; i<NDOSE1; i++){
    data.dose_scale_1[i] = log(prior_1[i]/(1-prior_1[i]));
  }

  for(int i=0; i<NDOSE2; i++){
    data.dose_scale_2[i] = log(prior_2[i]/(1-prior_2[i]));
  }

  COHORT = 1;
  NCOHORT = *trial_end ? *pat_incl : *pat_incl + 1;

  data.pat_incl = *pat_incl;
  data.cdose1 = *cdose1;
  data.cdose2 = *cdose2;
  for(int i=0; i < *pat_incl; i++) {
    data.dose_adm1.push_back(dose_adm1[i]);
    data.dose_adm2.push_back(dose_adm2[i]);
    data.n[dose_adm1[i]][dose_adm2[i]]++;
    if(TITE) {
      data.time_ev.push_back(time_ev[i]);
      if(*trial_end && time_follow[i] < TIMEFULL)
        error("dfcomb : the final recommendation cannot be computed when "
              "all the patients have not been fully followed");
      data.time_follow.push_back(time_follow[i]);
      data.time_min.push_back(min(time_ev[i], min(time_follow[i], TIMEFULL)));
      data.delta.push_back(data.time_ev[i] <= min(time_follow[i], TIMEFULL));
    } else {
      data.delta.push_back(delta[i]);
    }
    data.y[dose_adm1[i]][dose_adm2[i]] += (int)data.delta[i];
  }

  // Startup
  if(!*trial_end)
    for(int d = 0; d < NDOSE1 && d < NDOSE2; d++) {
      if(data.y[d][d] > 0)
        break;
      if(data.n[d][d] == 0) {
        *cdose1 = *cdose2 = d;
        *inconc = 0;
        return;
      }
    }

  estimation(&data);
  *inconc = newdfrule(&data);

  if(!*inconc && *trial_end) {
    int recom1 = -1, recom2 = -1;
    double closest=0.0;
    for(int i=0; i<NDOSE1; i++){
      for(int j=0; j<NDOSE2; j++){
        if(data.n[i][j]!=0){
          double dis = data.ptox_targ[i][j];
          if(dis >= closest){
            closest = dis;
            recom1=i;
            recom2=j;
          }
        }
      }
    }
    if(recom1 == -1 || recom2 == -1)
      throw std::logic_error("Internal error: no recommended dose.");
    data.cdose1 = recom1;
    data.cdose2 = recom2;
  }

  *cdose1 = data.cdose1;
  *cdose2 = data.cdose2;

  for(int i=0; i<NDOSE1; i++){
    for(int j=0; j<NDOSE2; j++){
      pi[i + j*NDOSE1] = data.pi[i][j];
      ptox[i + j*NDOSE1] = data.ptox[i][j];
      ptox_inf_targ[i + j*NDOSE1] = data.ptox_inf_targ[i][j];
      ptox_targ[i + j*NDOSE1] = data.ptox_targ[i][j];
      ptox_sup_targ[i + j*NDOSE1] = data.ptox_sup_targ[i][j];
    }
  }

  }
  catch (std::logic_error &e) {
    error("Internal error in dfcomb (details: %s)", e.what());
  }
  catch (...) {
    error("Internal error in dfcomb");
  }
}
