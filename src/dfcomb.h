#include <R_ext/Rdynload.h>
#include <Rinternals.h>

extern R_NativePrimitiveArgType logistic_sim_args[];
extern const int logistic_sim_nargs;
void logistic_sim(int* tite, int* ndose1, int* ndose2,
                  double* timefull, double* week_incl,
                  double* piv,
                  double* target, double* target_max, double* target_min,
                  double* prior_1, double* prior_2,
                  int* ncohort, int* cohort, int* ntrial,
                  double* escp, double* desp, double* arret,
                  int* nmin, int* seed,

                  double* nrecdoserat, double* nptsdoserat, double* ntoxrat,
                  double* inconcrat);

extern R_NativePrimitiveArgType logistic_next_args[];
extern const int logistic_next_nargs;
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
                   double* ptox_targ, double* ptox_sup_targ);

extern R_NativePrimitiveArgType plateau_next_args[];
extern const int plateau_next_nargs;
void plateau_next(int* ndose1, int* ndose2,
                  double* targ_sup, double* eff_min,
                  double* dose1TV0, double* dose2TV0,
                  double* dose1EV0, double* dose2EV0,
                  double* time_full, double* cycle,
                  int* cohort_start, int* cohort,
                  double* c_tox, double* c_eff,

                  int* in_startup,
                  int* cdose1, int* cdose2,
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
                  double* proba_tau);

extern R_NativePrimitiveArgType plateau_sim_args[];
extern const int plateau_sim_nargs;
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
                 int* n_tox, int* n_prog,
                 double* duration);

