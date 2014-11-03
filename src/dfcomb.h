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
                   int* cdose1, int* cdose2, int* in_startup,
                   int* dose_adm1, int* dose_adm2,
                   double* time_ev, double* time_follow /* Only for tite */,
                   int* delta /* Only for non-tite */,

                   int* inconc,
                   double* pi, double* ptox, double* ptox_inf_targ,
                   double* ptox_targ, double* ptox_sup_targ);
