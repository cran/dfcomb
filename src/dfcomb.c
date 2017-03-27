#include "dfcomb.h"

static R_CMethodDef cMethods[] = {
  {"logistic_next", (DL_FUNC) &logistic_next,
   0, logistic_next_args},
  {"logistic_sim", (DL_FUNC) &logistic_sim,
   0, logistic_sim_args},
  {"plateau_next", (DL_FUNC) &plateau_next,
   0, plateau_next_args},
  {"plateau_sim", (DL_FUNC) &plateau_sim,
   0, plateau_sim_args},
  {NULL, NULL, 0}
};

void R_init_dfcomb(DllInfo *dll) {
  cMethods[0].numArgs = logistic_next_nargs;
  cMethods[1].numArgs = logistic_sim_nargs;
  cMethods[2].numArgs = plateau_next_nargs;
  cMethods[3].numArgs = plateau_sim_nargs;
  R_registerRoutines(dll, cMethods, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
