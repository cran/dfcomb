#include "dfcomb.h"

static R_CMethodDef cMethods[] = {
  {"logistic_next", (DL_FUNC) &logistic_next,
   0, logistic_next_args},
  {"logistic_sim", (DL_FUNC) &logistic_sim,
   0, logistic_sim_args},
  {NULL, NULL, 0}
};

void R_init_dfcomb(DllInfo *info) {
  cMethods[0].numArgs = logistic_next_nargs;
  cMethods[1].numArgs = logistic_sim_nargs;
  R_registerRoutines(info, cMethods, NULL, NULL, NULL);
}
