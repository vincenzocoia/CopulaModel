#include <R.h>
#include <Rmath.h>

/* unif_rand() is the RNG in R */
/* This is used in several r*.c files */
double urand() { return(unif_rand()); }

