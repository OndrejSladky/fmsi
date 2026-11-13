/* Compiled by openmp.mk to find out whether the C compiler supports OpenMP. */
#include <omp.h>

int main(void) { return omp_get_max_threads() < 1; }
