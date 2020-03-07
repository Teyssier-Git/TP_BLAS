#include <stdio.h>
#include <x86intrin.h>

#include "mnblas.h"
#include "complexe.h"

#include "flop.h"

#define    NB_FOIS        4194304

int main (int argc, char **argv)
{
  complexe_float_t c1= {1.0, 2.0} ;
  complexe_float_t c2= {3.0, 6.0} ;
  complexe_float_t c3 ;
  complexe_float_t c3m ;
  complexe_float_t c3d ;

  complexe_double_t cd1 ;
  complexe_double_t cd2 ;
  complexe_double_t cd3 ;
  complexe_double_t cd3m ;
  complexe_double_t cd3d ;

  unsigned long long int start, end ;
  int i ;

  init_flop () ;

  c3 = add_complexe_float (c1, c2) ;

  printf ("c3.r %f c3.i %f\n", c3.real, c3.imaginary) ;

  cd1 = (complexe_double_t) {10.0, 7.0} ;
  cd2 = (complexe_double_t) {25.0, 32.0} ;

  cd3 = add_complexe_double (cd1, cd2) ;

  printf ("cd3.r %f cd3.i %f\n", cd3.real, cd3.imaginary) ;

  start =_rdtsc () ;

  for (i = 0 ; i < NB_FOIS; i++)
  {
    cd3 = add_complexe_double (cd1, cd3) ;
  }

  end = _rdtsc () ;

  printf ("apres boucle cd3.real %f cd3.imaginary %f %lld cycles \n", cd3.real, cd3.imaginary, end-start) ;

  calcul_flop ("calcul complexe ", NB_FOIS*2, end-start) ;
  c3m = mult_complexe_float (c1, c2) ;

  printf ("c3m.r %f c3m.i %f\n", c3m.real, c3m.imaginary) ;
  start =_rdtsc () ;

  for (i = 0 ; i < NB_FOIS; i++)
  {
    cd3m = mult_complexe_double (cd1, cd3m) ;
  }

  end = _rdtsc () ;

  printf ("apres boucle cd3m.real %f cd3m.imaginary %f %lld cycles \n", cd3m.real, cd3m.imaginary, end-start) ;

  calcul_flop ("calcul complexe ", NB_FOIS*6, end-start) ;

  c3d = div_complexe_float (c1, c2) ;

  printf ("c3d.r %f c3d.i %f\n", c3d.real, c3d.imaginary) ;
  start =_rdtsc () ;

  for (i = 0 ; i < NB_FOIS; i++)
  {
    cd3d = div_complexe_double (cd3d, cd2) ;
  }

  end = _rdtsc () ;

  printf ("apres boucle cd3d.real %f cd3d.imaginary %f %lld cycles \n", cd3d.real, cd3d.imaginary, end-start) ;

  calcul_flop ("calcul complexe ", NB_FOIS*11, end-start) ;
  exit (0) ;
}
