#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

double a;
void f(double *x, double *y, gsl_rng *r);
int main(){
  
  const gsl_rng_type *T;
  gsl_rng * r;
  T = gsl_rng_mt19937;
  r = gsl_rng_alloc(T);
  /*gsl_rng_set(r, 0);*/
  
  FILE *out = fopen("a.dat", "a");
  fprintf(out, "Hola\n");
  fclose(out);
 
  double x = 0;
  double y = 0;
  int i;
  for(i=0;i<10;i++){
    f(&x,&y, r);
    /* printf("%f\n",y);*/
  }
  return 0;
}

void f(double *x, double *y, gsl_rng *r){
  /* *x+=gsl_rng_uniform(r);*/
  *y=gsl_ran_gaussian(r, 1);
  FILE *out = fopen("a.dat", "a");
  fprintf(out, "%f\n",1.0);
  fclose(out);
}
