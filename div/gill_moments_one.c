#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>

#define k_r 0.01 /* 1/s */
#define tau_r 120.0  /* s */
#define tau_p 3600.0 /* s */
#define g_r log(2)/tau_r
#define g_p log(2)/tau_p
#define b 20.0
#define k_p b*g_r

void step(double *t, double *r, double *p, gsl_rng *R);
void cell(double t_total, double *r_f, double *p_f, gsl_rng *R);
void cells(double n_cells, double t_total, gsl_rng *R);
int main(int argc, char **argv){
  
  /*
  double k_r = 0.01;
  double tau_r = 120.0;
  double tau_p = 3600.0;
  double g_r = log(2)/tau_r;
  double g_p = log(2)/tau_p;
  double b = 20.0;
  double k_p = b*g_r;
  */
  
  const gsl_rng_type *T;
  gsl_rng *R;

  T = gsl_rng_mt19937;
  R = gsl_rng_alloc(T);
  gsl_rng_set(R,atoi(argv[1]));

  cells(100, 10.0*tau_p, R);
  
  return 0;
}

void step(double *t, double *r, double *p, gsl_rng *R){
  
  double K = k_r + g_r**r + k_p**r + g_p**p;
  
  double Dt = -log(gsl_rng_uniform(R))/K;
  double which = gsl_rng_uniform(R);

  *t+=Dt;

  if(which < k_r/K){
    *r+=1;
  } else if(which >= k_r/K && which < (k_r+g_r**r)/K){
    *r-=1;
  } else if(which >= (k_r+g_r**r)/K && which < (k_r+g_r**r+k_p**r)/K){
    *p+=1;
  } else {
    *p-=1;
  }
}

void cell(double t_total, double *r_f, double *p_f, gsl_rng *R){
  
  double t = 0;
  double r = 0;
  double p = 0;
  
  while(t<t_total){
    step(&t, &r, &p, R);
  }
  *r_f = r;
  *p_f = p;
}

void cells(double n_cells, double t_total, gsl_rng *R){
  int i;
  double ave_r = 0;
  double ave_r2 = 0;

  double r_f;
  double p_f;

  for(i=0;i<n_cells;i++){
    cell(t_total, &r_f, &p_f, R);
    ave_r += r_f;
    ave_r2 += r_f*r_f;
  }
  
  ave_r/=n_cells;
  ave_r2/=n_cells;

  printf("%f %f\n", ave_r, ave_r2 - ave_r*ave_r);
}
