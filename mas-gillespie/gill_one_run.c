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
void cell(double t_total, gsl_rng *R);
int main(){
  
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
  gsl_rng_set(R,0);

  cell(10.0*tau_p, R);
  
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

void cell(double t_total, gsl_rng *R){
  
  double t = 0;
  double r = 0;
  double p = 0;
  
  printf("%f %f %f\n",t,r,p);
  
  while(t<t_total){
    step(&t, &r, &p, R);
    printf("%f %f %f\n",t,r,p);
  }
}
