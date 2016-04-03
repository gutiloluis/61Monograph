#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>

int i;
int n_cells = 10000;
int t_total;

double p_unrep =1200.0;
double k_r_max;
double Kd = 800.0;
double n = 2.0;

double tau_r = 120.0;
double tau_p = 3600.0;
double g_r;
double g_p;
double b = 10.0;
double k_p;

void step_rep(double *t, double *r, double *p, gsl_rng *R);
void step_unrep(double *t, double *r, double *p, gsl_rng *R);
void cell(double t_total, double *r_f, double *p_f, gsl_rng *R);
void cells(double n_cells, double t_total, gsl_rng *R, char *filename);
int main(int argc, char **argv){
  
  *&t_total = 10.0*tau_p;

  *&g_r = log(2)/tau_r;
  *&g_p = log(2)/tau_p;
  *&k_p = b*g_r;

  *&k_r_max = p_unrep*g_r*g_p/k_p;

  const gsl_rng_type *T;
  gsl_rng *R;

  T = gsl_rng_mt19937;
  R = gsl_rng_alloc(T);
  gsl_rng_set(R,atoi(argv[1]));
  
  cells(n_cells, t_total, R, "h_rep.dat");

  return 0;
}

void step_unrep(double *t, double *r, double *p, gsl_rng *R){
  

  double k_r = k_r_max;
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

void step_rep(double *t, double *r, double *p, gsl_rng *R){
  
  double k_r = k_r_max/(1.0+pow(*p/Kd,n));
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
    step_rep(&t, &r, &p, R);
  }
  *r_f = r;
  *p_f = p;
}

void cells(double n_cells, double t_total, gsl_rng *R, char *filename){
  int i;

  double r_f;
  double p_f;

  FILE *out = fopen(filename, "w");
  
  for(i=0;i<n_cells;i++){
    cell(t_total, &r_f, &p_f, R);
    fprintf(out, "%f\n", p_f);
  }
  
  fclose(out);
}
