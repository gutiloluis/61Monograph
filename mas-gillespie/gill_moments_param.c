#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>

int i;
int n_cells = 10000;
int n_elements = 8;
int t_total;

double k_r = 0.01;
double tau_r = 120.0;
double tau_p = 3600.0;
double g_r;
double g_p;
double b = 20.0;
double k_p;

void step(double *t, double *r, double *p, gsl_rng *R);
void cell(double t_total, double *r_f, double *p_f, gsl_rng *R);
void cells(double n_cells, double t_total, gsl_rng *R, char *filename);
void param_B(gsl_rng *R);
void param_K_r(gsl_rng *R);
void param_Tau_p(gsl_rng *R);
int main(int argc, char **argv){
  
  *&t_total = 10.0*tau_p;

  *&g_r = log(2)/tau_r;
  *&g_p = log(2)/tau_p;
  *&k_p = b*g_r;
  
  const gsl_rng_type *T;
  gsl_rng *R;

  T = gsl_rng_mt19937;
  R = gsl_rng_alloc(T);
  gsl_rng_set(R,atoi(argv[1]));
  
  param_B(R);
  param_K_r(R);
  param_Tau_p(R);

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

void cells(double n_cells, double t_total, gsl_rng *R, char *filename){
  int i;
  double ave_p = 0;
  double ave_p2 = 0;

  double r_f;
  double p_f;

  FILE *out = fopen(filename, "a");
  
  for(i=0;i<n_cells;i++){
    cell(t_total, &r_f, &p_f, R);
    ave_p += p_f;
    ave_p2 += p_f*p_f;
  }
  
  ave_p/=n_cells;
  ave_p2/=n_cells;
  
  fprintf(out, "%f %f\n", ave_p, (ave_p2 - ave_p*ave_p)/ave_p);
  fclose(out);
}

void param_B(gsl_rng *R){
  
  char filename[10] = "b.dat";
  for(i=0;i<n_elements;i++){
    *&b = 5.0*(i+1);
    *&k_p = b*g_r;
    cells(n_cells, t_total, R, filename);
  }
  *&b = 20.0;
  *&k_p = b*g_r;
}

void param_K_r(gsl_rng *R){
  
  char filename[10] = "kr.dat";
  for(i=0;i<n_elements;i++){
    *&k_r = 0.0025*(i+1);
    cells(n_cells, t_total, R, filename);
  }
  *&k_r = 0.01;
}

void param_Tau_p(gsl_rng *R){
  
  char filename[10] = "taup.dat";
  for(i=0;i<n_elements;i++){
    *&tau_p = 900.0*(i+1);
    *&g_p = log(2)/tau_p;
    cells(n_cells, t_total, R, filename);
  }
  *&tau_p = 3600.0;
  *&g_p = log(2)/tau_p;
}
