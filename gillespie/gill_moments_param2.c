#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>

int i;
int n_cells = 1000;
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

void step(double *t, double *r, double *p, gsl_rng *R);
void cell(double t_total, double *r_f, double *p_f, gsl_rng *R);
void cells(double n_cells, double t_total, gsl_rng *R, char *filename);
void param_n(gsl_rng *R);
void param_Kd(gsl_rng *R);
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
  
  param_n(R);
  param_Kd(R);

  return 0;
}

void step(double *t, double *r, double *p, gsl_rng *R){
  
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

void param_n(gsl_rng *R){
  
  char filename[10] = "n.dat";
  for(i=0;i<21;i++){
    *&n = i;
    cells(n_cells, t_total, R, filename);
  }
  *&n = 2.0;
}

void param_Kd(gsl_rng *R){
  
  char filename[10] = "Kd.dat";
  for(i=0;i<21;i++){
    *&Kd =100*i;
    cells(n_cells, t_total, R, filename);
  }
  *&Kd =3000;
  cells(n_cells, t_total, R, filename);
  *&Kd =4000;
  cells(n_cells, t_total, R, filename);
  *&Kd =5000;
  cells(n_cells, t_total, R, filename);
  *&Kd = 800.0;
}
