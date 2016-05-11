#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>

const gsl_rng_type *T;
gsl_rng *R;

int i;
int n_cells = 10000;
int t_total;

double k_r = 1.0;
double g_r = 1.0/5.0;
double g_p = 1.0/30.0;
double k_p = 50.0;

void step(double *t, double *r, double *p);
void cell(double t_total, double *r_f, double *p_f);
void cells(double n_cells, double t_total, char *filename);
int main(int argc, char **argv){
  
  *&t_total = 10.0/g_p;
  
  T = gsl_rng_mt19937;
  R = gsl_rng_alloc(T);
  gsl_rng_set(R,atoi(argv[1]));
  
  double r_f, p_f;
  
  cell(t_total, &r_f, &p_f);
  
  return 0;
}

void step(double *t, double *r, double *p){
  
  double K = k_r + g_r**r + k_p**r + g_p**p;
  
  double Dt = -log(gsl_rng_uniform(R))/K;
  double which = gsl_rng_uniform(R);
  
  *t+=Dt;
  
  if(which < k_r/K){
    *r+=1;
  }else if(which >= k_r/K && which < (k_r+g_r**r)/K){
    *r-=1;
  }else if(which >= (k_r+g_r**r)/K && which < (k_r+g_r**r+k_p**r)/K) {
    *p+=1;
  } else {
    *p-=1;
  }
}

void cell(double t_total, double *r_f, double *p_f){
  
  double t = 0;
  double r = 0;
  double p = 0;
  
  while(t<t_total){
    step(&t, &r, &p);
    printf("%f %f %f\n", t, r, p);
  }
  *r_f = r;
  *p_f = p;
}

void cells(double n_cells, double t_total, char *filename){
  int i;
  double ave_p = 0;
  double ave_p2 = 0;

  double r_f;
  double p_f;

  /*FILE *out = fopen(filename, "a");*/
  
  for(i=0;i<n_cells;i++){
    cell(t_total, &r_f, &p_f);
    ave_p += p_f;
    ave_p2 += p_f*p_f;
  }
  
  ave_p/=n_cells;
  ave_p2/=n_cells;
  
  /*fprintf(out, "%f %f\n", ave_p, (ave_p2 - ave_p*ave_p)/ave_p);
    fclose(out);*/
}
