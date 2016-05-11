#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>

const gsl_rng_type *T;
gsl_rng *R;

int i;
int n_cells = 10000;
int t_total;
int t_div;

double k_r = 1.0;
double g_r = 1.0/5.0;
double g_p = 1.0/30.0;
double k_p = 50.0;

void step(double *t, double *t_gen, double *r, double *p);
void divide_ind(double *p, double r_daughter);
void cell(double t_total, double *r_f, double *p_f);
void cells(double n_cells, double t_total, char *filename);
int main(int argc, char **argv){
  
  *&t_total = 15.0/g_p;
  *&t_div = log(2)/g_p;

  T = gsl_rng_mt19937;
  R = gsl_rng_alloc(T);
  gsl_rng_set(R,atoi(argv[1]));
  
  double r_f, p_f;
  
  cell(t_total, &r_f, &p_f);
  
  return 0;
}

void step(double *t, double *t_gen, double *r, double *p){
  
  double K = k_r + g_r**r + k_p**r;
  
  double Dt = -log(gsl_rng_uniform(R))/K;
  double which = gsl_rng_uniform(R);
  
  *t+=Dt;
  *t_gen+=Dt;
  
  if(which < k_r/K){
    *r+=1;
  }else if(which >= k_r/K && which < (k_r+g_r**r)/K){
    *r-=1;
  }else {
    *p+=1;
  }
}

void divide_ind(double *p, double r_daughter){
 
  double count = 0.0;
  
  for(i=0;i<*p;i++){
    if(gsl_rng_uniform(R) < r_daughter){
      count++;
    }
  }
  *p = count;
}
    
void cell(double t_total, double *r_f, double *p_f){
  
  double t = 0;
  double t_gen = 0;
  double r = 0;
  double p = 0;
  
  while(t<t_total){
    while(t_gen < t_div){ 
      step(&t, &t_gen, &r, &p);
      printf("%f %f %f\n", t, r, p);
    }
    t_gen = 0;
    divide_ind(&p, 1.0/2.0);
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
