#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>

const gsl_rng_type *T;
gsl_rng *R;

int n1_max = 4;
double k1p = 0.15;
double k1m = 0.05;
double tau1;
double k2 = 1.0;
double tau2 = 5.0;
double k3 = 60.0;
double tau3 = 30.0;

void gill_step(double *t, int *n1, int *n2, int *n3);
void cell(double t_total);
int main(){

  T = gsl_rng_mt19937;
  R = gsl_rng_alloc(T);
  gsl_rng_set(R,0);

  *&tau1 = 1.0/(k1p+k1m);
  
  cell(4.0);

  return 0;
}
void gill_step(double *t, int *n1, int *n2, int *n3){

  double K = k1p*(n1_max-*n1) + k1m*(*n1) + k2**n1 + *n2/tau2 + k3**n2 + *n3/tau3;

  double which = gsl_rng_uniform(R);
  double Dt = -log(gsl_rng_uniform(R))/K;
  
  *t+=Dt;

  if(which < k1p*(n1_max-*n1)){
    *n1+=1;
  } else if(which >= k1p*(n1_max-*n1) && which < k1p*(n1_max-*n1) + k1m*(*n1)){
    *n1-=1;
  } else if(which >= k1p*(n1_max-*n1) + k1m*(*n1) && which < k1p*(n1_max-*n1) + k1m*(*n1) + k2**n1){
    *n2+=1;
  } else if(which >= k1p*(n1_max-*n1) + k1m*(*n1) + k2**n1 && which < k1p*(n1_max-*n1) + k1m*(*n1) + k2**n1 + *n2/tau2){
    *n2-=1;
  } else if(which >= k1p*(n1_max-*n1) + k1m*(*n1) + k2**n1 + *n2/tau2 && which < k1p*(n1_max-*n1) + k1m*(*n1) + k2**n1 + *n2/tau2 + k3**n2){
    *n3+=1;
  } else {
    *n3-=1;
  }
}

void cell(double t_total){

  FILE *out = fopen("tave.dat", "w");
  
  double t = 0.0;
  int n1 = 0;
  int n2 = 0;
  int n3 = 0;


  fprintf(out, "%f %d %d %d\n", t, n1, n2, n3);
  while(t<t_total){
    gill_step(&t, &n1, &n2, &n3);
    printf("%f\n",t);
    fprintf(out, "%f %d %d %d\n", t, n1, n2, n3);
  }  
  fclose(out);  
}
