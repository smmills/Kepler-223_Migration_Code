//  Written by Cezary Migazewski - Last Edit 1/1/16

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef double my_float;

#define k         0.01720209895
#define MAXN      10
#define year      365.242190419

struct params_str
{
 int dim;
};

my_float period_part;
int N, lp;
my_float m[MAXN+1], mu[MAXN+1], beta[MAXN+1];
my_float k2 = k*k;
my_float tau0, T0, ll, T, Kappa;
my_float mass_ref = 3.00406467116955e-06, ll2;
int Nfixedpoint=20;
my_float aa[5][5], bb[5], sqrt21;

#include "equations2d_migration.h"
#include "lobatto_params.h"
#include "lobattoIIIAs5.h"
#include "kepler_equation.h"
#include "integrate.h"
#include "kepl2xy_M_astro.h"
#include "xy2kepl_2d_astro_M.h"

main()
{

 char name1[100];
 FILE *file1;
 int i, bla, Nstep, j, ii;
 my_float kepl_P[6*MAXN], x0[6*MAXN]; 
 my_float X[4*MAXN], X1[4*MAXN];
 my_float h, t, tend, Amin;
 struct params_str params;

 sprintf(name1, "start_migration.in");
 file1 = fopen(name1, "r");
 bla = fscanf(file1, "%d %d %le %d %le", &lp, &N, &tend, &Nstep, &period_part);
 tend = tend*year;
 bla = fscanf(file1, "%le", &m[0]);
 for (i=0; i<N; i++) bla = fscanf(file1, "%le %le %le %le %le", &m[i+1], &kepl_P[6*i], &kepl_P[6*i+1], &kepl_P[6*i+4], &kepl_P[6*i+5]);
 bla = fscanf(file1, "%le %le %le %le %le", &T, &tau0, &ll, &Kappa, &ll2);
 fclose(file1);
 
 T *= year;
 tau0 *= year;
 
 lobatto_params();
 
 for (i=0; i<N; i++) {
     kepl_P[6*i+2] = 0.0;
     kepl_P[6*i+3] = 0.0;
     kepl_P[6*i+4] *= M_PI/180.0;
     kepl_P[6*i+5] *= M_PI/180.0;
 }
 
 params.dim = 4*N;
 kepl2xy(kepl_P, x0, N);
 
 for (i=0; i<=N-1; i++) {
     X[4*i]   = x0[6*i];
     X[4*i+1] = x0[6*i+1];
     X[4*i+2] = x0[6*i+3];
     X[4*i+3] = x0[6*i+4];
 }
 
 for (i=1; i<=N; i++) {
     mu[i]     = k2*(m[0] + m[i]);
     beta[i]   = m[0]*m[i]/(m[0] + m[i]);
 }

 t = 0.0;
 ii = 0;

 xy2kepl_2d_astro(x0, kepl_P, N);

 Amin = kepl_P[0];
 for (i=1; i<=N-1; i++) 
     if (kepl_P[6*i] < Amin)  Amin = kepl_P[6*i];
 h = 2*M_PI*sqrt(Amin*Amin*Amin/(k2*m[0]))/period_part;

 sprintf(name1, "evolution_%i.out", lp);
 file1 = fopen(name1, "w");

 fprintf(file1, "%.14le ", t/year);
 for (i=0; i<=N-1; i++) {
     for (j=0; j<=5; j++) fprintf(file1, "%.14le ", kepl_P[6*i+j]);
 }
 fprintf(file1, "\n");
 

 while (t < tend) {
     integrate(X, X1, t, t+h, params);
     for (i=0; i<=params.dim-1; i++) X[i] = X1[i];
     t += h;
     for (i=0; i<=N-1; i++) {
         x0[6*i]   = X[4*i];
	 x0[6*i+1] = X[4*i+1];
	 x0[6*i+2] = 0.0;
	 x0[6*i+3] = X[4*i+2];
	 x0[6*i+4] = X[4*i+3];
	 x0[6*i+5] = 0.0;
     }
     xy2kepl_2d_astro(x0, kepl_P, N);
     Amin = kepl_P[0];
     for (i=1; i<=N-1; i++)
         if (kepl_P[6*i] < Amin)  Amin = kepl_P[6*i];
     h = 2*M_PI*sqrt(Amin*Amin*Amin/(k2*m[0]))/period_part;
     ii++;
     if (ii == Nstep) {
         fprintf(file1, "%.14le ", t/year);
         for (i=0; i<=N-1; i++) {
             for (j=0; j<=5; j++) fprintf(file1, "%.14le ", kepl_P[6*i+j]);
         }
         fprintf(file1, "\n");
         ii = 0;
     }
 }
 
 fclose(file1);



 return 0;
}
