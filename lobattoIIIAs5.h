//  Written by Cezary Migazewski - Last Edit 1/1/16

void lobattoIIIAs5 ( void (*f)( my_float x[], my_float t1, struct params_str params1, my_float xdot[] ), my_float x0[], my_float x1[], my_float t, my_float h, struct params_str params )
{
 int i, nn;
 
 nn = params.dim;
 my_float k1[nn], k2[nn], k3[nn], k4[nn], k5[nn], xp[nn];
 my_float pomh1, pomh2, pomh3;
 my_float pom[10], tt[5];
 int ii=0;
 my_float diff, pomocne;
 
 tt[0] = t;
 tt[1] = t + (0.5 - sqrt21/14.0)*h;
 tt[2] = t + 0.5*h;
 tt[3] = t + (0.5 + sqrt21/14.0)*h;
 tt[4] = t + h;
  
 pomh1 = 0.5*h;
 
 (*f)(x0, t, params, k1);
 for (i=0; i<=nn-1; i++) xp[i] = x0[i] + pomh1*k1[i];
 t += pomh1;
 (*f)(xp, t, params, k2);
 for (i=0; i<=nn-1; i++) xp[i] = x0[i] + pomh1*k2[i];
 (*f)(xp, t, params, k3);
 for (i=0; i<=nn-1; i++) xp[i] = x0[i] + h*k3[i];
 t += pomh1;
 (*f)(xp, t, params, k4);
 
 for (i=0; i<=nn-1; i++) {
     k5[i] = k4[i];
     k3[i] = 0.5*(k2[i] + k3[i]);
     k2[i] = k1[i] + (tt[1] - tt[0])*(k3[i] - k1[i])/(tt[2] - tt[0]);
     k4[i] = k3[i] + (tt[3] - tt[2])*(k5[i] - k3[i])/(tt[4] - tt[2]);
 }
 
 
 pomocne = sqrt(k2[0]*k2[0] + k3[0]*k3[0] + k4[0]*k4[0] + k5[0]*k5[0] + k2[1]*k2[1] + k3[1]*k3[1] + k4[1]*k4[1] + k5[1]*k5[1]);
 diff = 1.0;
 
 while ((ii < Nfixedpoint) && (diff > 1.e-16)) {
     pom[0] = k2[0];
     pom[1] = k3[0];
     pom[2] = k4[0];
     pom[3] = k5[0];
     pom[4] = k2[1];
     pom[5] = k3[1];
     pom[6] = k4[1];
     pom[7] = k5[1];
     for (i=0; i<nn; i++) xp[i] = x0[i] + h*(aa[1][0]*k1[i] + aa[1][1]*k2[i] + aa[1][2]*k3[i] + aa[1][3]*k4[i] + aa[1][4]*k5[i]);
     (*f)(xp, tt[1], params, k2);
     for (i=0; i<nn; i++) xp[i] = x0[i] + h*(aa[2][0]*k1[i] + aa[2][1]*k2[i] + aa[2][2]*k3[i] + aa[2][3]*k4[i] + aa[2][4]*k5[i]);
     (*f)(xp, tt[2], params, k3);
     for (i=0; i<nn; i++) xp[i] = x0[i] + h*(aa[3][0]*k1[i] + aa[3][1]*k2[i] + aa[3][2]*k3[i] + aa[3][3]*k4[i] + aa[3][4]*k5[i]);
     (*f)(xp, tt[3], params, k4);
     for (i=0; i<nn; i++) xp[i] = x0[i] + h*(aa[4][0]*k1[i] + aa[4][1]*k2[i] + aa[4][2]*k3[i] + aa[4][3]*k4[i] + aa[4][4]*k5[i]);
     (*f)(xp, tt[4], params, k5);
     diff = (fabs(pom[0] - k2[0]) + fabs(pom[1] - k3[0]) + fabs(pom[2] - k4[0]) + fabs(pom[3] - k5[0]) + fabs(pom[4] - k2[1]) + fabs(pom[5] - k3[1]) + fabs(pom[6] - k4[1]) + fabs(pom[7] - k5[1]))/pomocne;
     ii++;
 }

 for (i=0; i<=nn-1; i++)  x1[i] = x0[i] + h*(bb[0]*k1[i] + bb[1]*k2[i] + bb[2]*k3[i] + bb[3]*k4[i] + bb[4]*k5[i]);
 
}
