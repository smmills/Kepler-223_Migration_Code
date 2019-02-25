//  Written by Cezary Migazewski - Last Edit 1/1/16

void equations (my_float X[], my_float t, struct params_str params, my_float Xdot[])
{
 int i, j, ii;
 int nn = N+1;
 my_float x[nn], y[nn], z[nn], u[nn], v[nn], w[nn];
 my_float Omx[nn], Omy[nn], Omz[nn];
 my_float r2[nn], r[nn], r3[nn], v2[nn], h[nn];
 my_float rv[nn], hx[nn], hy[nn], hz[nn];
 my_float Om2[nn], Omr[nn], Omrd[nn], Omr0[nn], Omrd0[nn];
 my_float Omrx[nn], Omry[nn], Omrz[nn];
 my_float Omrx0[nn], Omry0[nn], Omrz0[nn];
 my_float xx, yy, zz, rr2, rr, rr3;
 my_float fx[nn][nn], fy[nn][nn], fz[nn][nn];
 my_float pom, pom1, pom2, pom3, pom4;
 my_float Fx[nn], Fy[nn], Fz[nn];
 my_float FOmx[nn], FOmy[nn], FOmz[nn];
 my_float vK, l, tau_a, tau_e, V, kappa;
 my_float pom1_in, pom2_in, tau_a_in, tau_e_in;


   for (i=1; i<=N; i++) {
       ii    = 4*(i-1);
       x[i]  = X[ii];
       y[i]  = X[ii+1];
       u[i]  = X[ii+2];
       v[i]  = X[ii+3];
       r2[i] = x[i]*x[i] + y[i]*y[i];
       r[i]  = sqrt(r2[i]);
       r3[i] = r[i]*r2[i];
       v2[i] = u[i]*u[i] + v[i]*v[i];
       rv[i] = x[i]*u[i] + y[i]*v[i];
       pom3  = x[i]*v[i] - y[i]*u[i];
       h[i]  = fabs(pom3);
   }
   
   for (i=1; i<=N-1; i++)
       for (j=i+1; j<=N; j++) {
           xx = x[i] - x[j];
	   yy = y[i] - y[j];
	   
	   rr2 = xx*xx + yy*yy;
	   rr  = sqrt(rr2);
	   rr3 = rr*rr2;
	   
	   // force acting on body i from body j
	   pom = -k2*m[j]/rr3;
	   fx[i][j]  = pom*xx;
	   fy[i][j]  = pom*yy;
       }

   for (i=2; i<=N; i++)
       for (j=1; j<=i-1; j++) {
           pom = m[j]/m[i];
           fx[i][j] = -fx[j][i]*pom;
	   fy[i][j] = -fy[j][i]*pom;
       }
       
   for (i=1; i<=N; i++) {
   
       Fx[i] = 0.0;
       Fy[i] = 0.0;
   
       for (j=1; j<=N; j++)
           if (j != i) {
	       Fx[i] += fx[i][j];
	       Fy[i] += fy[i][j];
	       pom = -k2*m[j]/r3[j];
	       Fx[i] += pom*x[j];
	       Fy[i] += pom*y[j];
	   }

       tau_a    = tau0*exp(t/T)*pow(r[i], ll)*pow(m[i]/mass_ref, ll2);
       tau_e    = tau_a/Kappa;
       vK       = sqrt(mu[i]/r[i]);
       V        = vK/(h[i]*r[i]);
       pom1     = 1.0/(2.0*tau_a) + 1.0/tau_e*(1.0 - V*r2[i]);
       pom2     = V*rv[i]/tau_e;
       Fx[i]   += -(pom1*u[i] + pom2*x[i]);
       Fy[i]   += -(pom1*v[i] + pom2*y[i]);

       pom = -mu[i]/r3[i];
       Fx[i] += pom*x[i];
       Fy[i] += pom*y[i];
   }


   for (i=1; i<=N; i++) {
       ii = 4*(i-1);
       Xdot[ii]   = u[i];
       Xdot[ii+1] = v[i];
       Xdot[ii+2] = Fx[i];
       Xdot[ii+3] = Fy[i];
   }


}
