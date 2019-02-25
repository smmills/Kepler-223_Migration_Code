//  Written by Cezary Migazewski - Last Edit 1/1/16

void xy2kepl_2d_astro ( my_float xy_astro[], my_float kepl_P[], int N )
{
 my_float x[N], y[N], z[N];
 my_float xd[N], yd[N], zd[N];
 my_float vx[N], vy[N], vz[N];
 my_float mu[N], beta[N], h[N], r[N], v2[N];
 my_float Vx, Vy, Vz, M;
 my_float a[N], e[N], I[N], Om[N], om[N], f[N];
 my_float sinOm, cosOm, sinomf, cosomf, sinf, cosf, rd, omf;
 my_float hx[N], hy[N], hz[N];
 int i, ii;
 my_float k2;
 my_float rv[N], E[N];
 my_float cosE, sinE, cosw, sinw, x0, y0, Mmean[N];
 
 k2 = k*k;

 for (i=0; i<=N-1; i++) {
    ii = 6*i;
    x[i]  = xy_astro[ii];
    y[i]  = xy_astro[ii+1];
    z[i]  = xy_astro[ii+2];
    xd[i] = xy_astro[ii+3];
    yd[i] = xy_astro[ii+4];
    zd[i] = xy_astro[ii+5];
    mu[i] = k2*(m[0]+m[i+1]);
    beta[i] = m[i+1]*m[0]/(m[0]+m[i+1]);
 }

 for (i=0; i<=N-1; i++) {
    vx[i] = xd[i];
    vy[i] = yd[i];
    vz[i] = zd[i];
 }
 
 for (i=0; i<=N-1; i++) {
    hx[i] = y[i]*vz[i] - z[i]*vy[i];
    hy[i] = z[i]*vx[i] - x[i]*vz[i];
    hz[i] = x[i]*vy[i] - y[i]*vx[i];
    h[i] = sqrt(hx[i]*hx[i] + hy[i]*hy[i] + hz[i]*hz[i]);
    r[i] = sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]);
    v2[i] = vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
    rv[i] = x[i]*vx[i] + y[i]*vy[i] + z[i]*vz[i];
 }
 
 for (i=0; i<=N-1; i++) {
    a[i] = 1/(2./r[i] - v2[i]/mu[i]);
    e[i] = sqrt(1. - h[i]*h[i]/(mu[i]*a[i]));
    cosE = (1.0 - r[i]/a[i])/e[i];
    sinE = rv[i]/(e[i]*sqrt(mu[i]*a[i]));
    E[i] = atan2(sinE, cosE);
    Mmean[i] = E[i] - e[i]*sinE;
    x0 = a[i]*(cosE - e[i]);
    y0 = a[i]*sqrt(1 - e[i]*e[i])*sinE;
    cosw = (x[i]*x0 + y[i]*y0)/(r[i]*r[i]);
    sinw = (x0*y[i] - y0*x[i])/(r[i]*r[i]);
    om[i] = atan2(sinw, cosw);
 }
 
 for (i=0; i<=N-1; i++) {
    ii =6*i;
    kepl_P[ii]   = a[i];
    kepl_P[ii+1] = e[i];
    kepl_P[ii+2] = 0.0;
    kepl_P[ii+3] = 0.0;
    kepl_P[ii+4] = om[i];
    kepl_P[ii+5] = Mmean[i];
 }

}
