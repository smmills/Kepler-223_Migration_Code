//  Written by Cezary Migazewski - Last Edit 1/1/16

void kepl2xy ( my_float kepl_P[], my_float xy_astro[], int N )
{
 my_float a[N], e[N], I[N], Om[N], om[N], f[N];
 my_float mu[N], beta[N];
 my_float P11, P12, P21, P22, P31, P32;
 my_float r, h, rd, X, Y, Xd, Yd;
 my_float x[N], y[N], z[N];
 my_float vx[N], vy[N], vz[N];
 my_float xd[N], yd[N], zd[N];
 my_float Vx, Vy, Vz, k2;
 int i, ii;
 

 k2 = k*k;

 for (i=0; i<=N-1; i++) {
    ii = 6*i;
    a[i]  = kepl_P[ii];
    e[i]  = kepl_P[ii+1];
    I[i]  = kepl_P[ii+2];
    Om[i] = kepl_P[ii+3];
    om[i] = kepl_P[ii+4];
    f[i]  = kepl_P[ii+5];  // Mmean
    f[i]  = kepler_equation(e[i], f[i]);
    mu[i] = k2*(m[0]+m[i+1]);
    beta[i] = m[i+1]*m[0]/(m[0]+m[i+1]);
 }
 
 
 for (i=0; i<=N-1; i++) {
 
    ii = 6*i;

    P11 = cos(Om[i])*cos(om[i]) - sin(Om[i])*cos(I[i])*sin(om[i]);
    P12 = -cos(Om[i])*sin(om[i]) - sin(Om[i])*cos(I[i])*cos(om[i]);
    P21 = sin(Om[i])*cos(om[i]) + cos(Om[i])*cos(I[i])*sin(om[i]);
    P22 = -sin(Om[i])*sin(om[i]) + cos(Om[i])*cos(I[i])*cos(om[i]);
    P31 = sin(I[i])*sin(om[i]);
    P32 = sin(I[i])*cos(om[i]);
    
    r = a[i]*(1.0 - e[i]*e[i]) / (1 + e[i]*cos(f[i]));
    h = sqrt(mu[i]*a[i]*(1.0 - e[i]*e[i]));
    rd = h*e[i]*sin(f[i])/(a[i]*(1.0 - e[i]*e[i]));
    
    X = r*cos(f[i]);
    Y = r*sin(f[i]);
    Xd = rd*cos(f[i]) - (h/r)*sin(f[i]);
    Yd = rd*sin(f[i]) + (h/r)*cos(f[i]);
    
    x[i] = P11*X + P12*Y;
    y[i] = P21*X + P22*Y;
    z[i] = P31*X + P32*Y;
    
    vx[i] = P11*Xd + P12*Yd;
    vy[i] = P21*Xd + P22*Yd;
    vz[i] = P31*Xd + P32*Yd;
    
    xy_astro[ii]   =  x[i];
    xy_astro[ii+1] =  y[i];
    xy_astro[ii+2] =  z[i];
    xy_astro[ii+3] = vx[i];
    xy_astro[ii+4] = vy[i];
    xy_astro[ii+5] = vz[i];
    
 }
 

}
