//  Written by Cezary Migazewski - Last Edit 1/1/16

my_float kepler_equation (my_float e, my_float M)
{
 my_float E, E0, f0, f1, f2, f3, d1, d2, d3, delta, nu;
 int i;

  E = M;
  
  delta=1;
 
  i=0;
  while ((delta>1e-8) && (i<100)) {
     ++i;
     E0=E;
     f0=E-e*sin(E)-M;
     f1=1-e*cos(E);
     f2=e*sin(E);
     f3=e*cos(E);
     d1=-f0/f1;
     d2=-f0/(f1+d1*f2/2);
     d3=-f0/(f1+d2*f2/2+d2*d2*f3/6);
     E=E+d3;
     delta=fabs(E0-E);
  }
 
  nu = 2*atan2(sqrt(1 + e)*sin(E/2.), sqrt(1 - e)*cos(E/2.));

 return nu;

}
