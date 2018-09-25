void xxx(x,crea_par)
double x[1];
double crea_par[2];
{
  double t1;
  double t2;
  double t4;
  double t5;
  double t6;
  double t7;
  double t9;
  {
    t1 = a[0];
    t2 = a[1];
    t6 = x[0];
    t4 = a[3]*t6;
    t5 = a[2];
    t7 = (t4+t5)*t6;
    t9 = (t2+t7)*t6;
    crea_par[0] = (t1+t9)*t6;
    crea_par[1] = ((2.0*t4+t5)*t6+t2+t7)*t6+t1+t9;
    return;
  }
}

