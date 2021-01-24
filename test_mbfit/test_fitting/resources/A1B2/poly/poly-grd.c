void xxx(x,crea_par)
double x[3];
double crea_par[4];
{
  double t10;
  double t11;
  double t13;
  double t14;
  double t15;
  double t2;
  double t3;
  double t4;
  double t5;
  double t6;
  double t7;
  double t8;
  double t9;
  {
    t4 = x[2];
    t2 = a[5]*t4;
    t3 = a[1];
    t6 = a[3];
    t5 = x[1];
    t7 = t6*t5;
    t8 = a[2];
    t9 = t8*t4;
    t10 = a[0];
    t11 = x[0];
    t13 = t6*t11;
    t14 = a[4];
    t15 = t14*t5;
    crea_par[0] = (t2+t3)*t4+(t7+t9+t10)*t5+(t13+t15+t9+t10)*t11;
    crea_par[1] = 2.0*t13+t15+t9+t10;
    crea_par[2] = t11*t14+t10+2.0*t7+t9;
    crea_par[3] = t11*t8+t5*t8+2.0*t2+t3;
    return;
  }
}

