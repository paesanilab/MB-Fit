void xxx(x,crea_par)
double x[3];
double crea_par[4];
{
  double t1;
  double t10;
  double t12;
  double t13;
  double t14;
  double t15;
  double t16;
  double t2;
  double t22;
  double t3;
  double t4;
  double t5;
  double t6;
  double t7;
  double t8;
  {
    t1 = a[2];
    t4 = x[2];
    t2 = t1*t4;
    t3 = a[0];
    t5 = x[1];
    t6 = t1*t5;
    t7 = a[3];
    t8 = t7*t4;
    t10 = x[0];
    t12 = a[4]*t10;
    t13 = a[5];
    t14 = t13*t5;
    t15 = t13*t4;
    t16 = a[1];
    t22 = t13*t10;
    crea_par[0] = (t2+t3)*t4+(t6+t8+t3)*t5+(t12+t14+t15+t16)*t10;
    crea_par[1] = 2.0*t12+t14+t15+t16;
    crea_par[2] = t22+2.0*t6+t8+t3;
    crea_par[3] = t7*t5+2.0*t2+t22+t3;
    return;
  }
}

