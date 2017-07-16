
#include "dispersion.h"

x2b_disp::x2b_disp() {
  xyz1 = new double[3];
  xyz2 = new double[3];
  grd1 = new double[3];
  grd2 = new double[3];
}
x2b_disp::~x2b_disp() {
  delete[] xyz1;
  delete[] xyz2;
  delete[] grd1;
  delete[] grd2;
}

x2b_disp::x2b_disp(double * c1, double * c2, size_t n1, size_t n2) {
  xyz1 = new double[3*n1];
  xyz2 = new double[3*n2];
  grd1 = new double[3*n1];
  grd2 = new double[3*n2];
  std::copy(c1, c1 + 3*n1, xyz1);
  std::copy(c2, c2 + 3*n2, xyz2);
}

x2b_disp::x2b_disp(double * c1, double * c2, 
                   double * g1, double * g2,
                   size_t n1, size_t n2) {
  xyz1 = new double[3*n1];
  xyz2 = new double[3*n2];
  grd1 = new double[3*n1];
  grd2 = new double[3*n2];
  std::copy(c1, c1 + 3*n1, xyz1);
  std::copy(c2, c2 + 3*n2, xyz2);
  std::copy(g1, g1 + 3*n1, grd1);
  std::copy(g2, g2 + 3*n2, grd2);
}

double x2b_disp::get_dispersion() {

  double disp = 0.0;
  const double* A_1_a= xyz1 + 0;
  const double* B_1_a= xyz1 + 3;
  const double* B_2_a= xyz1 + 6;

  const double* D_1_b= xyz2 + 0;
  const double* E_1_b= xyz2 + 3;
  const double* E_2_b= xyz2 + 6;

  disp += x6(m_C6_AD, m_d6_AD, m_C8, m_d8, A_1_a, D_1_b);
  disp += x6(m_C6_AE, m_d6_AE, m_C8, m_d8, A_1_a, E_1_b);
  disp += x6(m_C6_AE, m_d6_AE, m_C8, m_d8, A_1_a, E_2_b);

  disp += x6(m_C6_BD, m_d6_BD, m_C8, m_d8, B_1_a, D_1_b);
  disp += x6(m_C6_BE, m_d6_BE, m_C8, m_d8, B_1_a, E_1_b);
  disp += x6(m_C6_BE, m_d6_BE, m_C8, m_d8, B_1_a, E_2_b);

  disp += x6(m_C6_BD, m_d6_BD, m_C8, m_d8, B_2_a, D_1_b);
  disp += x6(m_C6_BE, m_d6_BE, m_C8, m_d8, B_2_a, E_1_b);
  disp += x6(m_C6_BE, m_d6_BE, m_C8, m_d8, B_2_a, E_2_b);


  return disp;
}

double x2b_disp::get_dispersion(double * grd) {

  double disp = 0.0;
  const double* A_1_a= xyz1 + 0;
  const double* B_1_a= xyz1 + 3;
  const double* B_2_a= xyz1 + 6;

  const double* D_1_b= xyz2 + 0;
  const double* E_1_b= xyz2 + 3;
  const double* E_2_b= xyz2 + 6;

  double* A_1_a_g= grd1 + 0;
  double* B_1_a_g= grd1 + 3;
  double* B_2_a_g= grd1 + 6;

  double* D_1_b_g= grd2 + 0;
  double* E_1_b_g= grd2 + 3;
  double* E_2_b_g= grd2 + 6;

  disp += x6(m_C6_AD, m_d6_AD, m_C8, m_d8, 
             A_1_a, D_1_b, A_1_a_g, D_1_b_g);
  disp += x6(m_C6_AE, m_d6_AE, m_C8, m_d8, 
             A_1_a, E_1_b, A_1_a_g, E_1_b_g);
  disp += x6(m_C6_AE, m_d6_AE, m_C8, m_d8, 
             A_1_a, E_2_b, A_1_a_g, E_2_b_g);

  disp += x6(m_C6_BD, m_d6_BD, m_C8, m_d8, 
             B_1_a, D_1_b, B_1_a_g, D_1_b_g);
  disp += x6(m_C6_BE, m_d6_BE, m_C8, m_d8, 
             B_1_a, E_1_b, B_1_a_g, E_1_b_g);
  disp += x6(m_C6_BE, m_d6_BE, m_C8, m_d8, 
             B_1_a, E_2_b, B_1_a_g, E_2_b_g);

  disp += x6(m_C6_BD, m_d6_BD, m_C8, m_d8, 
             B_2_a, D_1_b, B_2_a_g, D_1_b_g);
  disp += x6(m_C6_BE, m_d6_BE, m_C8, m_d8, 
             B_2_a, E_1_b, B_2_a_g, E_1_b_g);
  disp += x6(m_C6_BE, m_d6_BE, m_C8, m_d8, 
             B_2_a, E_2_b, B_2_a_g, E_2_b_g);


  return disp;
}
