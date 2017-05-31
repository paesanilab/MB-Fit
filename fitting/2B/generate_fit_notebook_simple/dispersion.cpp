
#include "dispersion.h"

x2b_disp::x2b_disp() {}
x2b_disp::~x2b_disp() {
  delete[] xyz1;
  delete[] xyz2;
}

x2b_disp::x2b_disp(double * c1, double * c2, size_t n1, size_t n2) {
  xyz1 = new double[3*n1];
  xyz2 = new double[3*n2];
  std::copy(c1, c1 + 3*n1, xyz1);
  std::copy(c2, c2 + 3*n2, xyz2);
}

double x2b_disp::get_dispersion() {

  double disp = 0.0;
  const double* A_1_a= xyz1 + 0;
  const double* B_1_a= xyz1 + 3;
  const double* B_2_a= xyz1 + 6;

  const double* A_1_b= xyz2 + 0;
  const double* B_1_b= xyz2 + 3;
  const double* B_2_b= xyz2 + 6;

  disp += x6(m_C6_AA, m_d6_AA, m_C8, m_d8, A_1_a, A_1_b);
  disp += x6(m_C6_AB, m_d6_AB, m_C8, m_d8, A_1_a, B_1_b);
  disp += x6(m_C6_AB, m_d6_AB, m_C8, m_d8, A_1_a, B_2_b);

  disp += x6(m_C6_AB, m_d6_AB, m_C8, m_d8, B_1_a, A_1_b);
  disp += x6(m_C6_BB, m_d6_BB, m_C8, m_d8, B_1_a, B_1_b);
  disp += x6(m_C6_BB, m_d6_BB, m_C8, m_d8, B_1_a, B_2_b);

  disp += x6(m_C6_AB, m_d6_AB, m_C8, m_d8, B_2_a, A_1_b);
  disp += x6(m_C6_BB, m_d6_BB, m_C8, m_d8, B_2_a, B_1_b);
  disp += x6(m_C6_BB, m_d6_BB, m_C8, m_d8, B_2_a, B_2_b);


  return disp;
}
