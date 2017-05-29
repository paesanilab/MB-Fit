#include "poly_2b_A1B2_A1B2.h"

namespace mb_system_fit {

void poly_model::eval(const double x[15], double a[134])
{
    double p[134];
    p[0] = x[11] + x[14] + x[12] + x[13];
    p[1] = x[9] + x[7] + x[10] + x[8];
    p[2] = x[6];

    p[3] = x[12]*x[2] + x[11]*x[2] + x[14]*x[2] + x[12]*x[5] + x[13]*x[2] + x[13]*x[5] + x[14]*x[5] + x[11]*x[5];
    p[4] = x[6]*x[9] + x[6]*x[7] + x[10]*x[6] + x[6]*x[8];
    p[5] = x[4]*x[9] + x[1]*x[8] + x[0]*x[7] + x[3]*x[9] + x[10]*x[3] + x[0]*x[8] + x[10]*x[4] + x[1]*x[7];
    p[6] = x[14]*x[3] + x[11]*x[1] + x[0]*x[14] + x[12]*x[1] + x[11]*x[4] + x[0]*x[13] + x[13]*x[4] + x[12]*x[3];
    p[7] = x[10]*x[5] + x[2]*x[8] + x[2]*x[7] + x[5]*x[9];
    p[8] = x[5]*x[6] + x[2]*x[6];
    p[9] = x[10]*x[7] + x[7]*x[9] + x[10]*x[8] + x[8]*x[9];
    p[10] = x[14]*x[6] + x[11]*x[6] + x[12]*x[6] + x[13]*x[6];
    p[11] = x[7]*x[8] + x[10]*x[9];
    p[12] = x[10]*x[11] + x[14]*x[9] + x[10]*x[12] + x[11]*x[8] + x[13]*x[9] + x[13]*x[8] + x[14]*x[7] + x[12]*x[7];
    p[13] = x[13]*x[14] + x[12]*x[14] + x[11]*x[13] + x[11]*x[12];
    p[14] = x[4]*x[6] + x[0]*x[6] + x[3]*x[6] + x[1]*x[6];
    p[15] = x[7]*x[7] + x[8]*x[8] + x[9]*x[9] + x[10]*x[10];
    p[16] = x[12]*x[9] + x[11]*x[7] + x[12]*x[8] + x[13]*x[7] + x[11]*x[9] + x[14]*x[8] + x[10]*x[14] + x[10]*x[13];
    p[17] = x[10]*x[1] + x[4]*x[8] + x[0]*x[9] + x[3]*x[7];
    p[18] = x[5]*x[7] + x[2]*x[9] + x[5]*x[8] + x[10]*x[2];
    p[19] = x[11]*x[14] + x[12]*x[13];
    p[20] = x[11]*x[3] + x[0]*x[12] + x[0]*x[11] + x[13]*x[1] + x[12]*x[4] + x[13]*x[3] + x[14]*x[1] + x[14]*x[4];
    p[21] = x[1]*x[9] + x[4]*x[7] + x[3]*x[8] + x[0]*x[10];
    p[22] = x[13]*x[13] + x[11]*x[11] + x[14]*x[14] + x[12]*x[12];
    p[23] = x[6]*x[6];

    p[24] = x[6]*x[8]*x[9] + x[10]*x[6]*x[8] + x[6]*x[7]*x[9] + x[10]*x[6]*x[7];
    p[25] = x[12]*x[1]*x[1] + x[11]*x[1]*x[1] + x[14]*x[3]*x[3] + x[13]*x[4]*x[4] + x[0]*x[0]*x[13] + x[0]*x[0]*x[14] + x[11]*x[4]*x[4] + x[12]*x[3]*x[3];
    p[26] = x[14]*x[14]*x[8] + x[10]*x[13]*x[13] + x[12]*x[12]*x[8] + x[11]*x[11]*x[7] + x[11]*x[11]*x[9] + x[12]*x[12]*x[9] + x[13]*x[13]*x[7] + x[10]*x[14]*x[14];
    p[27] = x[12]*x[12]*x[1] + x[11]*x[11]*x[1] + x[0]*x[13]*x[13] + x[11]*x[11]*x[4] + x[0]*x[14]*x[14] + x[14]*x[14]*x[3] + x[13]*x[13]*x[4] + x[12]*x[12]*x[3];
    p[28] = x[14]*x[14]*x[4] + x[11]*x[11]*x[3] + x[0]*x[12]*x[12] + x[13]*x[13]*x[3] + x[13]*x[13]*x[1] + x[12]*x[12]*x[4] + x[14]*x[14]*x[1] + x[0]*x[11]*x[11];
    p[29] = x[14]*x[2]*x[8] + x[13]*x[2]*x[7] + x[11]*x[2]*x[7] + x[12]*x[2]*x[8] + x[11]*x[5]*x[9] + x[12]*x[5]*x[9] + x[10]*x[14]*x[5] + x[10]*x[13]*x[5];
    p[30] = x[10]*x[2]*x[9] + x[5]*x[7]*x[8];
    p[31] = x[4]*x[6]*x[7] + x[1]*x[6]*x[9] + x[3]*x[6]*x[8] + x[0]*x[10]*x[6];
    p[32] = x[11]*x[13]*x[3] + x[12]*x[14]*x[4] + x[13]*x[14]*x[1] + x[0]*x[11]*x[12];
    p[33] = x[1]*x[1]*x[6] + x[3]*x[3]*x[6] + x[0]*x[0]*x[6] + x[4]*x[4]*x[6];
    p[34] = x[12]*x[13]*x[2] + x[11]*x[14]*x[2] + x[11]*x[14]*x[5] + x[12]*x[13]*x[5];
    p[35] = x[10]*x[10]*x[6] + x[6]*x[8]*x[8] + x[6]*x[9]*x[9] + x[6]*x[7]*x[7];
    p[36] = x[10]*x[10]*x[3] + x[1]*x[8]*x[8] + x[0]*x[8]*x[8] + x[1]*x[7]*x[7] + x[10]*x[10]*x[4] + x[3]*x[9]*x[9] + x[0]*x[7]*x[7] + x[4]*x[9]*x[9];
    p[37] = x[10]*x[13]*x[7] + x[12]*x[8]*x[9] + x[11]*x[7]*x[9] + x[10]*x[14]*x[8];
    p[38] = x[12]*x[3]*x[8] + x[12]*x[1]*x[9] + x[0]*x[10]*x[13] + x[0]*x[10]*x[14] + x[13]*x[4]*x[7] + x[11]*x[4]*x[7] + x[14]*x[3]*x[8] + x[11]*x[1]*x[9];
    p[39] = x[10]*x[5]*x[9] + x[2]*x[7]*x[8];
    p[40] = x[14]*x[1]*x[5] + x[14]*x[2]*x[4] + x[12]*x[2]*x[4] + x[13]*x[1]*x[5] + x[13]*x[2]*x[3] + x[0]*x[12]*x[5] + x[0]*x[11]*x[5] + x[11]*x[2]*x[3];
    p[41] = x[5]*x[5]*x[9] + x[10]*x[5]*x[5] + x[2]*x[2]*x[8] + x[2]*x[2]*x[7];
    p[42] = x[10]*x[3]*x[8] + x[1]*x[7]*x[9] + x[0]*x[10]*x[8] + x[10]*x[4]*x[7] + x[1]*x[8]*x[9] + x[3]*x[8]*x[9] + x[4]*x[7]*x[9] + x[0]*x[10]*x[7];
    p[43] = x[12]*x[13]*x[6] + x[11]*x[14]*x[6];
    p[44] = x[13]*x[5]*x[8] + x[14]*x[5]*x[7] + x[14]*x[2]*x[9] + x[10]*x[12]*x[2] + x[10]*x[11]*x[2] + x[11]*x[5]*x[8] + x[13]*x[2]*x[9] + x[12]*x[5]*x[7];
    p[45] = x[5]*x[6]*x[8] + x[10]*x[2]*x[6] + x[5]*x[6]*x[7] + x[2]*x[6]*x[9];
    p[46] = x[14]*x[5]*x[8] + x[11]*x[5]*x[7] + x[10]*x[14]*x[2] + x[12]*x[5]*x[8] + x[12]*x[2]*x[9] + x[11]*x[2]*x[9] + x[13]*x[5]*x[7] + x[10]*x[13]*x[2];
    p[47] = x[0]*x[8]*x[9] + x[10]*x[1]*x[7] + x[10]*x[4]*x[8] + x[10]*x[1]*x[8] + x[10]*x[3]*x[7] + x[3]*x[7]*x[9] + x[4]*x[8]*x[9] + x[0]*x[7]*x[9];
    p[48] = x[10]*x[10]*x[2] + x[2]*x[9]*x[9] + x[5]*x[8]*x[8] + x[5]*x[7]*x[7];
    p[49] = x[14]*x[2]*x[7] + x[13]*x[5]*x[9] + x[13]*x[2]*x[8] + x[10]*x[12]*x[5] + x[10]*x[11]*x[5] + x[12]*x[2]*x[7] + x[11]*x[2]*x[8] + x[14]*x[5]*x[9];
    p[50] = x[11]*x[11]*x[12] + x[13]*x[13]*x[14] + x[11]*x[12]*x[12] + x[12]*x[14]*x[14] + x[11]*x[11]*x[13] + x[11]*x[13]*x[13] + x[12]*x[12]*x[14] + x[13]*x[14]*x[14];
    p[51] = x[0]*x[3]*x[7] + x[1]*x[3]*x[7] + x[0]*x[4]*x[8] + x[10]*x[1]*x[4] + x[10]*x[1]*x[3] + x[0]*x[4]*x[9] + x[1]*x[4]*x[8] + x[0]*x[3]*x[9];
    p[52] = x[10]*x[13]*x[9] + x[14]*x[7]*x[8] + x[13]*x[7]*x[8] + x[10]*x[14]*x[9] + x[10]*x[11]*x[9] + x[12]*x[7]*x[8] + x[11]*x[7]*x[8] + x[10]*x[12]*x[9];
    p[53] = x[14]*x[3]*x[7] + x[0]*x[14]*x[9] + x[0]*x[13]*x[9] + x[10]*x[12]*x[1] + x[13]*x[4]*x[8] + x[11]*x[4]*x[8] + x[12]*x[3]*x[7] + x[10]*x[11]*x[1];
    p[54] = x[14]*x[7]*x[9] + x[13]*x[8]*x[9] + x[10]*x[11]*x[8] + x[10]*x[12]*x[7];
    p[55] = x[2]*x[4]*x[6] + x[0]*x[5]*x[6] + x[1]*x[5]*x[6] + x[2]*x[3]*x[6];
    p[56] = x[0]*x[10]*x[9] + x[3]*x[7]*x[8] + x[10]*x[1]*x[9] + x[4]*x[7]*x[8];
    p[57] = x[2]*x[7]*x[7] + x[2]*x[8]*x[8] + x[10]*x[10]*x[5] + x[5]*x[9]*x[9];
    p[58] = x[12]*x[14]*x[3] + x[11]*x[13]*x[4] + x[0]*x[13]*x[14] + x[11]*x[12]*x[1];
    p[59] = x[14]*x[14]*x[6] + x[11]*x[11]*x[6] + x[12]*x[12]*x[6] + x[13]*x[13]*x[6];
    p[60] = x[2]*x[5]*x[6];
    p[61] = x[11]*x[1]*x[8] + x[14]*x[3]*x[9] + x[0]*x[13]*x[8] + x[13]*x[4]*x[9] + x[0]*x[14]*x[7] + x[12]*x[1]*x[7] + x[10]*x[12]*x[3] + x[10]*x[11]*x[4];
    p[62] = x[11]*x[1]*x[2] + x[11]*x[4]*x[5] + x[0]*x[14]*x[2] + x[13]*x[4]*x[5] + x[14]*x[3]*x[5] + x[12]*x[3]*x[5] + x[12]*x[1]*x[2] + x[0]*x[13]*x[2];
    p[63] = x[0]*x[11]*x[9] + x[10]*x[14]*x[1] + x[14]*x[4]*x[8] + x[12]*x[4]*x[8] + x[10]*x[13]*x[1] + x[0]*x[12]*x[9] + x[13]*x[3]*x[7] + x[11]*x[3]*x[7];
    p[64] = x[2]*x[6]*x[6] + x[5]*x[6]*x[6];
    p[65] = x[12]*x[14]*x[1] + x[0]*x[11]*x[13] + x[11]*x[12]*x[4] + x[11]*x[13]*x[1] + x[11]*x[12]*x[3] + x[13]*x[14]*x[3] + x[0]*x[12]*x[14] + x[13]*x[14]*x[4];
    p[66] = x[13]*x[7]*x[9] + x[12]*x[7]*x[9] + x[10]*x[11]*x[7] + x[10]*x[14]*x[7] + x[10]*x[13]*x[8] + x[11]*x[8]*x[9] + x[14]*x[8]*x[9] + x[10]*x[12]*x[8];
    p[67] = x[13]*x[6]*x[8] + x[11]*x[6]*x[8] + x[13]*x[6]*x[9] + x[10]*x[11]*x[6] + x[14]*x[6]*x[9] + x[10]*x[12]*x[6] + x[12]*x[6]*x[7] + x[14]*x[6]*x[7];
    p[68] = x[3]*x[3]*x[8] + x[4]*x[4]*x[7] + x[1]*x[1]*x[9] + x[0]*x[0]*x[10];
    p[69] = x[3]*x[4]*x[9] + x[0]*x[1]*x[7] + x[0]*x[1]*x[8] + x[10]*x[3]*x[4];
    p[70] = x[12]*x[4]*x[7] + x[0]*x[10]*x[12] + x[13]*x[3]*x[8] + x[13]*x[1]*x[9] + x[14]*x[4]*x[7] + x[0]*x[10]*x[11] + x[14]*x[1]*x[9] + x[11]*x[3]*x[8];
    p[71] = x[13]*x[1]*x[2] + x[12]*x[4]*x[5] + x[14]*x[1]*x[2] + x[0]*x[12]*x[2] + x[0]*x[11]*x[2] + x[13]*x[3]*x[5] + x[14]*x[4]*x[5] + x[11]*x[3]*x[5];
    p[72] = x[4]*x[6]*x[8] + x[0]*x[6]*x[9] + x[3]*x[6]*x[7] + x[10]*x[1]*x[6];
    p[73] = x[12]*x[13]*x[7] + x[10]*x[11]*x[14] + x[12]*x[13]*x[9] + x[11]*x[14]*x[9] + x[12]*x[13]*x[8] + x[11]*x[14]*x[8] + x[10]*x[12]*x[13] + x[11]*x[14]*x[7];
    p[74] = x[1]*x[3]*x[8] + x[1]*x[4]*x[7] + x[0]*x[10]*x[3] + x[0]*x[4]*x[7] + x[1]*x[4]*x[9] + x[1]*x[3]*x[9] + x[0]*x[10]*x[4] + x[0]*x[3]*x[8];
    p[75] = x[10]*x[5]*x[8] + x[10]*x[5]*x[7] + x[5]*x[7]*x[9] + x[2]*x[8]*x[9] + x[5]*x[8]*x[9] + x[10]*x[2]*x[8] + x[2]*x[7]*x[9] + x[10]*x[2]*x[7];
    p[76] = x[14]*x[1]*x[4] + x[0]*x[12]*x[4] + x[13]*x[1]*x[3] + x[0]*x[11]*x[3];
    p[77] = x[10]*x[14]*x[3] + x[0]*x[14]*x[8] + x[11]*x[4]*x[9] + x[0]*x[13]*x[7] + x[11]*x[1]*x[7] + x[10]*x[13]*x[4] + x[12]*x[1]*x[8] + x[12]*x[3]*x[9];
    p[78] = x[11]*x[13]*x[8] + x[10]*x[11]*x[12] + x[13]*x[14]*x[9] + x[12]*x[14]*x[7];
    p[79] = x[10]*x[10]*x[12] + x[14]*x[7]*x[7] + x[14]*x[9]*x[9] + x[13]*x[9]*x[9] + x[10]*x[10]*x[11] + x[12]*x[7]*x[7] + x[11]*x[8]*x[8] + x[13]*x[8]*x[8];
    p[80] = x[14]*x[1]*x[8] + x[0]*x[11]*x[7] + x[11]*x[3]*x[9] + x[0]*x[12]*x[8] + x[13]*x[1]*x[7] + x[10]*x[14]*x[4] + x[10]*x[13]*x[3] + x[12]*x[4]*x[9];
    p[81] = x[0]*x[13]*x[4] + x[11]*x[1]*x[4] + x[12]*x[1]*x[3] + x[0]*x[14]*x[3];
    p[82] = x[10]*x[10]*x[7] + x[8]*x[8]*x[9] + x[10]*x[10]*x[8] + x[10]*x[7]*x[7] + x[8]*x[9]*x[9] + x[7]*x[7]*x[9] + x[7]*x[9]*x[9] + x[10]*x[8]*x[8];
    p[83] = x[12]*x[14]*x[5] + x[11]*x[12]*x[2] + x[11]*x[13]*x[5] + x[13]*x[14]*x[2];
    p[84] = x[10]*x[10]*x[13] + x[10]*x[10]*x[14] + x[13]*x[7]*x[7] + x[11]*x[7]*x[7] + x[11]*x[9]*x[9] + x[12]*x[9]*x[9] + x[14]*x[8]*x[8] + x[12]*x[8]*x[8];
    p[85] = x[0]*x[2]*x[9] + x[10]*x[1]*x[2] + x[4]*x[5]*x[8] + x[3]*x[5]*x[7];
    p[86] = x[0]*x[13]*x[5] + x[14]*x[2]*x[3] + x[11]*x[1]*x[5] + x[11]*x[2]*x[4] + x[0]*x[14]*x[5] + x[12]*x[2]*x[3] + x[12]*x[1]*x[5] + x[13]*x[2]*x[4];
    p[87] = x[6]*x[6]*x[9] + x[6]*x[6]*x[8] + x[10]*x[6]*x[6] + x[6]*x[6]*x[7];
    p[88] = x[10]*x[4]*x[6] + x[10]*x[3]*x[6] + x[1]*x[6]*x[7] + x[0]*x[6]*x[8] + x[3]*x[6]*x[9] + x[4]*x[6]*x[9] + x[0]*x[6]*x[7] + x[1]*x[6]*x[8];
    p[89] = x[14]*x[6]*x[6] + x[13]*x[6]*x[6] + x[11]*x[6]*x[6] + x[12]*x[6]*x[6];
    p[90] = x[14]*x[3]*x[4] + x[13]*x[3]*x[4] + x[0]*x[14]*x[1] + x[0]*x[12]*x[1] + x[0]*x[11]*x[1] + x[11]*x[3]*x[4] + x[12]*x[3]*x[4] + x[0]*x[13]*x[1];
    p[91] = x[12]*x[13]*x[14] + x[11]*x[12]*x[14] + x[11]*x[13]*x[14] + x[11]*x[12]*x[13];
    p[92] = x[12]*x[2]*x[6] + x[14]*x[2]*x[6] + x[12]*x[5]*x[6] + x[11]*x[5]*x[6] + x[13]*x[2]*x[6] + x[11]*x[2]*x[6] + x[13]*x[5]*x[6] + x[14]*x[5]*x[6];
    p[93] = x[0]*x[5]*x[7] + x[2]*x[3]*x[9] + x[1]*x[5]*x[8] + x[0]*x[5]*x[8] + x[1]*x[5]*x[7] + x[10]*x[2]*x[4] + x[10]*x[2]*x[3] + x[2]*x[4]*x[9];
    p[94] = x[12]*x[1]*x[6] + x[11]*x[1]*x[6] + x[13]*x[4]*x[6] + x[11]*x[4]*x[6] + x[14]*x[3]*x[6] + x[12]*x[3]*x[6] + x[0]*x[13]*x[6] + x[0]*x[14]*x[6];
    p[95] = x[0]*x[11]*x[14] + x[12]*x[13]*x[3] + x[11]*x[14]*x[3] + x[12]*x[13]*x[4] + x[11]*x[14]*x[1] + x[0]*x[12]*x[13] + x[11]*x[14]*x[4] + x[12]*x[13]*x[1];
    p[96] = x[0]*x[14]*x[4] + x[13]*x[1]*x[4] + x[11]*x[1]*x[3] + x[12]*x[1]*x[4] + x[0]*x[12]*x[3] + x[0]*x[11]*x[4] + x[0]*x[13]*x[3] + x[14]*x[1]*x[3];
    p[97] = x[12]*x[14]*x[2] + x[11]*x[13]*x[2] + x[11]*x[12]*x[5] + x[13]*x[14]*x[5];
    p[98] = x[5]*x[6]*x[9] + x[2]*x[6]*x[8] + x[2]*x[6]*x[7] + x[10]*x[5]*x[6];
    p[99] = x[10]*x[8]*x[9] + x[7]*x[8]*x[9] + x[10]*x[7]*x[9] + x[10]*x[7]*x[8];
    p[100] = x[13]*x[5]*x[5] + x[11]*x[2]*x[2] + x[14]*x[2]*x[2] + x[12]*x[5]*x[5] + x[12]*x[2]*x[2] + x[14]*x[5]*x[5] + x[11]*x[5]*x[5] + x[13]*x[2]*x[2];
    p[101] = x[10]*x[3]*x[5] + x[0]*x[2]*x[8] + x[1]*x[2]*x[8] + x[4]*x[5]*x[9] + x[0]*x[2]*x[7] + x[10]*x[4]*x[5] + x[1]*x[2]*x[7] + x[3]*x[5]*x[9];
    p[102] = x[4]*x[4]*x[9] + x[10]*x[4]*x[4] + x[3]*x[3]*x[9] + x[1]*x[1]*x[8] + x[0]*x[0]*x[8] + x[1]*x[1]*x[7] + x[10]*x[3]*x[3] + x[0]*x[0]*x[7];
    p[103] = x[6]*x[7]*x[8] + x[10]*x[6]*x[9];
    p[104] = x[8]*x[8]*x[8] + x[9]*x[9]*x[9] + x[7]*x[7]*x[7] + x[10]*x[10]*x[10];
    p[105] = x[12]*x[4]*x[4] + x[13]*x[3]*x[3] + x[13]*x[1]*x[1] + x[14]*x[1]*x[1] + x[0]*x[0]*x[11] + x[14]*x[4]*x[4] + x[11]*x[3]*x[3] + x[0]*x[0]*x[12];
    p[106] = x[0]*x[1]*x[6] + x[3]*x[4]*x[6];
    p[107] = x[2]*x[4]*x[8] + x[2]*x[3]*x[7] + x[0]*x[5]*x[9] + x[10]*x[1]*x[5];
    p[108] = x[10]*x[12]*x[12] + x[12]*x[12]*x[7] + x[14]*x[14]*x[7] + x[14]*x[14]*x[9] + x[13]*x[13]*x[9] + x[13]*x[13]*x[8] + x[11]*x[11]*x[8] + x[10]*x[11]*x[11];
    p[109] = x[12]*x[14]*x[6] + x[11]*x[12]*x[6] + x[13]*x[14]*x[6] + x[11]*x[13]*x[6];
    p[110] = x[1]*x[5]*x[9] + x[0]*x[10]*x[5] + x[2]*x[3]*x[8] + x[2]*x[4]*x[7];
    p[111] = x[1]*x[4]*x[6] + x[0]*x[4]*x[6] + x[0]*x[3]*x[6] + x[1]*x[3]*x[6];
    p[112] = x[2]*x[5]*x[9] + x[2]*x[5]*x[8] + x[2]*x[5]*x[7] + x[10]*x[2]*x[5];
    p[113] = x[13]*x[13]*x[2] + x[12]*x[12]*x[5] + x[14]*x[14]*x[2] + x[14]*x[14]*x[5] + x[11]*x[11]*x[5] + x[13]*x[13]*x[5] + x[12]*x[12]*x[2] + x[11]*x[11]*x[2];
    p[114] = x[4]*x[6]*x[6] + x[1]*x[6]*x[6] + x[3]*x[6]*x[6] + x[0]*x[6]*x[6];
    p[115] = x[0]*x[10]*x[10] + x[1]*x[9]*x[9] + x[4]*x[7]*x[7] + x[3]*x[8]*x[8];
    p[116] = x[10]*x[2]*x[2] + x[5]*x[5]*x[8] + x[5]*x[5]*x[7] + x[2]*x[2]*x[9];
    p[117] = x[10]*x[11]*x[3] + x[13]*x[3]*x[9] + x[0]*x[11]*x[8] + x[14]*x[4]*x[9] + x[14]*x[1]*x[7] + x[13]*x[1]*x[8] + x[0]*x[12]*x[7] + x[10]*x[12]*x[4];
    p[118] = x[12]*x[13]*x[13] + x[12]*x[12]*x[13] + x[11]*x[11]*x[14] + x[11]*x[14]*x[14];
    p[119] = x[10]*x[1]*x[1] + x[0]*x[0]*x[9] + x[3]*x[3]*x[7] + x[4]*x[4]*x[8];
    p[120] = x[11]*x[2]*x[5] + x[13]*x[2]*x[5] + x[14]*x[2]*x[5] + x[12]*x[2]*x[5];
    p[121] = x[2]*x[2]*x[6] + x[5]*x[5]*x[6];
    p[122] = x[10]*x[10]*x[9] + x[10]*x[9]*x[9] + x[7]*x[8]*x[8] + x[7]*x[7]*x[8];
    p[123] = x[12]*x[12]*x[12] + x[11]*x[11]*x[11] + x[14]*x[14]*x[14] + x[13]*x[13]*x[13];
    p[124] = x[4]*x[5]*x[6] + x[1]*x[2]*x[6] + x[3]*x[5]*x[6] + x[0]*x[2]*x[6];
    p[125] = x[13]*x[3]*x[6] + x[11]*x[3]*x[6] + x[14]*x[4]*x[6] + x[13]*x[1]*x[6] + x[0]*x[11]*x[6] + x[0]*x[12]*x[6] + x[12]*x[4]*x[6] + x[14]*x[1]*x[6];
    p[126] = x[10]*x[13]*x[14] + x[12]*x[14]*x[8] + x[11]*x[13]*x[7] + x[11]*x[12]*x[9];
    p[127] = x[3]*x[4]*x[8] + x[3]*x[4]*x[7] + x[0]*x[1]*x[9] + x[0]*x[10]*x[1];
    p[128] = x[10]*x[4]*x[9] + x[10]*x[3]*x[9] + x[1]*x[7]*x[8] + x[0]*x[7]*x[8];
    p[129] = x[10]*x[10]*x[1] + x[4]*x[8]*x[8] + x[3]*x[7]*x[7] + x[0]*x[9]*x[9];
    p[130] = x[13]*x[14]*x[8] + x[11]*x[13]*x[9] + x[11]*x[12]*x[7] + x[11]*x[12]*x[8] + x[12]*x[14]*x[9] + x[10]*x[12]*x[14] + x[13]*x[14]*x[7] + x[10]*x[11]*x[13];
    p[131] = x[1]*x[2]*x[9] + x[3]*x[5]*x[8] + x[4]*x[5]*x[7] + x[0]*x[10]*x[2];
    p[132] = x[11]*x[6]*x[7] + x[10]*x[14]*x[6] + x[12]*x[6]*x[8] + x[14]*x[6]*x[8] + x[10]*x[13]*x[6] + x[11]*x[6]*x[9] + x[13]*x[6]*x[7] + x[12]*x[6]*x[9];
    p[133] = x[6]*x[6]*x[6];

    for(int i = 0; i < 134; ++i)
        a[i] = p[i];

}
} // namespace mb_system
