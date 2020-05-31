#include "poly-model.h"

namespace mb_system {

double poly_model::eval(const double a[24], const double x[15])
{
    const double t1 = a[13];
    const double t3 = a[6];
    const double t2 = x[14];
    const double t4 = t3*t2;
    const double t5 = a[0];
    const double t9 = a[10];
    const double t6 = x[13];
    const double t10 = t9*t6;
    const double t14 = a[17];
    const double t27 = a[3];
    const double t7 = x[10];
    const double t28 = t27*t7;
    const double t29 = a[4];
    const double t8 = x[11];
    const double t30 = t29*t8;
    const double t35 = t29*t7;
    const double t36 = t27*t8;
    const double t41 = a[15];
    const double t43 = a[18];
    const double t11 = x[7];
    const double t44 = t43*t11;
    const double t12 = x[8];
    const double t45 = t43*t12;
    const double t46 = a[7];
    const double t13 = x[9];
    const double t47 = t46*t13;
    const double t48 = a[12];
    const double t49 = t48*t7;
    const double t50 = t48*t8;
    const double t51 = a[23];
    const double t15 = x[12];
    const double t52 = t51*t15;
    const double t53 = t51*t6;
    const double t54 = a[9];
    const double t55 = t54*t2;
    const double t56 = a[1];
    const double t60 = a[19];
    const double t62 = t51*t7;
    const double t63 = t51*t8;
    const double t64 = t48*t15;
    const double t65 = t48*t6;
    const double t16 = x[5];
    const double t18 = x[6];
    const double t66 = t41*t16+t60*t18+t44+t45+t47+t55+t56+t62+t63+t64+t65;
    const double t69 = a[22];
    const double t70 = t69*t16;
    const double t71 = t69*t18;
    const double t72 = a[16];
    const double t74 = a[5];
    const double t76 = t54*t13;
    const double t77 = t46*t2;
    const double t20 = x[4];
    const double t78 = t72*t11+t74*t12+t41*t20+t49+t53+t56+t63+t64+t70+t71+t76+t77;
    const double t24 = x[3];
    const double t84 = t74*t11+t72*t12+t60*t20+t41*t24+t50+t52+t56+t62+t65+t70+t71+t76+t77;
    const double t86 = t15+t6;
    const double t90 = t43*t20;
    const double t91 = t43*t24;
    const double t101 = a[14];
    const double t104 = a[8];
    const double t111 = a[21];
    const double t113 = a[20];
    const double t33 = x[1];
    const double t37 = x[2];
    const double t81 = x[0];
    const double t120 = t101*t11+t101*t12+t101*t33+t101*t37+t104*t16+t104*t18+t104*t20+t104*
t24+t111*t13+t111*t2+t113*t15+t113*t6+t113*t7+t113*t8+a[11]*t81+a[2];
    return((t1*t6+t4+t5)*t6+(t1*t15+t10+t4+t5)*t15+(t1*t8+t14*t15+t10+t4+t5)*t8+(t1*t7+t14*t6+t9*t15+t9*t8+t4+t5)*t7+t3*(t6+t8+t7+t15)*t13+(t27*t15+t29*t6+t28
+t30)*t12+(t29*t15+t27*t6+t35+t36)*t11+(t41*t18+t44+t45+t47+t49+t50+t52+t53+t55
+t56)*t18+t66*t16+t78*t20+t84*t24+(t72*t16+t74*t18+t29*t86+t28+t36+t90+t91)*t37
+(t74*t16+t72*t18+t27*t86+t30+t35+t90+t91)*t33+t120*t81);

}

} // namespace mb_system