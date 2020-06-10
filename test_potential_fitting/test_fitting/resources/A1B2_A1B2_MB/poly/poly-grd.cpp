#include "poly-model.h"

namespace mb_system {

double poly_model::eval(const double a[24], const double x[15],
                        double g[15])
{
    const double t1 = a[13];
    const double t6 = x[13];
    const double t2 = t6*t1;
    const double t3 = a[6];
    const double t7 = x[14];
    const double t4 = t3*t7;
    const double t5 = a[0];
    const double t11 = x[12];
    const double t8 = t11*t1;
    const double t9 = a[10];
    const double t10 = t6*t9;
    const double t12 = x[11];
    const double t13 = t12*t1;
    const double t14 = a[17];
    const double t15 = t11*t14;
    const double t16 = x[10];
    const double t18 = t1*t16;
    const double t19 = t9*t12;
    const double t20 = t9*t11;
    const double t21 = t14*t6;
    const double t25 = t3*(t6+t12+t16+t11);
    const double t27 = a[3];
    const double t28 = t27*t16;
    const double t29 = a[4];
    const double t30 = t29*t12;
    const double t31 = t27*t11;
    const double t32 = t29*t6;
    const double t35 = t29*t16;
    const double t36 = t27*t12;
    const double t37 = t29*t11;
    const double t38 = t27*t6;
    const double t41 = a[15];
    const double t22 = x[6];
    const double t42 = t22*t41;
    const double t43 = a[18];
    const double t23 = x[7];
    const double t44 = t23*t43;
    const double t24 = x[8];
    const double t45 = t43*t24;
    const double t46 = a[7];
    const double t26 = x[9];
    const double t47 = t46*t26;
    const double t48 = a[12];
    const double t49 = t16*t48;
    const double t50 = t48*t12;
    const double t51 = a[23];
    const double t52 = t51*t11;
    const double t53 = t6*t51;
    const double t54 = a[9];
    const double t55 = t54*t7;
    const double t56 = a[1];
    const double t33 = x[5];
    const double t59 = t33*t41;
    const double t60 = a[19];
    const double t61 = t60*t22;
    const double t62 = t16*t51;
    const double t63 = t51*t12;
    const double t64 = t48*t11;
    const double t65 = t6*t48;
    const double t66 = t59+t61+t44+t45+t47+t62+t63+t64+t65+t55+t56;
    const double t34 = x[4];
    const double t68 = t34*t41;
    const double t69 = a[22];
    const double t70 = t33*t69;
    const double t71 = t69*t22;
    const double t72 = a[16];
    const double t73 = t72*t23;
    const double t74 = a[5];
    const double t75 = t74*t24;
    const double t76 = t54*t26;
    const double t77 = t46*t7;
    const double t78 = t68+t70+t71+t73+t75+t76+t49+t63+t64+t53+t77+t56;
    const double t39 = x[3];
    const double t80 = t39*t41;
    const double t81 = t60*t34;
    const double t82 = t74*t23;
    const double t83 = t72*t24;
    const double t84 = t80+t81+t70+t71+t82+t83+t76+t62+t50+t52+t65+t77+t56;
    const double t86 = t11+t6;
    const double t87 = t29*t86;
    const double t88 = t74*t22;
    const double t89 = t72*t33;
    const double t90 = t43*t34;
    const double t91 = t43*t39;
    const double t94 = t27*t86;
    const double t95 = t72*t22;
    const double t96 = t74*t33;
    const double t40 = x[0];
    const double t100 = t40*a[11];
    const double t101 = a[14];
    const double t58 = x[1];
    const double t102 = t58*t101;
    const double t67 = x[2];
    const double t103 = t67*t101;
    const double t104 = a[8];
    const double t105 = t39*t104;
    const double t106 = t34*t104;
    const double t107 = t33*t104;
    const double t108 = t104*t22;
    const double t109 = t23*t101;
    const double t110 = t101*t24;
    const double t111 = a[21];
    const double t112 = t26*t111;
    const double t113 = a[20];
    const double t114 = t16*t113;
    const double t115 = t12*t113;
    const double t116 = t11*t113;
    const double t117 = t6*t113;
    const double t118 = t111*t7;
    const double t119 = a[2];
    const double t120 = t100+t102+t103+t105+t106+t107+t108+t109+t110+t112+t114+t115+t116+
t117+t118+t119;
    const double t122 = (t2+t4+t5)*t6+(t8+t10+t4+t5)*t11+(t13+t15+t10+t4+t5)*t12+(t18+t19+
t20+t21+t4+t5)*t16+t25*t26+(t28+t30+t31+t32)*t24+(t35+t36+t37+t38)*t23+(t42+t44
+t45+t47+t49+t50+t52+t53+t55+t56)*t22+t66*t33+t78*t34+t84*t39+(t87+t36+t28+t88+
t89+t90+t91)*t67+(t30+t94+t35+t95+t96+t90+t91)*t58+t120*t40;
    const double t124 = 2.0*t100+t102+t103+t105+t106+t107+t108+t109+t110+t112+t114+t115+t116
+t117+t118+t119;
    const double t125 = t101*t40;
    const double t128 = t40*t104;
    const double t129 = t58*t43;
    const double t130 = t67*t43;
    const double t132 = t128+t129+t130+2.0*t80+t81+t70+t71+t82+t83+t76+t62+t50+t52+t65+t77+
t56;
    const double t135 = t60*t39+t128+t129+t130+t49+t53+t56+t63+t64+2.0*t68+t70+t71+t73+t75+
t76+t77;
    const double t138 = t39*t69;
    const double t139 = t34*t69;
    const double t141 = t74*t58+t72*t67+t128+t138+t139+t44+t45+t47+t55+t56+2.0*t59+t61+t62+
t63+t64+t65;
    const double t146 = t60*t33+t72*t58+t74*t67+t128+t138+t139+2.0*t42+t44+t45+t47+t49+t50+
t52+t53+t55+t56;
    const double t149 = t33*t43;
    const double t150 = t22*t43;
    const double t159 = t111*t40;
    const double t161 = t40*t113;
    const double t162 = t58*t29;
    const double t163 = t67*t27;
    const double t164 = t39*t51;
    const double t165 = t34*t48;
    const double t166 = t33*t51;
    const double t167 = t22*t48;
    const double t168 = t23*t29;
    const double t169 = t24*t27;
    const double t170 = t26*t3;
    const double t172 = t161+t162+t163+t164+t165+t166+t167+t168+t169+t170+2.0*t18+t19+t20+
t21+t4+t5;
    const double t173 = t39*t48;
    const double t174 = t34*t51;
    const double t175 = t23*t27;
    const double t176 = t24*t29;
    const double t177 = t16*t9;
    const double t179 = t161+t162+t163+t173+t174+t166+t167+t175+t176+t170+t177+2.0*t13+t15+
t10+t4+t5;
    const double t180 = t58*t27;
    const double t181 = t67*t29;
    const double t182 = t33*t48;
    const double t183 = t22*t51;
    const double t186 = t14*t12+t10+t161+t164+t165+t168+t169+t170+t177+t180+t181+t182+t183+
t4+t5+2.0*t8;
    const double t189 = t14*t16+t161+t170+t173+t174+t175+t176+t180+t181+t182+t183+t19+2.0*t2
+t20+t4+t5;
    g[0] = t124;
    g[1] = t30+t94+t35+t95+t96+t90+t91+t125;
    g[2] = t87+t36+t28+t88+t89+t90+t91+t125;
    g[3] = t132;
    g[4] = t135;
    g[5] = t141;
    g[6] = t146;
    g[7] = t72*t34+t74*t39+t125+t149+t150+t35+t36+t37+t38;
    g[8] = t74*t34+t72*t39+t125+t149+t150+t28+t30+t31+t32;
    g[9] = t46*t22+t46*t33+t54*t34+t54*t39+t159+t25;
    g[10] = t172;
    g[11] = t179;
    g[12] = t186;
    g[13] = t189;
    g[14] = t3*t11+t3*t12+t3*t16+t54*t22+t3*t6+t54*t33+t46*t34+t46*t39+
t159;
    return t122;

}

} // namespace mb_system