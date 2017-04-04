#include "poly-1b-v1x.h"

namespace h3o {

double poly_1b_v1x::eval(const double a[50], const double x[6]){
    const double t1 = a[1];
    const double t2 = a[7];
    const double t3 = a[23];
    const double t5 = a[16];
    const double t12 = a[2];
    const double t13 = a[35];
    const double t4 = x[5];
    const double t14 = t13*t4;
    const double t15 = a[9];
    const double t19 = (t12+(t14+t15)*t4)*t4;
    const double t20 = a[28];
    const double t23 = (t20*t4+t15)*t4;
    const double t32 = a[37]*t4;
    const double t16 = x[4];
    const double t36 = t13*t16;
    const double t51 = a[0];
    const double t52 = a[3];
    const double t53 = a[32];
    const double t55 = a[12];
    const double t59 = (t52+(t53*t4+t55)*t4)*t4;
    const double t60 = a[41];
    const double t61 = t60*t4;
    const double t62 = a[15];
    const double t64 = (t61+t62)*t4;
    const double t65 = t53*t16;
    const double t70 = a[4];
    const double t71 = a[38];
    const double t72 = t71*t4;
    const double t73 = a[14];
    const double t75 = (t72+t73)*t4;
    const double t76 = t71*t16;
    const double t78 = a[34]*t4;
    const double t81 = a[22];
    const double t83 = a[36];
    const double t84 = t83*t16;
    const double t85 = t83*t4;
    const double t86 = a[19];
    const double t91 = a[5];
    const double t92 = a[43];
    const double t94 = a[17];
    const double t96 = (t92*t4+t94)*t4;
    const double t97 = t92*t16;
    const double t98 = a[44];
    const double t99 = t98*t4;
    const double t102 = a[48];
    const double t104 = a[46];
    const double t105 = t104*t16;
    const double t106 = t104*t4;
    const double t107 = a[13];
    const double t110 = a[27];
    const double t112 = a[24];
    const double t114 = a[45];
    const double t115 = t114*t16;
    const double t116 = t114*t4;
    const double t117 = a[8];
    const double t29 = x[3];
    const double t131 = t53*t29;
    const double t136 = a[6];
    const double t137 = a[49];
    const double t139 = a[18];
    const double t142 = a[30];
    const double t143 = t142*t16;
    const double t144 = a[31];
    const double t145 = t144*t4;
    const double t146 = a[10];
    const double t149 = t142*t29;
    const double t150 = a[25];
    const double t154 = a[33];
    const double t30 = x[2];
    const double t155 = t154*t30;
    const double t156 = a[29];
    const double t157 = t156*t29;
    const double t158 = a[39];
    const double t159 = t158*t16;
    const double t160 = a[42];
    const double t161 = t160*t4;
    const double t162 = a[11];
    const double t170 = t92*t29;
    const double t173 = a[47];
    const double t174 = t173*t30;
    const double t175 = t158*t29;
    const double t176 = t156*t16;
    const double t180 = t114*t29;
    const double t194 = (t85+t73)*t4;
    const double t199 = t60*t16;
    const double t208 = (t142*t4+t146)*t4;
    const double t212 = t144*t16;
    const double t213 = t150*t4;
    const double t216 = t160*t16;
    const double t217 = t158*t4;
    const double t228 = a[40]*t30;
    const double t229 = a[26];
    const double t37 = x[1];
    const double t236 = t154*t37;
    const double t237 = t160*t29;
    const double t250 = t156*t4;
    const double t261 = x[0];
    return((t1+(t2+(t3*t4+t5)*t4)*t4)*t4+(t1+t19+(t2+t23+(t3*t16+t14+t5)*t16)*t16)*t16+(t1+t19+(t12+(t32+a[20])*t4+(t36+t32+t15)*t16)*t16+(t2+t23+(t20*t16+
t15+t32)*t16+(t3*t29+t14+t36+t5)*t29)*t29)*t29+(t51+t59+(t52+t64+(t65+t61+t55)*
t16)*t16+(t70+t75+(t76+t78+t73)*t16+(t81*t29+t84+t85+t86)*t29)*t29+(t91+t96+(
t97+t99+t94)*t16+(t102*t29+t105+t106+t107)*t29+(t110*t30+t112*t29+t115+t116+
t117)*t30)*t30)*t30+(t51+t59+(t70+t75+(t81*t16+t85+t86)*t16)*t16+(t52+t64+(t84+
t78+t73)*t16+(t131+t76+t61+t55)*t29)*t29+(t136+(t137*t4+t139)*t4+(t143+t145+
t146)*t16+(t150*t16+t145+t146+t149)*t29+(t155+t157+t159+t161+t162)*t30)*t30+(
t91+t96+(t102*t16+t106+t107)*t16+(t170+t105+t99+t94)*t29+(t174+t175+t176+t161+
t162)*t30+(t110*t37+t112*t16+t116+t117+t155+t180)*t37)*t37)*t37+(t51+(t70+(t81*
t4+t86)*t4)*t4+(t52+t194+(t65+t72+t55)*t16)*t16+(t52+t194+(t199+t78+t62)*t16+(
t131+t199+t72+t55)*t29)*t29+(t136+t208+(t137*t16+t139+t145)*t16+(t149+t212+t213
+t146)*t29+(t155+t157+t216+t217+t162)*t30)*t30+(t136+t208+(t143+t213+t146)*t16+
(t137*t29+t139+t145+t212)*t29+(t229*t16+t229*t29+t229*t4+t228+a[21])*t30+(t236+
t228+t237+t176+t217+t162)*t37)*t37+(t91+(t102*t4+t107)*t4+(t97+t106+t94)*t16+(
t98*t16+t106+t170+t94)*t29+(t174+t175+t216+t250+t162)*t30+(t173*t37+t159+t162+
t228+t237+t250)*t37+(t110*t261+t112*t4+t115+t117+t155+t180+t236)*t261)*t261)*
t261);
}

} // namespace mb_system
