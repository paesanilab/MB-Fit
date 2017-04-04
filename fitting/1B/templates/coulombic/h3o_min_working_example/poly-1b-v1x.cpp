#include "poly-1b-v1x.h"

namespace h3o {

double poly_1b_v1x::eval(const double a[50], const double x[6],
                        double g[6])
{
    const double t1 = a[1];
    const double t2 = a[7];
    const double t3 = a[23];
    const double t6 = x[5];
    const double t4 = t3*t6;
    const double t5 = a[16];
    const double t7 = (t4+t5)*t6;
    const double t9 = (t2+t7)*t6;
    const double t12 = a[2];
    const double t13 = a[35];
    const double t14 = t13*t6;
    const double t15 = a[9];
    const double t17 = (t14+t15)*t6;
    const double t19 = (t12+t17)*t6;
    const double t20 = a[28];
    const double t21 = t20*t6;
    const double t23 = (t21+t15)*t6;
    const double t22 = x[4];
    const double t24 = t3*t22;
    const double t26 = (t24+t14+t5)*t22;
    const double t28 = (t2+t23+t26)*t22;
    const double t31 = a[37];
    const double t32 = t31*t6;
    const double t33 = a[20];
    const double t35 = (t32+t33)*t6;
    const double t36 = t13*t22;
    const double t38 = (t36+t32+t15)*t22;
    const double t40 = (t12+t35+t38)*t22;
    const double t41 = t20*t22;
    const double t43 = (t41+t32+t15)*t22;
    const double t39 = x[3];
    const double t44 = t3*t39;
    const double t46 = (t44+t36+t14+t5)*t39;
    const double t48 = (t2+t23+t43+t46)*t39;
    const double t51 = a[0];
    const double t52 = a[3];
    const double t53 = a[32];
    const double t54 = t53*t6;
    const double t55 = a[12];
    const double t57 = (t54+t55)*t6;
    const double t59 = (t52+t57)*t6;
    const double t60 = a[41];
    const double t61 = t60*t6;
    const double t62 = a[15];
    const double t64 = (t61+t62)*t6;
    const double t65 = t53*t22;
    const double t67 = (t65+t61+t55)*t22;
    const double t69 = (t52+t64+t67)*t22;
    const double t70 = a[4];
    const double t71 = a[38];
    const double t72 = t71*t6;
    const double t73 = a[14];
    const double t75 = (t72+t73)*t6;
    const double t76 = t71*t22;
    const double t77 = a[34];
    const double t78 = t77*t6;
    const double t80 = (t76+t78+t73)*t22;
    const double t81 = a[22];
    const double t82 = t81*t39;
    const double t83 = a[36];
    const double t84 = t83*t22;
    const double t85 = t83*t6;
    const double t86 = a[19];
    const double t88 = (t82+t84+t85+t86)*t39;
    const double t90 = (t70+t75+t80+t88)*t39;
    const double t91 = a[5];
    const double t92 = a[43];
    const double t93 = t92*t6;
    const double t94 = a[17];
    const double t96 = (t93+t94)*t6;
    const double t97 = t92*t22;
    const double t98 = a[44];
    const double t99 = t98*t6;
    const double t101 = (t97+t99+t94)*t22;
    const double t102 = a[48];
    const double t103 = t102*t39;
    const double t104 = a[46];
    const double t105 = t104*t22;
    const double t106 = t104*t6;
    const double t107 = a[13];
    const double t109 = (t103+t105+t106+t107)*t39;
    const double t110 = a[27];
    const double t95 = x[2];
    const double t111 = t110*t95;
    const double t112 = a[24];
    const double t113 = t112*t39;
    const double t114 = a[45];
    const double t115 = t114*t22;
    const double t116 = t114*t6;
    const double t117 = a[8];
    const double t119 = (t111+t113+t115+t116+t117)*t95;
    const double t121 = (t91+t96+t101+t109+t119)*t95;
    const double t124 = t81*t22;
    const double t126 = (t124+t85+t86)*t22;
    const double t128 = (t70+t75+t126)*t22;
    const double t130 = (t84+t78+t73)*t22;
    const double t131 = t53*t39;
    const double t133 = (t131+t76+t61+t55)*t39;
    const double t135 = (t52+t64+t130+t133)*t39;
    const double t136 = a[6];
    const double t137 = a[49];
    const double t138 = t137*t6;
    const double t139 = a[18];
    const double t141 = (t138+t139)*t6;
    const double t142 = a[30];
    const double t143 = t142*t22;
    const double t144 = a[31];
    const double t145 = t144*t6;
    const double t146 = a[10];
    const double t148 = (t143+t145+t146)*t22;
    const double t149 = t142*t39;
    const double t150 = a[25];
    const double t151 = t150*t22;
    const double t153 = (t149+t151+t145+t146)*t39;
    const double t154 = a[33];
    const double t155 = t154*t95;
    const double t156 = a[29];
    const double t157 = t156*t39;
    const double t158 = a[39];
    const double t159 = t158*t22;
    const double t160 = a[42];
    const double t161 = t160*t6;
    const double t162 = a[11];
    const double t164 = (t155+t157+t159+t161+t162)*t95;
    const double t166 = (t136+t141+t148+t153+t164)*t95;
    const double t167 = t102*t22;
    const double t169 = (t167+t106+t107)*t22;
    const double t170 = t92*t39;
    const double t172 = (t170+t105+t99+t94)*t39;
    const double t173 = a[47];
    const double t174 = t173*t95;
    const double t175 = t158*t39;
    const double t176 = t156*t22;
    const double t178 = (t174+t175+t176+t161+t162)*t95;
    const double t165 = x[1];
    const double t179 = t110*t165;
    const double t180 = t114*t39;
    const double t181 = t112*t22;
    const double t183 = (t179+t155+t180+t181+t116+t117)*t165;
    const double t185 = (t91+t96+t169+t172+t178+t183)*t165;
    const double t188 = t81*t6;
    const double t190 = (t188+t86)*t6;
    const double t192 = (t70+t190)*t6;
    const double t194 = (t85+t73)*t6;
    const double t196 = (t65+t72+t55)*t22;
    const double t198 = (t52+t194+t196)*t22;
    const double t199 = t60*t22;
    const double t201 = (t199+t78+t62)*t22;
    const double t203 = (t131+t199+t72+t55)*t39;
    const double t205 = (t52+t194+t201+t203)*t39;
    const double t206 = t142*t6;
    const double t208 = (t206+t146)*t6;
    const double t209 = t137*t22;
    const double t211 = (t209+t145+t139)*t22;
    const double t212 = t144*t22;
    const double t213 = t150*t6;
    const double t215 = (t149+t212+t213+t146)*t39;
    const double t216 = t160*t22;
    const double t217 = t158*t6;
    const double t219 = (t155+t157+t216+t217+t162)*t95;
    const double t221 = (t136+t208+t211+t215+t219)*t95;
    const double t223 = (t143+t213+t146)*t22;
    const double t224 = t137*t39;
    const double t226 = (t224+t212+t145+t139)*t39;
    const double t227 = a[40];
    const double t228 = t227*t95;
    const double t229 = a[26];
    const double t230 = t229*t39;
    const double t231 = t229*t22;
    const double t232 = t229*t6;
    const double t233 = a[21];
    const double t235 = (t228+t230+t231+t232+t233)*t95;
    const double t236 = t154*t165;
    const double t237 = t160*t39;
    const double t239 = (t236+t228+t237+t176+t217+t162)*t165;
    const double t241 = (t136+t208+t223+t226+t235+t239)*t165;
    const double t242 = t102*t6;
    const double t244 = (t242+t107)*t6;
    const double t246 = (t97+t106+t94)*t22;
    const double t247 = t98*t22;
    const double t249 = (t170+t247+t106+t94)*t39;
    const double t250 = t156*t6;
    const double t252 = (t174+t175+t216+t250+t162)*t95;
    const double t253 = t173*t165;
    const double t255 = (t253+t228+t237+t159+t250+t162)*t165;
    const double t243 = x[0];
    const double t256 = t110*t243;
    const double t257 = t112*t6;
    const double t259 = (t256+t236+t155+t180+t115+t257+t117)*t243;
    const double t261 = (t91+t244+t246+t249+t252+t255+t259)*t243;
    const double t279 = t154*t243;
    const double t291 = 2.0*t155;
    const double t294 = 2.0*t174;
    const double t301 = t227*t165;
    const double t324 = 2.0*t131;
    const double t327 = t156*t95;
    const double t328 = 2.0*t149;
    const double t331 = t114*t165;
    const double t332 = t158*t95;
    const double t333 = 2.0*t170;
    const double t342 = t160*t165;
    const double t343 = t229*t95;
    const double t347 = t114*t243;
    const double t361 = t13*t39;
    const double t367 = 2.0*t65;
    const double t370 = t83*t39;
    const double t374 = t114*t95;
    const double t375 = t104*t39;
    const double t376 = 2.0*t97;
    const double t384 = t71*t39;
    const double t388 = t150*t39;
    const double t389 = 2.0*t143;
    const double t400 = t60*t39;
    const double t404 = t160*t95;
    const double t405 = t144*t39;
    const double t409 = t156*t165;
    const double t412 = t158*t165;
    const double t413 = t98*t39;
    const double t426 = (2.0*t14+t15)*t6;
    const double t427 = 2.0*t21;
    const double t432 = t31*t22;
    const double t442 = (2.0*t54+t55)*t6;
    const double t443 = 2.0*t61;
    const double t446 = t77*t22;
    const double t447 = 2.0*t72;
    const double t450 = 2.0*t93;
    const double t469 = 2.0*t85;
    const double t474 = 2.0*t206;
    g[0] = ((2.0*t256+t236+t155+t180+t115+t257+t117)*t243+t91+t244+t246+
t249+t252+t255+t259)*t243+t51+t192+t198+t205+t221+t241+t261;
    g[1] = ((2.0*t179+t155+t180+t181+t116+t117)*t165+t91+t96+t169+t172+
t178+t183)*t165+t51+t59+t128+t135+t166+t185+((2.0*t236+t228+t237+t176+t217+t162
)*t165+t136+t208+t223+t226+t235+t239+(t279+2.0*t253+t228+t237+t159+t250+t162)*
t243)*t243;
    g[2] = ((2.0*t111+t113+t115+t116+t117)*t95+t91+t96+t101+t109+t119)*
t95+t51+t59+t69+t90+t121+((t291+t157+t159+t161+t162)*t95+t136+t141+t148+t153+
t164+(t236+t294+t175+t176+t161+t162)*t165)*t165+((t291+t157+t216+t217+t162)*t95
+t136+t208+t211+t215+t219+(t301+2.0*t228+t230+t231+t232+t233)*t165+(t279+t301+
t294+t175+t216+t250+t162)*t243)*t243;
    g[3] = ((2.0*t44+t36+t14+t5)*t39+t2+t23+t43+t46)*t39+t1+t19+t40+t48+
((2.0*t82+t84+t85+t86)*t39+t70+t75+t80+t88+(t112*t95+2.0*t103+t105+t106+t107)*
t95)*t95+((t324+t76+t61+t55)*t39+t52+t64+t130+t133+(t327+t328+t151+t145+t146)*
t95+(t331+t332+t333+t105+t99+t94)*t165)*t165+((t324+t199+t72+t55)*t39+t52+t194+
t201+t203+(t327+t328+t212+t213+t146)*t95+(t342+t343+2.0*t224+t212+t145+t139)*
t165+(t347+t342+t332+t333+t247+t106+t94)*t243)*t243;
    g[4] = ((2.0*t24+t14+t5)*t22+t2+t23+t26)*t22+t1+t19+t28+((2.0*t36+
t32+t15)*t22+t12+t35+t38+(t361+2.0*t41+t32+t15)*t39)*t39+((t367+t61+t55)*t22+
t52+t64+t67+(t370+2.0*t76+t78+t73)*t39+(t374+t375+t376+t99+t94)*t95)*t95+((2.0*
t124+t85+t86)*t22+t70+t75+t126+(t384+2.0*t84+t78+t73)*t39+(t332+t388+t389+t145+
t146)*t95+(t112*t165+t106+t107+2.0*t167+t327+t375)*t165)*t165+((t367+t72+t55)*
t22+t52+t194+t196+(t400+2.0*t199+t78+t62)*t39+(t404+t405+2.0*t209+t145+t139)*
t95+(t409+t343+t405+t389+t213+t146)*t165+(t347+t412+t404+t413+t376+t106+t94)*
t243)*t243;
    g[5] = ((2.0*t4+t5)*t6+t2+t7)*t6+t1+t9+(t426+t12+t17+(t36+t427+t15)*
t22)*t22+(t426+t12+t17+(t432+2.0*t32+t33)*t22+(t361+t432+t427+t15)*t39)*t39+(
t442+t52+t57+(t199+t443+t62)*t22+(t370+t446+t447+t73)*t39+(t374+t375+t247+t450+
t94)*t95)*t95+(t442+t52+t57+(t84+t447+t73)*t22+(t400+t446+t443+t62)*t39+(t404+
t405+t212+2.0*t138+t139)*t95+(t331+t404+t413+t105+t450+t94)*t165)*t165+((2.0*
t188+t86)*t6+t70+t190+(t76+t469+t73)*t22+(t384+t446+t469+t73)*t39+(t332+t388+
t212+t474+t146)*t95+(t412+t343+t405+t151+t474+t146)*t165+(t112*t243+t105+t107+
2.0*t242+t327+t375+t409)*t243)*t243;
    return (t1+t9)*t6+(t1+t19+t28)*t22+(t1+t19+t40+t48)*t39+(t51+t59+t69
+t90+t121)*t95+(t51+t59+t128+t135+t166+t185)*t165+(t51+t192+t198+t205+t221+t241
+t261)*t243;

}

} // namespace mb_system
