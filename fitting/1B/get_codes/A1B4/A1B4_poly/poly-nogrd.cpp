#include "poly-model.h"

namespace mb_system {

double poly_model::eval(const double a[82], const double x[10])
{
    const double t1 = a[0];
    const double t2 = a[7];
    const double t3 = a[45];
    const double t5 = a[26];
    const double t12 = a[2];
    const double t13 = a[31];
    const double t4 = x[9];
    const double t14 = t13*t4;
    const double t15 = a[20];
    const double t19 = (t12+(t14+t15)*t4)*t4;
    const double t20 = a[29];
    const double t23 = (t20*t4+t15)*t4;
    const double t31 = a[47];
    const double t32 = t31*t4;
    const double t33 = a[11];
    const double t35 = (t32+t33)*t4;
    const double t16 = x[8];
    const double t36 = t13*t16;
    const double t41 = t20*t16;
    const double t51 = a[34];
    const double t52 = t51*t4;
    const double t53 = a[21];
    const double t55 = (t52+t53)*t4;
    const double t60 = a[5];
    const double t61 = a[51];
    const double t62 = t61*t4;
    const double t63 = a[13];
    const double t65 = (t62+t63)*t4;
    const double t66 = t61*t16;
    const double t67 = a[30];
    const double t68 = t67*t4;
    const double t71 = a[66];
    const double t21 = x[7];
    const double t72 = t71*t21;
    const double t73 = a[57];
    const double t74 = t73*t16;
    const double t75 = t73*t4;
    const double t76 = a[27];
    const double t83 = a[79];
    const double t94 = t71*t16;
    const double t100 = (t74+t68+t63)*t16;
    const double t101 = t13*t21;
    const double t106 = t73*t21;
    const double t107 = a[56];
    const double t24 = x[6];
    const double t111 = t13*t24;
    const double t119 = t20*t21;
    const double t122 = t20*t24;
    const double t123 = t61*t21;
    const double t133 = t71*t4;
    const double t139 = (t75+t63)*t4;
    const double t144 = t51*t16;
    const double t151 = t31*t16;
    const double t154 = t67*t16;
    const double t155 = t107*t4;
    const double t164 = t31*t21;
    const double t167 = t51*t24;
    const double t26 = x[5];
    const double t171 = t13*t26;
    const double t195 = a[1];
    const double t196 = a[4];
    const double t197 = a[54];
    const double t199 = a[24];
    const double t203 = (t196+(t197*t4+t199)*t4)*t4;
    const double t204 = a[36];
    const double t205 = t204*t4;
    const double t206 = a[18];
    const double t208 = (t205+t206)*t4;
    const double t209 = t197*t16;
    const double t214 = a[3];
    const double t215 = a[35];
    const double t216 = t215*t4;
    const double t217 = a[12];
    const double t219 = (t216+t217)*t4;
    const double t220 = t215*t16;
    const double t221 = a[77];
    const double t222 = t221*t4;
    const double t225 = a[65];
    const double t226 = t225*t21;
    const double t227 = a[61];
    const double t228 = t227*t16;
    const double t229 = t227*t4;
    const double t230 = a[14];
    const double t235 = t204*t16;
    const double t236 = a[81];
    const double t237 = t236*t4;
    const double t240 = a[53];
    const double t241 = t240*t21;
    const double t242 = a[68];
    const double t243 = t242*t16;
    const double t244 = t242*t4;
    const double t245 = a[17];
    const double t248 = t197*t24;
    const double t249 = a[48];
    const double t250 = t249*t21;
    const double t255 = t249*t16;
    const double t258 = a[32];
    const double t259 = t258*t21;
    const double t260 = a[62];
    const double t261 = t260*t16;
    const double t262 = a[33];
    const double t263 = t262*t4;
    const double t264 = a[9];
    const double t267 = t215*t24;
    const double t268 = t260*t21;
    const double t271 = t225*t26;
    const double t272 = t227*t24;
    const double t273 = t240*t16;
    const double t278 = t249*t4;
    const double t280 = (t278+t245)*t4;
    const double t283 = t262*t16;
    const double t284 = t260*t4;
    const double t287 = t221*t16;
    const double t290 = t258*t26;
    const double t291 = t262*t24;
    const double t292 = a[71];
    const double t39 = x[4];
    const double t296 = t225*t39;
    const double t297 = t240*t4;
    const double t302 = a[6];
    const double t303 = a[76];
    const double t305 = a[16];
    const double t307 = (t303*t4+t305)*t4;
    const double t308 = t303*t16;
    const double t309 = a[41];
    const double t310 = t309*t4;
    const double t313 = a[46];
    const double t314 = t313*t21;
    const double t315 = a[40];
    const double t316 = t315*t16;
    const double t317 = t315*t4;
    const double t318 = a[22];
    const double t321 = t303*t24;
    const double t322 = a[44];
    const double t323 = t322*t21;
    const double t324 = t309*t16;
    const double t327 = t313*t26;
    const double t328 = t315*t24;
    const double t329 = a[50];
    const double t330 = t329*t21;
    const double t331 = t322*t16;
    const double t334 = t313*t39;
    const double t336 = t322*t4;
    const double t339 = a[73];
    const double t341 = a[38];
    const double t342 = t341*t39;
    const double t343 = t341*t26;
    const double t344 = a[74];
    const double t345 = t344*t24;
    const double t346 = t341*t21;
    const double t347 = t344*t16;
    const double t348 = t344*t4;
    const double t349 = a[28];
    const double t356 = t225*t16;
    const double t363 = t197*t21;
    const double t368 = t258*t16;
    const double t373 = t225*t24;
    const double t380 = t204*t21;
    const double t383 = t242*t21;
    const double t386 = t197*t26;
    const double t393 = t215*t21;
    const double t396 = t258*t24;
    const double t400 = t215*t26;
    const double t401 = t221*t21;
    const double t404 = t227*t26;
    const double t405 = t227*t21;
    const double t410 = a[8];
    const double t411 = a[58];
    const double t413 = a[10];
    const double t416 = a[49];
    const double t417 = t416*t16;
    const double t418 = a[39];
    const double t419 = t418*t4;
    const double t420 = a[23];
    const double t423 = t416*t21;
    const double t424 = a[70];
    const double t425 = t424*t16;
    const double t428 = t416*t24;
    const double t429 = a[42];
    const double t430 = t429*t21;
    const double t431 = a[43];
    const double t432 = t431*t16;
    const double t435 = t416*t26;
    const double t436 = t424*t24;
    const double t437 = t431*t21;
    const double t438 = t429*t16;
    const double t441 = a[52];
    const double t443 = a[64];
    const double t444 = t443*t26;
    const double t445 = t443*t24;
    const double t446 = t443*t21;
    const double t447 = t443*t16;
    const double t448 = a[78];
    const double t449 = t448*t4;
    const double t450 = a[15];
    const double t453 = a[59];
    const double t43 = x[3];
    const double t454 = t453*t43;
    const double t455 = a[69];
    const double t456 = t455*t39;
    const double t457 = a[55];
    const double t458 = t457*t26;
    const double t459 = a[37];
    const double t460 = t459*t24;
    const double t461 = t457*t21;
    const double t462 = t459*t16;
    const double t463 = a[63];
    const double t464 = t463*t4;
    const double t465 = a[19];
    const double t470 = t313*t16;
    const double t473 = t303*t21;
    const double t476 = t313*t24;
    const double t477 = t329*t16;
    const double t480 = t303*t26;
    const double t481 = t309*t21;
    const double t484 = t315*t26;
    const double t485 = t329*t24;
    const double t486 = t315*t21;
    const double t489 = a[67];
    const double t490 = t489*t43;
    const double t491 = t459*t26;
    const double t492 = t457*t24;
    const double t493 = t459*t21;
    const double t494 = t457*t16;
    const double t498 = t344*t26;
    const double t499 = t341*t24;
    const double t500 = t344*t21;
    const double t501 = t341*t16;
    const double t512 = (t214+(t225*t4+t230)*t4)*t4;
    const double t514 = (t229+t217)*t4;
    const double t525 = t258*t4;
    const double t527 = (t525+t264)*t4;
    const double t540 = t292*t4;
    const double t548 = (t297+t245)*t4;
    const double t558 = t197*t39;
    const double t565 = (t4*t416+t420)*t4;
    const double t569 = t418*t16;
    const double t570 = t424*t4;
    const double t573 = t431*t4;
    const double t577 = t448*t16;
    const double t578 = t443*t4;
    const double t581 = t416*t39;
    const double t582 = t429*t4;
    const double t585 = t457*t39;
    const double t586 = t455*t26;
    const double t587 = t463*t16;
    const double t588 = t459*t4;
    const double t599 = t448*t21;
    const double t602 = t418*t21;
    const double t608 = a[72];
    const double t609 = t608*t43;
    const double t610 = a[60];
    const double t611 = t610*t39;
    const double t612 = t610*t26;
    const double t613 = t610*t24;
    const double t614 = a[75];
    const double t615 = t614*t21;
    const double t616 = t614*t16;
    const double t617 = t614*t4;
    const double t618 = a[25];
    const double t57 = x[2];
    const double t621 = t453*t57;
    const double t622 = t455*t24;
    const double t623 = t463*t21;
    const double t630 = (t313*t4+t318)*t4;
    const double t635 = t329*t4;
    const double t640 = t303*t39;
    const double t643 = t459*t39;
    const double t644 = t457*t4;
    const double t647 = t489*t57;
    const double t651 = t344*t39;
    const double t652 = t341*t4;
    const double t681 = t204*t24;
    const double t694 = t204*t26;
    const double t710 = t418*t24;
    const double t716 = t463*t24;
    const double t717 = t455*t21;
    const double t732 = t418*t26;
    const double t733 = t431*t24;
    const double t734 = t424*t21;
    const double t737 = t614*t26;
    const double t738 = t614*t24;
    const double t739 = t610*t21;
    const double t740 = t610*t16;
    const double t743 = t463*t26;
    const double t744 = t455*t16;
    const double t763 = t614*t39;
    const double t764 = t610*t4;
    const double t767 = t608*t57;
    const double t64 = x[1];
    const double t772 = t453*t64;
    const double t773 = t463*t39;
    const double t774 = t455*t4;
    const double t785 = t309*t24;
    const double t69 = x[0];
    const double t799 = t339*t69+t345+t346+t349+t454+t498+t501+t621+t651+t652+t772;
    const double t801 = t302+t630+(t470+t635+t318)*t16+(t314+t477+t635+t318)*t21+(t321+t323+
t316+t317+t305)*t24+(t480+t785+t486+t331+t317+t305)*t26+(t26*t309+t305+t316+
t336+t486+t640+t785)*t39+(t490+t643+t491+t716+t717+t494+t644+t465)*t43+(t647+
t609+t643+t743+t460+t461+t744+t644+t465)*t57+(t489*t64+t460+t461+t465+t491+t494
+t609+t767+t773+t774)*t64+t799*t69;
    const double t803 = t195+t512+(t214+t527+(t356+t525+t230)*t16)*t16+(t214+t527+(t368+t540
+t264)*t16+(t226+t368+t525+t230)*t21)*t21+(t196+t514+(t228+t263+t217)*t16+(t241
+t261+t284+t245)*t21+(t248+t250+t220+t216+t199)*t24)*t24+(t196+t514+(t273+t284+
t245)*t16+(t405+t261+t263+t217)*t21+(t681+t383+t243+t222+t206)*t24+(t386+t681+
t393+t255+t216+t199)*t26)*t26+(t196+t548+(t228+t284+t217)*t16+(t405+t283+t284+
t217)*t21+(t681+t383+t287+t244+t206)*t24+(t236*t24+t206+t243+t244+t401+t694)*
t26+(t558+t694+t681+t393+t220+t278+t199)*t39)*t39+(t410+t565+(t417+t573+t420)*
t16+(t21*t441+t447+t450+t578)*t21+(t24*t411+t413+t419+t569+t599)*t24+(t435+t710
+t446+t438+t570+t420)*t26+(t26*t431+t420+t425+t446+t581+t582+t710)*t39+(t454+
t585+t458+t716+t717+t462+t588+t465)*t43)*t43+(t410+t565+(t16*t441+t450+t578)*
t16+(t423+t447+t573+t420)*t21+(t428+t430+t447+t570+t420)*t24+(t26*t411+t413+
t419+t577+t602+t710)*t26+(t581+t732+t733+t734+t447+t582+t420)*t39+(t609+t611+
t737+t738+t739+t740+t617+t618)*t43+(t621+t609+t585+t743+t492+t493+t744+t588+
t465)*t57)*t57+(t410+(t4*t441+t450)*t4+(t417+t578+t420)*t16+(t423+t432+t578+
t420)*t21+(t428+t430+t425+t578+t420)*t24+(t435+t733+t734+t438+t578+t420)*t26+(
t39*t411+t413+t449+t569+t602+t710+t732)*t39+(t609+t763+t612+t738+t739+t616+t764
+t618)*t43+(t43*a[80]+t613+t615+t618+t737+t740+t763+t764+t767)*t57+(t772+t767+
t609+t773+t458+t492+t493+t462+t774+t465)*t64)*t64+t801*t69;
    return((t1+(t2+(t3*t4+t5)*t4)*t4)*t4+(t1+t19+(t2+t23+(t16*t3+t14+t5)*t16)*t16)*t16+(t1+t19+(t12+t35+(t36+t32+t15)*t16)*t16+(t2+t23+(t41+t32+t15)*t16+(t21
*t3+t14+t36+t5)*t21)*t21)*t21+(t1+t19+(t12+t55+(t36+t52+t15)*t16)*t16+(t60+t65+
(t66+t68+t63)*t16+(t72+t74+t75+t76)*t21)*t21+(t2+t23+(t41+t52+t15)*t16+(t21*t83
+t74+t75+t76)*t21+(t24*t3+t14+t36+t5+t72)*t24)*t24)*t24+(t1+t19+(t60+t65+(t94+
t75+t76)*t16)*t16+(t12+t55+t100+(t101+t66+t52+t15)*t21)*t21+(t12+t35+t100+(t107
*t16+t106+t63+t68)*t21+(t111+t106+t66+t32+t15)*t24)*t24+(t2+t23+(t16*t83+t75+
t76)*t16+(t119+t74+t52+t15)*t21+(t122+t123+t74+t32+t15)*t24+(t26*t3+t101+t111+
t14+t5+t94)*t26)*t26)*t26+(t1+(t60+(t133+t76)*t4)*t4+(t12+t139+(t36+t62+t15)*
t16)*t16+(t12+t139+(t144+t68+t53)*t16+(t101+t144+t62+t15)*t21)*t21+(t12+t139+(
t151+t68+t33)*t16+(t106+t154+t155+t63)*t21+(t111+t106+t151+t62+t15)*t24)*t24+(
t12+t139+(t74+t155+t63)*t16+(t164+t154+t68+t33)*t21+(t21*t67+t154+t167+t53+t68)
*t24+(t171+t167+t164+t74+t62+t15)*t26)*t26+(t2+(t4*t83+t76)*t4+(t41+t75+t15)*
t16+(t119+t144+t75+t15)*t21+(t122+t123+t151+t75+t15)*t24+(t20*t26+t15+t164+t167
+t66+t75)*t26+(t3*t39+t101+t111+t133+t171+t36+t5)*t39)*t39)*t39+(t195+t203+(
t196+t208+(t209+t205+t199)*t16)*t16+(t214+t219+(t220+t222+t217)*t16+(t226+t228+
t229+t230)*t21)*t21+(t196+t208+(t235+t237+t206)*t16+(t241+t243+t244+t245)*t21+(
t248+t250+t235+t205+t199)*t24)*t24+(t214+t219+(t255+t244+t245)*t16+(t259+t261+
t263+t264)*t21+(t267+t268+t243+t222+t217)*t24+(t271+t272+t259+t273+t229+t230)*
t26)*t26+(t214+t280+(t220+t244+t217)*t16+(t259+t283+t284+t264)*t21+(t267+t268+
t287+t244+t217)*t24+(t21*t292+t261+t264+t284+t290+t291)*t26+(t296+t290+t272+
t259+t228+t297+t230)*t39)*t39+(t302+t307+(t308+t310+t305)*t16+(t314+t316+t317+
t318)*t21+(t321+t323+t324+t310+t305)*t24+(t327+t328+t330+t331+t317+t318)*t26+(
t26*t329+t316+t318+t328+t330+t334+t336)*t39+(t339*t43+t342+t343+t345+t346+t347+
t348+t349)*t43)*t43)*t43+(t195+t203+(t214+t219+(t356+t229+t230)*t16)*t16+(t196+
t208+(t228+t222+t217)*t16+(t363+t220+t205+t199)*t21)*t21+(t214+t219+(t368+t263+
t264)*t16+(t250+t261+t244+t245)*t21+(t373+t241+t368+t229+t230)*t24)*t24+(t196+
t208+(t273+t244+t245)*t16+(t380+t243+t237+t206)*t21+(t272+t383+t261+t222+t217)*
t24+(t386+t267+t380+t255+t205+t199)*t26)*t26+(t214+t280+(t368+t284+t264)*t16+(
t393+t283+t244+t217)*t21+(t16*t292+t264+t268+t284+t396)*t24+(t400+t291+t401+
t261+t244+t217)*t26+(t296+t404+t396+t405+t368+t297+t230)*t39)*t39+(t410+(t4*
t411+t413)*t4+(t417+t419+t420)*t16+(t423+t425+t419+t420)*t21+(t428+t430+t432+
t419+t420)*t24+(t435+t436+t437+t438+t419+t420)*t26+(t39*t441+t444+t445+t446+
t447+t449+t450)*t39+(t454+t456+t458+t460+t461+t462+t464+t465)*t43)*t43+(t302+
t307+(t470+t317+t318)*t16+(t473+t316+t310+t305)*t21+(t476+t323+t477+t317+t318)*
t24+(t480+t328+t481+t331+t310+t305)*t26+(t334+t484+t485+t486+t477+t336+t318)*
t39+(t490+t456+t491+t492+t493+t494+t464+t465)*t43+(t339*t57+t342+t348+t349+t454
+t498+t499+t500+t501)*t57)*t57)*t57+(t195+t512+(t196+t514+(t209+t216+t199)*t16)
*t16+(t196+t514+(t235+t222+t206)*t16+(t363+t235+t216+t199)*t21)*t21+(t214+t527+
(t220+t263+t217)*t16+(t250+t243+t284+t245)*t21+(t373+t241+t228+t525+t230)*t24)*
t24+(t214+t527+(t255+t284+t245)*t16+(t393+t243+t263+t217)*t21+(t396+t268+t261+
t540+t264)*t24+(t271+t396+t405+t273+t525+t230)*t26)*t26+(t196+t548+(t235+t244+
t206)*t16+(t16*t236+t206+t244+t380)*t21+(t272+t383+t287+t284+t217)*t24+(t404+
t291+t401+t243+t284+t217)*t26+(t558+t400+t267+t380+t235+t278+t199)*t39)*t39+(
t410+t565+(t16*t411+t413+t419)*t16+(t423+t569+t570+t420)*t21+(t428+t430+t569+
t573+t420)*t24+(t26*t441+t445+t446+t450+t577+t578)*t26+(t581+t444+t436+t437+
t569+t582+t420)*t39+(t454+t585+t586+t460+t461+t587+t588+t465)*t43)*t43+(t410+
t565+(t417+t570+t420)*t16+(t21*t411+t413+t419+t569)*t21+(t24*t441+t447+t450+
t578+t599)*t24+(t435+t445+t602+t438+t573+t420)*t26+(t26*t424+t420+t432+t445+
t581+t582+t602)*t39+(t609+t611+t612+t613+t615+t616+t617+t618)*t43+(t621+t609+
t585+t491+t622+t623+t494+t588+t465)*t57)*t57+(t302+t630+(t308+t317+t305)*t16+(
t473+t324+t317+t305)*t21+(t476+t323+t316+t635+t318)*t24+(t327+t485+t486+t331+
t635+t318)*t26+(t640+t484+t328+t481+t324+t336+t305)*t39+(t490+t643+t586+t492+
t493+t587+t644+t465)*t43+(t647+t609+t643+t458+t622+t623+t462+t644+t465)*t57+(
t339*t64+t343+t347+t349+t454+t499+t500+t621+t651+t652)*t64)*t64)*t64+t803*t69);

}

} // namespace mb_system