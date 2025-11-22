/*
 * delta2star3D_mex.c - MEX implementation of delta2star3D for maximum performance
 * 
 * This is a C implementation of the delta2star3D function for use as a MEX file.
 * Compilation: mex -O delta2star3D_mex.c -output delta2star3D_mex
 * 
 * Usage from MATLAB: E = delta2star3D_mex(s300, s400, ..., s022)
 */

#include "mex.h"
#include <math.h>

void delta2star3D_compute(const double *s, double *E) {
    /* Extract inputs for readability */
    double s300=s[0], s400=s[1], s110=s[2], s210=s[3], s310=s[4];
    double s120=s[5], s220=s[6], s030=s[7], s130=s[8], s040=s[9];
    double s101=s[10], s201=s[11], s301=s[12], s102=s[13], s202=s[14];
    double s003=s[15], s103=s[16], s004=s[17], s011=s[18], s111=s[19];
    double s211=s[20], s021=s[21], s121=s[22], s031=s[23], s012=s[24];
    double s112=s[25], s013=s[26], s022=s[27];
    
    /* Pre-compute powers */
    double t41 = s003*s003, t42 = s011*s011, t43 = t42*s011;
    double t44 = s012*s012, t45 = s021*s021, t46 = s030*s030;
    double t47 = s101*s101, t48 = t47*s101, t49 = s102*s102;
    double t50 = s110*s110, t51 = t50*s110, t52 = s111*s111;
    double t53 = s120*s120, t54 = s201*s201, t55 = s210*s210;
    double t56 = s300*s300;
    
    /* Level-1 products */
    double t2 = s003*s012, t3 = s003*s021, t4 = s012*s021;
    double t5 = s012*s030, t6 = s021*s030, t7 = s003*s102;
    double t8 = s011*s101, t9 = s003*s111, t10 = s012*s102;
    double t11 = s011*s110, t12 = s012*s111, t13 = s021*s102;
    double t14 = s012*s120, t15 = s021*s111, t16 = s021*s120;
    double t17 = s030*s111, t18 = s030*s120, t19 = s003*s201;
    double t20 = s101*s110, t21 = s012*s201, t22 = s102*s111;
    double t23 = s012*s210, t24 = s021*s201, t25 = s102*s120;
    double t26 = s021*s210, t27 = s111*s120, t28 = s030*s210;
    double t29 = s102*s201, t30 = s102*s210, t31 = s111*s201;
    double t32 = s111*s210, t33 = s120*s201, t34 = s120*s210;
    double t35 = s102*s300, t36 = s111*s300, t37 = s201*s210;
    double t38 = s120*s300, t39 = s201*s300, t40 = s210*s300;
    
    /* Special 3-way products */
    double t58 = s003*s011*s030, t63 = s003*s011*s120;
    double t68 = s011*s030*s102, t73 = s003*s011*s210;
    double t74 = s003*s101*s120, t85 = s011*s030*s201;
    double t87 = s030*s102*s110, t91 = s003*s101*s210;
    double t109 = s030*s110*s201, t112 = s003*s101*s300;
    double t114 = s012*s101*s300, t117 = s012*s110*s300;
    double t118 = s021*s101*s300, t121 = s021*s110*s300;
    double t124 = s030*s110*s300;
    
    /* Negations */
    double t132 = -s013, t133 = -s031, t134 = -s103;
    double t135 = -s112, t136 = -s121, t137 = -s130;
    double t138 = -s211, t139 = -s301, t140 = -s310;
    
    /* Level-2 products (with s011, s101, s110) */
    double t57 = s011*t3, t59 = s011*t4, t60 = s011*t5;
    double t61 = s011*t9, t62 = s011*t10, t64 = s011*t12;
    double t65 = s011*t13, t66 = s011*t14, t67 = s011*t15;
    double t69 = s011*t16, t70 = s011*t17, t71 = s101*t9;
    double t72 = s101*t10, t75 = s011*t21, t76 = s101*t13;
    double t77 = s011*t23, t78 = s011*t24, t79 = s011*t25;
    double t80 = s101*t14, t81 = s110*t12, t82 = s101*t15;
    double t83 = s110*t13, t84 = s011*t26, t86 = s110*t14;
    double t88 = s110*t16, t89 = s110*t17, t90 = s101*t19;
    double t92 = s101*t21, t93 = s101*t22, t94 = s011*t30;
    double t95 = s011*t31, t96 = s101*t23, t97 = s110*t21;
    double t98 = s101*t24, t99 = s101*t25, t100 = s110*t22;
    double t101 = s011*t32, t102 = s011*t33, t103 = s110*t23;
    double t104 = s101*t26, t105 = s110*t24, t106 = s101*t27;
    double t107 = s110*t25, t108 = s110*t26, t110 = s110*t27;
    double t111 = s110*t28, t113 = s101*t29, t115 = s101*t30;
    double t116 = s101*t31, t119 = s101*t33, t120 = s110*t30;
    double t122 = s110*t32, t123 = s110*t33, t125 = s110*t34;
    double t126 = s101*t35, t127 = s101*t36, t128 = s101*t37;
    double t129 = s110*t36, t130 = s110*t37, t131 = s110*t38;
    
    /* Triple products */
    double t141 = t9*t11, t142 = t3*t20, t143 = t8*t12;
    double t144 = t10*t11, t145 = t8*t13, t146 = s003*s120*t11;
    double t147 = s003*s030*t20, t148 = t8*t14, t149 = t11*t13;
    double t150 = s030*s102*t8, t151 = t4*t20, t152 = t11*t14;
    double t153 = t8*t16, t154 = t11*t15, t155 = t8*t17;
    double t156 = t5*t20, t157 = t11*t19, t158 = t9*t20;
    double t159 = t8*t21, t160 = t8*t22, t161 = t10*t20;
    double t162 = s003*s210*t11, t163 = s003*s120*t20, t164 = t8*t23;
    double t165 = t11*t21, t166 = t8*t24, t167 = t8*t25;
    double t168 = t11*t22, t169 = t12*t20, t170 = t13*t20;
    double t171 = t11*t23, t172 = t8*t26, t173 = t11*t24;
    double t174 = s030*s201*t8, t175 = t8*t27, t176 = t11*t25;
    double t177 = t14*t20, t178 = t15*t20, t179 = s030*s102*t20;
    double t180 = t11*t26, t181 = t8*t28, t182 = t11*t27;
    double t183 = t16*t20, t184 = t17*t20, t185 = s003*s300*t11;
    double t186 = s003*s210*t20, t187 = s012*s300*t8, t188 = t8*t30;
    double t189 = t11*t29, t190 = t20*t21, t191 = s012*s300*t11;
    double t192 = s021*s300*t8, t193 = t8*t32, t194 = t8*t33;
    double t195 = t11*t30, t196 = t11*t31, t197 = t20*t23;
    double t198 = t20*t24, t199 = t20*t25, t200 = s021*s300*t11;
    double t201 = s030*s300*t8, t202 = t8*t34, t203 = t11*t33;
    double t204 = t20*t26, t205 = s030*s201*t20, t206 = t8*t36;
    double t207 = t8*t37, t208 = t11*t35, t209 = t20*t30;
    double t210 = t20*t31, t211 = t8*t38, t212 = t11*t36;
    double t213 = t11*t37, t214 = t20*t32, t215 = t20*t33;
    
    /* Simple negations */
    double t216=-t3, t217=-t5, t218=-t8, t219=-t10, t220=-t11;
    double t221=-t12, t222=-t15, t223=-t16, t224=-t19, t225=-t20;
    double t226=-t22, t227=-t23, t228=-t24, t229=-t25, t230=-t27;
    double t231=-t28, t232=-t31, t233=-t32, t234=-t35, t235=-t37;
    double t236=-t38;
    
    /* More complex expressions */
    double t237 = s011*t44, t238 = s013*t42, t239 = s011*t45;
    double t240 = s022*t42, t241 = s031*t42, t242 = s101*t8;
    double t243 = s011*t8, t244 = t8*t47, t245 = t8*t42;
    double t246 = s013*t47, t247 = s103*t42, t248 = s110*t11;
    double t250 = t11*t50, t251 = s011*t52, t252 = t11*t42;
    double t253 = s013*t50, t254 = s022*t47, t255 = s112*t42;
    double t256 = s022*t50, t257 = s031*t47, t258 = s121*t42;
    double t259 = s031*t50, t260 = s130*t42, t261 = s101*t49;
    double t262 = s103*t47, t263 = s110*t20, t265 = t20*t50;
    double t266 = s101*t52, t267 = t20*t47, t268 = s103*t50;
    double t269 = s202*t42, t270 = s112*t47, t271 = s110*t52;
    double t272 = s211*t42, t273 = s121*t47, t274 = s112*t50;
    double t275 = s110*t53, t276 = s220*t42, t277 = s130*t47;
    double t278 = s121*t50, t279 = s130*t50, t280 = s101*t54;
    double t281 = s202*t47, t282 = s301*t42, t283 = s211*t47;
    double t284 = s202*t50, t285 = s110*t55, t286 = s310*t42;
    double t287 = s220*t47, t288 = s211*t50, t289 = s220*t50;
    double t290 = s301*t47, t291 = s310*t47, t292 = s301*t50;
    double t293 = s310*t50;
    
    double t294 = s110*t8*2.0;
    double t295 = t2*t50, t296 = t3*t50, t297 = t4*t47;
    double t298 = t4*t50, t299 = t5*t47, t300 = t6*t47;
    double t301 = t7*t50, t302 = t8*t50, t303 = t8*t20;
    double t304 = t8*t11, t306 = t11*t49, t307 = t20*t44;
    double t308 = t9*t50, t309 = t10*t50, t312 = t8*t53;
    double t314 = t20*t45, t315 = t14*t47, t316 = t12*t50;
    double t317 = t13*t50, t318 = t15*t47, t319 = t25*t42;
    double t320 = t16*t47, t322 = t17*t47, t324 = t18*t47;
    double t325 = t19*t50, t326 = t29*t42, t327 = t8*t55;
    double t328 = t11*t54, t330 = t23*t47, t331 = t21*t50;
    double t332 = t22*t50, t333 = t30*t42, t334 = t31*t42;
    double t335 = t26*t47, t336 = t24*t50, t337 = t32*t42;
    double t338 = t33*t42, t339 = t27*t47, t340 = t28*t47;
    double t341 = t34*t42, t342 = t29*t50, t343 = t35*t42;
    double t344 = t36*t42, t345 = t37*t42, t348 = t38*t42;
    double t349 = t34*t47, t350 = t39*t42, t351 = t40*t42;
    
    /* High-level neg */
    double t364 = -t43, t365 = -t44, t366 = -t45;
    double t367 = -t48, t368 = -t49, t369 = -t51;
    double t370 = -t52, t371 = -t53, t372 = -t54;
    double t373 = -t55, t378 = -t63, t381 = -t68;
    double t396 = -t91, t410 = -t109, t413 = -t114;
    double t415 = -t121, t423 = t8*t8, t424 = t11*t11;
    double t425 = t20*t20;
    
    double t429 = s013*s110*t8*-2.0;
    double t436 = s031*s110*t8*-2.0;
    double t437 = s103*s110*t8*-2.0;
    double t443 = s110*s112*t8*-2.0;
    double t449 = s110*s121*t8*-2.0;
    double t450 = s110*s130*t8*-2.0;
    double t462 = s110*s211*t8*-2.0;
    double t469 = s110*s301*t8*-2.0;
    double t470 = s110*s310*t8*-2.0;
    
    double t479 = t42*t135, t481 = t42*t136;
    double t486 = t47*t135, t487 = t42*t138;
    double t488 = t47*t136, t489 = t50*t135;
    double t492 = t50*t136, t495 = t47*t138;
    double t499 = t50*t138;
    
    double t353 = s022*t294, t356 = s112*t294, t357 = s121*t294;
    double t359 = s202*t294, t360 = s211*t294, t361 = s220*t294;
    
    double t374=-t57, t375=-t60, t376=-t61, t377=-t62, t379=-t64;
    double t380=-t67, t382=-t69, t383=-t70, t384=-t71, t385=-t72;
    double t386=-t294, t387=-t77, t388=-t78, t389=-t80, t390=-t81;
    double t391=-t82, t392=-t83, t393=-t88, t394=-t89, t395=-t90;
    double t397=-t93, t398=-t94, t399=-t95, t400=-t97, t401=-t98;
    double t402=-t99, t403=-t100, t404=-t101, t405=-t102, t406=-t103;
    double t407=-t104, t408=-t106, t409=-t107, t411=-t110, t412=-t111;
    double t414=-t116, t416=-t122, t417=-t126, t418=-t127, t419=-t128;
    double t420=-t129, t421=-t130, t422=-t131;
    
    double t426 = t302*2.0, t427 = t303*2.0, t428 = t304*2.0;
    
    double t430=-t146, t431=-t147, t432=-t148, t433=-t149, t434=-t150;
    double t435=-t151, t438=-t165, t439=-t166, t440=-t168, t441=-t169;
    double t442=-t170, t444=-t171, t445=-t172, t446=-t175, t447=-t177;
    double t448=-t178, t451=-t185, t452=-t186, t453=-t187, t454=-t188;
    double t455=-t189, t456=-t190, t457=-t193, t458=-t194, t459=-t195;
    double t460=-t196, t461=-t199, t463=-t200, t464=-t201, t465=-t202;
    double t466=-t203, t467=-t204, t468=-t205;
    
    double t471=-t237, t472=-t239, t473=-t240, t474=s101*t218;
    double t475=s011*t218, t476=s110*t220, t477=s011*t220;
    double t478=-t254, t480=-t256, t482=-t261, t483=s110*t225;
    double t484=s101*t225, t485=-t269, t490=-t275, t491=-t276;
    double t493=-t280, t494=-t281, t496=-t284, t497=-t285, t498=-t287;
    double t500=-t289, t501=-t295, t502=-t297, t503=-t298, t504=-t300;
    double t505=-t301, t506=t52*t218, t507=-t308, t508=t47*t221;
    double t509=t42*t226, t510=t52*t220, t511=-t315, t512=-t317;
    double t513=t50*t222, t514=-t322, t515=t42*t230, t516=-t324;
    double t517=-t326, t518=t52*t225;
    
    /* FINAL STAGE */
    double t519 = -t331, t520 = -t333, t521 = -t335;
    double t522 = -t338, t523 = -t341, t524 = -t342;
    double t525 = -t344, t526 = t47*t233, t527 = t50*t232;
    double t528 = -t349, t529 = -t350, t530 = -t351;
    double t531 = t8*t263*2.0, t532 = t8*t248*2.0;
    double t536 = s110*t423*-2.0, t534 = -t531, t535 = -t532;
    
    double t537 = t42+t47+t50+t386-1.0;
    double t542 = s110+t9+t14+t30+t135+t162+t163+t164+t167+t168+t169+t255+t270+t274+t369+t378+t379+t396+t397+t406+t409+t426+t443+t477+t484+t507+t511+t520;
    double t543 = s101+t13+t17+t33+t136+t173+t174+t175+t176+t178+t179+t258+t273+t278+t367+t380+t381+t401+t402+t410+t411+t427+t449+t475+t483+t512+t514+t522;
    double t544 = s011+t21+t26+t36+t138+t191+t192+t193+t196+t197+t198+t272+t283+t288+t364+t387+t388+t413+t414+t415+t416+t428+t462+t474+t476+t519+t521+t525;
    double t545 = s011+t2+t4+t22+t132+t141+t142+t143+t144+t145+t238+t246+t253+t307+t364+t374+t384+t385+t390+t392+t428+t429+t471+t474+t476+t501+t502+t509;
    double t546 = s011+t4+t6+t27+t133+t152+t153+t154+t155+t156+t241+t257+t259+t314+t364+t375+t389+t391+t393+t394+t428+t436+t472+t474+t476+t503+t504+t515;
    double t547 = s101+t7+t12+t29+t134+t157+t158+t159+t160+t161+t247+t262+t268+t306+t367+t376+t377+t395+t400+t403+t427+t437+t475+t482+t483+t505+t508+t517;
    double t548 = s110+t15+t18+t34+t137+t180+t181+t182+t183+t184+t260+t277+t279+t312+t369+t382+t383+t407+t408+t412+t426+t450+t477+t484+t490+t513+t516+t523;
    double t549 = s101+t29+t32+t39+t139+t206+t207+t208+t209+t210+t282+t290+t292+t328+t367+t398+t399+t417+t420+t421+t427+t469+t475+t483+t493+t524+t526+t529;
    double t550 = s110+t31+t34+t40+t140+t211+t212+t213+t214+t215+t286+t291+t293+t327+t369+t404+t405+t418+t419+t422+t426+t470+t477+t484+t497+t527+t528+t530;
    double t551 = s112+t64+t65+t92+t93+t105+t218+t219+t222+t232+t244+t245+t271+t302+t309+t318+t334+t356+t438+t439+t440+t441+t442+t479+t486+t489+t506+t536;
    
    double t538 = 1.0/t537;
    double t539 = s022+t58+t59+t74+t76+t86+t87+t216+t217+t229+t296+t299+t319+t353+t430+t431+t432+t433+t434+t435+t473+t478+t480+t537;
    double t540 = s202+t73+t75+t112+t113+t117+t120+t224+t227+t234+t325+t330+t343+t359+t451+t452+t453+t454+t455+t456+t485+t494+t496+t537;
    double t541 = s220+t84+t85+t118+t119+t124+t125+t228+t231+t236+t336+t340+t348+t361+t463+t464+t465+t466+t467+t468+t491+t498+t500+t537;
    double t552 = s121+t66+t67+t96+t108+t110+t220+t221+t223+t233+t250+t252+t266+t303+t316+t320+t337+t357+t444+t445+t446+t447+t448+t481+t488+t492+t510+t535;
    double t553 = s211+t79+t115+t116+t122+t123+t225+t226+t230+t235+t251+t265+t267+t304+t332+t339+t345+t360+t457+t458+t459+t460+t461+t487+t495+t499+t518+t534;
    
    double t554 = t538*t539, t555 = t538*t540, t556 = t538*t541;
    double t560 = t538*t542, t561 = t538*t543, t562 = t538*t544;
    double t563 = t538*t545, t564 = t538*t546, t565 = t538*t547;
    double t566 = t538*t548, t567 = t538*t549, t568 = t538*t550;
    double t569 = t538*t551, t570 = t538*t552, t571 = t538*t553;
    
    double t557 = -t554, t558 = -t555, t559 = -t556;
    double t572 = -t569, t573 = -t570, t574 = -t571;
    
    /* Compute matrix elements directly */
    double t50_sq = t50*t50;
    double t47_sq = t47*t47;
    double t42_sq = t42*t42;
    
    E[0] = -t538*(s400-t56+t372+t373+t537+s011*t37*2.0+s101*t39*2.0+s110*t40*2.0-s400*t42-s400*t47-s400*t50+s400*t294-t8*t40*2.0-t11*t39*2.0-t20*t37*2.0+t42*t56+t47*t55+t50*t54);
    E[1] = t568;
    E[2] = t567;
    E[3] = t559;
    E[4] = t562;
    E[5] = t558;
    E[6] = t568;
    E[7] = -t538*(s220-t50+t125*2.0-t202*2.0+t361+t370+t371+t373+t424+t425+t491+t498+t500+s011*t27*2.0+s101*t32*2.0-t11*t32*2.0-t20*t27*2.0-t8*t51*2.0+t42*t55+t47*t53+t50*t52+t50_sq);
    E[8] = t574;
    E[9] = t566;
    E[10] = t573;
    E[11] = t560;
    E[12] = t567;
    E[13] = t574;
    E[14] = -t538*(s202-t47+t113*2.0-t189*2.0+t359+t368+t370+t372+t423+t425+t485+t494+t496+s011*t22*2.0+s110*t31*2.0-t8*t31*2.0-t20*t22*2.0+t42*t54+t47*t52+t49*t50-t20*t242*2.0+t47_sq);
    E[15] = t561;
    E[16] = t572;
    E[17] = t565;
    E[18] = t559;
    E[19] = t566;
    E[20] = t561;
    E[21] = -t538*(s040-t46+t366+t371+t537+s011*t6*2.0-s040*t42-s040*t47-s040*t50+s101*t16*2.0+s110*t18*2.0+s040*t294-t6*t20*2.0-t8*t18*2.0-t11*t16*2.0+t46*t47+t42*t53+t45*t50);
    E[22] = t564;
    E[23] = t557;
    E[24] = t562;
    E[25] = t573;
    E[26] = t572;
    E[27] = t564;
    E[28] = -t538*(s022-t42+t59*2.0-t151*2.0+t353+t365+t366+t370+t423+t424+t473+t478+t480+s101*t12*2.0+s110*t15*2.0-t8*t15*2.0-t11*t12*2.0+t45*t47+t42*t52+t44*t50-t11*t243*2.0+t42_sq);
    E[29] = t563;
    E[30] = t558;
    E[31] = t560;
    E[32] = t565;
    E[33] = t557;
    E[34] = t563;
    E[35] = -t538*(s004-t41+t365+t368+t537+s011*t2*2.0-s004*t42-s004*t47-s004*t50+s101*t7*2.0+s110*t10*2.0+s004*t294-t7*t11*2.0-t8*t10*2.0-t2*t20*2.0+t41*t50+t42*t49+t44*t47);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double *s, *E;
    
    /* Check for proper number of arguments */
    if (nrhs != 28) {
        mexErrMsgIdAndTxt("delta2star3D:nrhs", "28 inputs required.");
    }
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("delta2star3D:nlhs", "One output required.");
    }
    
    /* Check that all inputs are scalar */
    for (int i = 0; i < 28; i++) {
        if (mxGetNumberOfElements(prhs[i]) != 1) {
            mexErrMsgIdAndTxt("delta2star3D:notScalar", "All inputs must be scalars.");
        }
    }
    
    /* Get input scalars */
    double s_array[28];
    for (int i = 0; i < 28; i++) {
        s_array[i] = mxGetScalar(prhs[i]);
    }
    
    /* Create output matrix */
    plhs[0] = mxCreateDoubleMatrix(6, 6, mxREAL);
    E = mxGetPr(plhs[0]);
    
    /* Call computation function */
    delta2star3D_compute(s_array, E);
}

