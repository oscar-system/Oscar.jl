export abelian_d10_pi6
export abelian_d15_pi21_quintic_1
export abelian_d15_pi21_quintic_3
export bielliptic_d10_pi6
export bielliptic_d15_pi21
export bordiga
export castelnuovo
export cubic_scroll
export elliptic_d7_pi6
export elliptic_d8_pi7
export elliptic_d9_pi7
export elliptic_d10_pi9
export elliptic_d10_pi10
export elliptic_d11_pi12
export elliptic_d12_pi13
export elliptic_d12_pi14_ss_0
export elliptic_d12_pi14_ss_inf
export enriques_d9_pi6
export enriques_d10_pi8
export enriques_d11_pi10
export enriques_d13_pi16
export enriques_d13_pi16_two
export k3_d7_pi5
export k3_d8_pi6
export k3_d9_pi8
export k3_d10_pi9_quart_1
export k3_d10_pi9_quart_2
export k3_d11_pi11_ss_0
export k3_d11_pi11_ss_1
export k3_d11_pi11_ss_2
export k3_d11_pi11_ss_3
export k3_d11_pi12
export k3_d12_pi14
export k3_d13_pi16
export k3_d14_pi19
export quintic_elliptic_scroll
export rational_d7_pi4
export rational_d8_pi5
export rational_d8_pi6
export rational_d9_pi6
export rational_d9_pi7
export rational_d10_pi8
export rational_d10_pi9_quart_1
export rational_d10_pi9_quart_2
export rational_d11_pi11_ss_0
export rational_d11_pi11_ss_1
export rational_d11_pi11_ss_inf
export veronese

function surface(n::String)
  n = joinpath(oscardir, "data", "Surfaces", "$n" * ".mrdi")
  I = load(n)
  S = base_ring(I)
  R = grade(S)[1]
  return proj(R, ideal(R, map(R, gens(I))))
end

###############################

@doc raw"""
    abelian_d10_pi6()

Return a smooth abelian surface in $\mathbb P^4$ with degree `10` and sectional genus `6`.

The returned surface is defined over a prime field of characteristic 31991.
"""
abelian_d10_pi6() = surface("abelian_d10_pi6")


@doc raw"""
    abelian_d15_pi21_quintic_3()

Return a smooth abelian surface in $\mathbb P^4$ with degree `15` and sectional genus `21` which is contained in a net of quintics.

The returned surface is defined over a prime field of characteristic 31991.
"""
abelian_d15_pi21_quintic_3() = surface("abelian_d15_pi21_quintic_3")



@doc raw"""
    abelian_d15_pi21_quintic_1()

Return a smooth abelian surface in $\mathbb P^4$ with degree `15` and sectional genus `21` which is contained in precisely one quintic.

The returned surface is defined over a prime field of characteristic 31991.
"""
abelian_d15_pi21_quintic_1() = surface("abelian_d15_pi21_quintic_1")
   F = GF(31991)
   R, (x,y,z,u,v) = graded_polynomial_ring(F, ["x", "y", "z", "u", "v"])
   return proj(R, ideal([F(5)/8*x^3*y^3+2992*x^4*y*z-F(121)/98*x*y^2*z^3-F(96)/29*x^2*z^4+F(121)/89*x*y^3*z*u-4353*x^2*y*z^2*u+F(1)/2*z^5*u-F(8)/99*x^2*y^2*u^2-1983*x^3*z*u^2-117*y*z^3*u^2+F(6)/35*y^2*z*u^3-F(2)/5*x*z^2*u^3-1983*x*y*u^4+2992*x*y^4*v+F(96)/79*x^2*y^2*z*v+F(1)/2*x^3*z^2*v-F(5)/8*y*z^4*v+F(121)/89*x^3*y*u*v+F(19)/62*y^2*z^2*u*v-6369*x*z^3*u*v-1983*y^3*u^2*v-F(17)/42*x*y*z*u^2*v+F(6)/35*x^2*u^3*v-F(23)/11*z*u^4*v+F(1)/2*y^3*z*v^2+F(79)/60*x*y*z^2*v^2-4353*x*y^2*u*v^2+F(19)/62*x^2*z*u*v^2+14197*z^2*u^2*v^2-F(2)/5*y*u^3*v^2-F(121)/98*x^2*y*v^3+F(120)/29*z^3*v^3-6369*y*z*u*v^3-117*x*u^2*v^3-F(96)/29*y^2*v^4-F(5)/8*x*z*v^4+F(1)/2*u*v^5,F(1)/2*y^5*z-117*x*y^3*z^2+F(6)/35*x^2*y*z^3-F(5)/8*x*y^4*u+F(19)/62*x^2*y^2*z*u-1983*x^3*z^2*u-F(23)/11*y*z^4*u+F(1)/2*x^3*y*u^2+14197*y^2*z^2*u^2-F(2)/5*x*z^3*u^2+F(120)/29*y^3*u^3-6369*x*y*z*u^3-F(96)/29*x^2*u^4+F(1)/2*z*u^5-F(121)/98*x^2*y^3*v+F(121)/89*x^3*y*z*v-F(2)/5*y^2*z^3*v-1983*x*z^4*v+2992*x^4*u*v-6369*y^3*z*u*v-F(17)/42*x*y*z^2*u*v+F(79)/60*x*y^2*u^2*v-4353*x^2*z*u^2*v-117*z^2*u^3*v-F(5)/8*y*u^4*v-F(96)/29*y^4*v^2-4353*x*y^2*z*v^2-F(8)/99*x^2*z^2*v^2+F(96)/79*x^2*y*u*v^2+F(6)/35*z^3*u*v^2+F(19)/62*y*z*u^2*v^2-F(121)/98*x*u^3*v^2+F(5)/8*x^3*v^3-1983*y*z^2*v^3+F(1)/2*y^2*u*v^3+F(121)/89*x*z*u*v^3+2992*x*y*v^4,-x^2*y^3*z+1611*x^3*y*z^2+F(65)/99*y^2*z^4-F(21)/121*x*z^5+F(5)/4*x^3*y^2*u-F(19)/31*x^4*z*u-F(21)/121*y^3*z^2*u-10148*x*y*z^3*u-F(46)/11*x*y^2*z*u^2+3471*x^2*z^2*u^2-F(16)/99*x^2*y*u^3+F(97)/93*z^3*u^3+F(12)/35*y*z*u^4-F(21)/121*x*u^5-x^4*y*v-F(52)/99*y^4*z*v+F(122)/99*x*y^2*z^2*v-F(16)/99*x^2*z^3*v+F(13)/59*x*y^3*u*v+F(124)/101*x^2*y*z*u*v+F(12)/35*z^4*u*v+1611*x^3*u^2*v+9919*y*z^2*u^2*v+F(3)/7*y^2*u^3*v-10148*x*z*u^3*v+234*x^2*y^2*v^2+F(5)/4*x^3*z*v^2+F(3)/7*y*z^3*v^2+F(11)/4*y^2*z*u*v^2-F(46)/11*x*z^2*u*v^2+F(122)/99*x*y*u^2*v^2+F(65)/99*u^4*v^2-F(12)/35*y^3*v^3+F(13)/59*x*y*z*v^3-x^2*u*v^3-F(21)/121*z*u^2*v^3-F(52)/99*y*u*v^4,y^4*z^2-234*x*y^2*z^3+F(12)/35*x^2*z^4-F(90)/121*x*y^3*z*u+57*x^2*y*z^2*u-F(46)/11*z^5*u-F(99)/52*x^2*y^2*u^2-650*x^3*z*u^2-F(109)/4*y*z^3*u^2-2683*y^2*z*u^3-4373*x*z^2*u^3+F(16)/65*x*y*u^4+1611*u^6-F(72)/49*x^2*y^2*z*v+10613*x^3*z^2*v+6526*y*z^4*v-F(90)/121*x^3*y*u*v-3237*y^2*z^2*u*v-7261*x*z^3*u*v-650*y^3*u^2*v-F(4)/103*x*y*z*u^2*v-2683*x^2*u^3*v+F(11)/126*z*u^4*v+x^4*v^2+10613*y^3*z*v^2+3249*x*y*z^2*v^2+57*x*y^2*u*v^2-3237*x^2*z*u*v^2-F(70)/57*z^2*u^2*v^2-4373*y*u^3*v^2-234*x^2*y*v^3+9744*z^3*v^3-7261*y*z*u*v^3-F(109)/4*x*u^2*v^3+F(12)/35*y^2*v^4+6526*x*z*v^4-F(46)/11*u*v^5,-F(4)/5*y^6-F(113)/79*x*y^4*z-11797*x^2*y^2*z^2+F(96)/29*x^3*z^3+13879*x^2*y^3*u-15044*x^3*y*z*u-12643*y^2*z^3*u-F(1)/2*x*z^4*u+F(21)/121*x^4*u^2-F(121)/23*y^3*z*u^2+F(17)/38*x*y*z^2*u^2+3639*x*y^2*u^3+5300*x^2*z*u^3-F(5)/8*z^2*u^4+F(93)/109*y*u^5-12643*x^3*y^2*v-F(1)/2*x^4*z*v+13879*y^3*z^2*v-15044*x*y*z^3*v-F(69)/31*y^4*u*v-9826*x*y^2*z*u*v-F(22)/29*x^2*z^2*u*v-F(100)/23*x^2*y*u^2*v+9933*z^3*u^2*v-F(7)/80*y*z*u^3*v+3908*x*u^4*v-F(121)/23*x*y^3*v^2+F(17)/38*x^2*y*z*v^2+F(21)/121*z^4*v^2+9933*x^3*u*v^2-F(100)/23*y*z^2*u*v^2+F(49)/53*y^2*u^2*v^2+2111*x*z*u^2*v^2+3639*y^2*z*v^3+5300*x*z^2*v^3-F(7)/80*x*y*u*v^3-12076*u^3*v^3-F(5)/8*x^2*v^4+3908*z*u*v^4+F(93)/109*y*v^5,-F(5)/8*x^2*y^4+9933*x^3*y^2*z+F(21)/121*x^4*z^2-12076*y^3*z^3+3908*x*y*z^4-F(1)/2*x^4*y*u+3908*y^4*z*u+2111*x*y^2*z^2*u+5300*x^2*z^3*u+5300*x*y^3*u^2-F(22)/29*x^2*y*z*u^2-F(5)/8*z^4*u^2+F(96)/29*x^3*u^3+9933*y*z^2*u^3+F(21)/121*y^2*u^4-F(1)/2*x*z*u^4+F(93)/109*y^5*v-F(7)/80*x*y^3*z*v-F(100)/23*x^2*y*z^2*v+F(93)/109*z^5*v+F(17)/38*x^2*y^2*u*v-15044*x^3*z*u*v-F(7)/80*y*z^3*u*v-F(100)/23*y^2*z*u^2*v+F(17)/38*x*z^2*u^2*v-15044*x*y*u^3*v-12643*x^3*y*v^2+F(49)/53*y^2*z^2*v^2+3639*x*z^3*v^2+3639*y^3*u*v^2-9826*x*y*z*u*v^2-11797*x^2*u^2*v^2-12643*z*u^3*v^2-F(121)/23*x*y^2*v^3+13879*x^2*z*v^3-F(121)/23*z^2*u*v^3+13879*y*u^2*v^3-F(69)/31*y*z*v^4-F(113)/79*x*u*v^4-F(4)/5*v^6,F(33)/122*x^3*y^2*z+4022*x^4*z^2+F(52)/99*y^3*z^3-F(122)/117*x*y*z^4-F(52)/99*x^4*y*u-F(122)/117*y^4*z*u+F(83)/26*x*y^2*z^2*u+F(89)/58*x^2*z^3*u+F(89)/58*x*y^3*u^2+F(99)/70*x^2*y*z*u^2-F(12)/35*x^3*u^3+F(33)/122*y*z^2*u^3+4022*y^2*u^4-F(52)/99*x*z*u^4-F(89)/58*y^5*v-F(100)/71*x*y^3*z*v+F(78)/101*x^2*y*z^2*v-F(89)/58*z^5*v+F(63)/94*x^2*y^2*u*v-F(17)/31*x^3*z*u*v-F(100)/71*y*z^3*u*v+F(78)/101*y^2*z*u^2*v+F(63)/94*x*z^2*u^2*v-F(17)/31*x*y*u^3*v-2602*x^3*y*v^2+11793*y^2*z^2*v^2+8051*x*z^3*v^2+8051*y^3*u*v^2+12945*x*y*z*u*v^2+F(108)/109*x^2*u^2*v^2-2602*z*u^3*v^2-2520*x*y^2*v^3-F(95)/52*x^2*z*v^3-2520*z^2*u*v^3-F(95)/52*y*u^2*v^3+8051*y*z*v^4-F(65)/83*x*u*v^4-F(89)/58*v^6,-F(89)/58*y^6-F(65)/83*x*y^4*z+F(108)/109*x^2*y^2*z^2-F(12)/35*x^3*z^3-F(95)/52*x^2*y^3*u-F(17)/31*x^3*y*z*u-2602*y^2*z^3*u-F(52)/99*x*z^4*u+4022*x^4*u^2-2520*y^3*z*u^2+F(63)/94*x*y*z^2*u^2+8051*x*y^2*u^3+F(89)/58*x^2*z*u^3-F(89)/58*y*u^5-2602*x^3*y^2*v-F(52)/99*x^4*z*v-F(95)/52*y^3*z^2*v-F(17)/31*x*y*z^3*v+8051*y^4*u*v+12945*x*y^2*z*u*v+F(99)/70*x^2*z^2*u*v+F(78)/101*x^2*y*u^2*v+F(33)/122*z^3*u^2*v-F(100)/71*y*z*u^3*v-F(122)/117*x*u^4*v-2520*x*y^3*v^2+F(63)/94*x^2*y*z*v^2+4022*z^4*v^2+F(33)/122*x^3*u*v^2+F(78)/101*y*z^2*u*v^2+11793*y^2*u^2*v^2+F(83)/26*x*z*u^2*v^2+8051*y^2*z*v^3+F(89)/58*x*z^2*v^3-F(100)/71*x*y*u*v^3+F(52)/99*u^3*v^3-F(122)/117*z*u*v^4-F(89)/58*y*v^5,x^4*y^2+1611*y^5*z-F(77)/60*x*y^3*z^2+F(76)/69*x^2*y*z^3+1611*z^6-F(4)/5*x*y^4*u-15475*x^2*y^2*z*u+F(86)/19*x^3*z^2*u-F(77)/60*y*z^4*u-11038*x^3*y*u^2+F(58)/25*y^2*z^2*u^2+F(112)/43*x*z^3*u^2-F(21)/121*y^3*u^3-11234*x*y*z*u^3+1611*z*u^5-1611*x^2*y^3*v+F(55)/73*x^3*y*z*v+F(112)/43*y^2*z^3*v-11235*x*z^4*v+F(21)/121*x^4*u*v-11234*y^3*z*u*v-F(43)/2*x*y*z^2*u*v-12177*x*y^2*u^2*v+F(18)/107*x^2*z*u^2*v-F(77)/60*z^2*u^3*v-F(4)/5*y*u^4*v+F(18)/107*x*y^2*z*v^2+10973*x^2*z^2*v^2+2311*x^2*y*u*v^2+F(76)/69*z^3*u*v^2-15475*y*z*u^2*v^2-1611*x*u^3*v^2-F(65)/99*x^3*v^3+F(86)/19*y*z^2*v^3-11038*y^2*u*v^3+F(55)/73*x*z*u*v^3+F(21)/121*x*y*v^4+u^2*v^4,-x*y^5+234*x^2*y^3*z-F(12)/35*x^3*y*z^2+F(4)/5*x^3*y^2*u+F(21)/121*x^4*z*u+F(121)/49*y^3*z^2*u+F(25)/123*x*y*z^3*u+11038*y^4*u^2+8706*x*y^2*z*u^2+F(16)/99*x^2*z^2*u^2+F(21)/121*x^2*y*u^3-F(5)/4*z^3*u^3-5984*y*z*u^4+F(46)/11*x^4*y*v+F(5)/4*y^4*z*v-F(19)/31*x*y^2*z^2*v+F(21)/121*x^2*z^3*v+F(91)/108*x*y^3*u*v+F(17)/21*x^2*y*z*u*v-5984*z^4*u*v-F(12)/35*x^3*u^2*v+12551*y*z^2*u^2*v-y^2*u^3*v+F(25)/123*x*z*u^3*v+3597*x^2*y^2*v^2+F(4)/5*x^3*z*v^2-y*z^3*v^2-F(79)/30*y^2*z*u*v^2+8706*x*z^2*u*v^2-F(19)/31*x*y*u^2*v^2+2198*y^3*v^3+F(91)/108*x*y*z*v^3+234*x^2*u*v^3+F(121)/49*z*u^2*v^3+11038*z^2*v^4+F(5)/4*y*u*v^4-x*v^5,-x^5-y^5+7884*x*y^3*z-F(99)/23*x^2*y*z^2-z^5-3570*x^2*y^2*u+F(3)/100*x^3*z*u+7884*y*z^3*u-F(99)/23*y^2*z*u^2-3570*x*z^2*u^2+F(3)/100*x*y*u^3-u^5+7884*x^3*y*v-3570*y^2*z^2*v+F(3)/100*x*z^3*v+F(3)/100*y^3*u*v+F(125)/81*x*y*z*u*v-F(99)/23*x^2*u^2*v+7884*z*u^3*v-F(99)/23*x*y^2*v^2-3570*x^2*z*v^2-F(99)/23*z^2*u*v^2-3570*y*u^2*v^2+F(3)/100*y*z*v^3+7884*x*u*v^3-v^5]))



@doc raw"""
    bielliptic_d10_pi6()

Return a smooth bielliptic surface in $\mathbb P^4$ with degree `10` and sectional genus `6`.

The returned surface is defined over a prime field of characteristic 911.
"""
bielliptic_d10_pi6() = surface("bielliptic_d10_pi6")



@doc raw"""
    bielliptic_d15_pi21()

Return a smooth bielliptic surface in $\mathbb P^4$ with degree `15` and sectional genus `21`.

The returned surface is defined over a prime field of characteristic 911.
"""
bielliptic_d15_pi21() = surface("bielliptic_d15_pi21")



@doc raw"""
    bordiga()

Return a smooth rational surface in $\mathbb P^4$ with degree `6` and sectional genus `3`.

The returned surface is defined over a prime field of characteristic 31991.
"""
bordiga() = surface("bordiga")



@doc raw"""
    castelnuovo()

Return a smooth rational surface in $\mathbb P^4$ with degree `5` and sectional genus `2`.

The returned surface is defined over a prime field of characteristic 31991.
"""
castelnuovo() = surface("castelnuovo")



@doc raw"""
    cubic_scroll()

Return a smooth rational surface in $\mathbb P^4$ with degree `3` and sectional genus `0`.

The returned surface is defined over a prime field of characteristic 31991.
"""
cubic_scroll() = surface("cubic_scroll")



@doc raw"""
    elliptic_d10_pi10()

Return a smooth elliptic surface in $\mathbb P^4$ with degree `10` and sectional genus `10`.

The returned surface is defined over a prime field of characteristic 31991.
"""
elliptic_d10_pi10() = surface("elliptic_d10_pi10")



@doc raw"""
    elliptic_d10_pi9()

Return a smooth elliptic surface in $\mathbb P^4$ with degree `10` and sectional genus `9`.

The returned surface is defined over a prime field of characteristic 31991.
"""
elliptic_d10_pi9() = surface("elliptic_d10_pi9")



@doc raw"""
    elliptic_d11_pi12()

Return a smooth elliptic surface in $\mathbb P^4$ with degree `11` and sectional genus `12`.

The returned surface is defined over a prime field of characteristic 31991.
"""
elliptic_d11_pi12() = surface("elliptic_d11_pi12")



@doc raw"""
    elliptic_d12_pi13()

Return a smooth elliptic surface in $\mathbb P^4$ with degree `12` and sectional genus `13`.

The returned surface is defined over a prime field of characteristic 31991.
"""
elliptic_d12_pi13() = surface("elliptic_d12_pi13")



@doc raw"""
    elliptic_d12_pi14_ss_0()

Return a smooth elliptic surface in $\mathbb P^4$ with degree `12`, sectional genus `14`, and no 6-secant.

The returned surface is defined over a prime field of characteristic 31991.
"""
elliptic_d12_pi14_ss_0() = surface("elliptic_d12_pi14_ss_0")



@doc raw"""
    elliptic_d12_pi14_ss_inf()

Return a smooth elliptic surface in $\mathbb P^4$ with degree `12`, sectional genus `14`, and infinitely many 6-secants.

The returned surface is defined over a prime field of characteristic 31991.
"""
elliptic_d12_pi14_ss_inf() = surface("elliptic_d12_pi14_ss_inf")



@doc raw"""
    elliptic_d7_pi6()

Return a smooth elliptic surface in $\mathbb P^4$ with degree `7` and sectional genus `6`.

The returned surface is defined over a prime field of characteristic 31991.
"""
elliptic_d7_pi6() = surface("elliptic_d7_pi6")



@doc raw"""
    elliptic_d8_pi7()

Return a smooth elliptic surface in $\mathbb P^4$ with degree `8` and sectional genus `7`.

The returned surface is defined over a prime field of characteristic 31991.
"""
elliptic_d8_pi7() = surface("elliptic_d8_pi7")



@doc raw"""
    elliptic_d9_pi7()

Return a smooth elliptic surface in $\mathbb P^4$ with degree `9` and sectional genus `7`.

The returned surface is defined over a prime field of characteristic 31991.
"""
elliptic_d9_pi7() = surface("elliptic_d9_pi7")



@doc raw"""
    quintic_elliptic_scroll()

Return a smooth ruled surface in $\mathbb P^4$ with degree `5` and sectional genus `1`.

The returned surface is defined over a prime field of characteristic 31991.
"""
 quintic_elliptic_scroll() = surface(" quintic_elliptic_scroll")



@doc raw"""
    enriques_d10_pi8()

Return a smooth Enriques surface in $\mathbb P^4$ with degree `10` and sectional genus `8`.

The returned surface is defined over a prime field of characteristic 31991.
"""
enriques_d10_pi8() = surface("enriques_d10_pi8")



@doc raw"""
    enriques_d11_pi10()

Return a smooth Enriques surface in $\mathbb P^4$ with degree `11` and sectional genus `10`.

The returned surface is defined over a prime field of characteristic 43.
"""
enriques_d11_pi10() = surface("enriques_d11_pi10")



@doc raw"""
    enriques_d13_pi16()

Return a smooth Enriques surface in $\mathbb P^4$ with degree `13` and sectional genus `16`.

The returned surface is defined over a prime field of characteristic 31991.
"""
enriques_d13_pi16() = surface("enriques_d13_pi16")



@doc raw"""
    enriques_d13_pi16_two()

Return a smooth Enriques surface in $\mathbb P^4$ with degree `13` and sectional genus `16`.

The returned surface is defined over a prime field of characteristic 31991.
"""
enriques_d13_pi16_two() = surface("enriques_d13_pi16_two")



@doc raw"""
    enriques_d9_pi6()

Return a smooth Enriques surface in $\mathbb P^4$ with degree `9` and sectional genus `6`.

The returned surface is defined over a prime field of characteristic 31991.
"""
enriques_d9_pi6() = surface("enriques_d9_pi6")



@doc raw"""
    k3_d10_pi9_quart_1()

Return a smooth K3 surface in $\mathbb P^4$ with degree `10` and sectional genus `9` which is contained in precisely one quartic.

The returned surface is defined over a prime field of characteristic 31991.
"""
k3_d10_pi9_quart_1() = surface("k3_d10_pi9_quart_1")



@doc raw"""
    k3_d10_pi9_quart_2()

Return a smooth K3 surface in $\mathbb P^4$ with degree `10` and sectional genus `9` which is contained in a pencil of quartics.

The returned surface is defined over a prime field of characteristic 31991.
"""
k3_d10_pi9_quart_2() = surface("k3_d10_pi9_quart_2")



@doc raw"""
    k3_d11_pi11_ss_0()

Return a smooth K3 surface in $\mathbb P^4$ with degree `11`, sectional genus `11`, and no 6-secant.

The returned surface is defined over a prime field of characteristic 31991.
"""
k3_d11_pi11_ss_0() = surface("k3_d11_pi11_ss_0")



@doc raw"""
    k3_d11_pi11_ss_1()

Return a smooth K3 surface in $\mathbb P^4$ with degree `11`, sectional genus `11`, and one 6-secant.

The returned surface is defined over a prime field of characteristic 31991.
"""
k3_d11_pi11_ss_1() = surface("k3_d11_pi11_ss_1")



@doc raw"""
    k3_d11_pi11_ss_2()

Return a smooth K3 surface in $\mathbb P^4$ with degree `11`, sectional genus `11`, and two 6-secants.

The returned surface is defined over a prime field of characteristic 31991.
"""
k3_d11_pi11_ss_2() = surface("k3_d11_pi11_ss_2")



@doc raw"""
    k3_d11_pi11_ss_3()

Return a smooth K3 surface in $\mathbb P^4$ with degree `11`, sectional genus `11`, and three 6-secants.

The returned surface is defined over a prime field of characteristic 31991.
"""
k3_d11_pi11_ss_3() = surface("k3_d11_pi11_ss_3")



@doc raw"""
    k3_d11_pi12()

Return a smooth K3 surface in $\mathbb P^4$ with degree `11` and sectional genus `12`.

The returned surface is defined over a prime field of characteristic 31991.
"""
k3_d11_pi12() = surface("k3_d11_pi12")



@doc raw"""
    k3_d12_pi14()

Return a smooth K3 surface in $\mathbb P^4$ with degree `12` and sectional genus `14`.

The returned surface is defined over a prime field of characteristic 31991.
"""
k3_d12_pi14() = surface("k3_d12_pi14")



@doc raw"""
    k3_d13_pi16()

Return a smooth K3 surface in $\mathbb P^4$ with degree `13` and sectional genus `16`.

The returned surface is defined over a prime field of characteristic 31991.
"""
k3_d13_pi16() = surface("k3_d13_pi16")



@doc raw"""
    k3_d14_pi19()

Return a smooth K3 surface in $\mathbb P^4$ with degree `14` and sectional genus `19`.

The returned surface is defined over a prime field of characteristic 31991.
"""
k3_d14_pi19() = surface("k3_d14_pi19")



@doc raw"""
    k3_d7_pi5()

Return a smooth K3 surface in $\mathbb P^4$ with degree `7` and sectional genus `5`.

The returned surface is defined over a prime field of characteristic 31991.
"""
k3_d7_pi5() = surface("k3_d7_pi5")



@doc raw"""
    k3_d8_pi6()

Return a smooth K3 surface in $\mathbb P^4$ with degree `8` and sectional genus `6`.

The returned surface is defined over a prime field of characteristic 31991.
"""
k3_d8_pi6() = surface("k3_d8_pi6")



@doc raw"""
    k3_d9_pi8()

Return a smooth K3 surface in $\mathbb P^4$ with degree `9` and sectional genus `8`.

The returned surface is defined over a prime field of characteristic 31991.
"""
k3_d9_pi8() = surface("k3_d9_pi8")



@doc raw"""
    rational_d10_pi8()

Return a smooth rational surface in $\mathbb P^4$ with degree `10` and sectional genus `8`.

The returned surface is defined over a prime field of characteristic 31991.
"""
rational_d10_pi8() = surface("rational_d10_pi8")



@doc raw"""
    rational_d10_pi9_quart_1()

Return a smooth rational surface in $\mathbb P^4$ with degree `10` and sectional genus `9` which is contained in precisely one quartic.

The returned surface is defined over a prime field of characteristic 31991.
"""
rational_d10_pi9_quart_1() = surface("rational_d10_pi9_quart_1")



@doc raw"""
    rational_d10_pi9_quart_2()

Return a smooth rational surface in $\mathbb P^4$ with degree `10` and sectional genus `9` which is contained in a pencil of quartics.

The returned surface is defined over a prime field of characteristic 31991.
"""
rational_d10_pi9_quart_2() = surface("rational_d10_pi9_quart_2")



@doc raw"""
    rational_d11_pi11_ss_0()

Return a smooth rational surface in $\mathbb P^4$ with degree `11`, sectional genus `11`, and no 6-secant.

The returned surface is defined over a prime field of characteristic 31991.
"""
rational_d11_pi11_ss_0() = surface("rational_d11_pi11_ss_0")



@doc raw"""
    rational_d11_pi11_ss_1()

Return a smooth rational surface in $\mathbb P^4$ with degree `11`, sectional genus `11`, and one 6-secant.

The returned surface is defined over a prime field of characteristic 31991.
"""
rational_d11_pi11_ss_1() = surface("rational_d11_pi11_ss_1")



@doc raw"""
    rational_d11_pi11_ss_inf()

Return a smooth rational surface in $\mathbb P^4$ with degree `11`, sectional genus `11`, and infinitely many 6-secants.

The returned surface is defined over a prime field of characteristic 31991.
"""
rational_d11_pi11_ss_inf() = surface("rational_d11_pi11_ss_inf")



@doc raw"""
    rational_d7_pi4()

Return a smooth rational surface in $\mathbb P^4$ with degree `7` and sectional genus `4`.

The returned surface is defined over a prime field of characteristic 31991.
"""
rational_d7_pi4() = surface("rational_d7_pi4")



@doc raw"""
    rational_d8_pi5()

Return a smooth rational surface in $\mathbb P^4$ with degree `8` and sectional genus `5`.

The returned surface is defined over a prime field of characteristic 31991.
"""
rational_d8_pi5() = surface("rational_d8_pi5")



@doc raw"""
    rational_d8_pi6()

Return a smooth rational surface in $\mathbb P^4$ with degree `8` and sectional genus `6`.

The returned surface is defined over a prime field of characteristic 31991.
"""
rational_d8_pi6() = surface("rational_d8_pi6")



@doc raw"""
    rational_d9_pi6()

Return a smooth rational surface in $\mathbb P^4$ with degree `9` and sectional genus `6`.

The returned surface is defined over a prime field of characteristic 31991.
"""
rational_d9_pi6() = surface("rational_d9_pi6")



@doc raw"""
    rational_d9_pi7()

Return a smooth rational surface in $\mathbb P^4$ with degree `9` and sectional genus `7`.

The returned surface is defined over a prime field of characteristic 31991.
"""
rational_d9_pi7() = surface("rational_d9_pi7")



@doc raw"""
    veronese()

Return a smooth rational surface in $\mathbb P^4$ with degree `4` and sectional genus `0`.

The returned surface is defined over a prime field of characteristic 31991.
"""
veronese() = surface("veronese")

