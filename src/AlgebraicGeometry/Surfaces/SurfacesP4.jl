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
  I = load("$n", "Surfaces", :Oscardb)
  S = base_ring(I)
  R = grade(S)[1]
  #return proj(R, ideal(R, map(R, gens(I))))
  return variety(ideal(R, map(R, gens(I))), check = false, is_radical = true)
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

