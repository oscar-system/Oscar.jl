# Nongeneral Type Surfaces in $\mathbb P^4$

Every smooth, projective surface can be embedded in $\mathbb P^5$, but there are constraints
on the numerical invariants of a smooth surface in $\mathbb P^4$: The invariants of each such
surface $S$ satisfy the double point formula

$d^2-5d-10(\pi-1)+2(\chi(\mathcal O_S)-K_S^2) = 0.$

Here, $d$ is the degree of $S$, $\pi$ its sectional genus, $\chi(\mathcal O_S)$ its Euler-Poincare
characteristic, and $K_S$ its canonical class.
The double point  formula is a key ingredient in the proof of a
theorem of Ellingsrud and Peskine which states that there are only finitely many families of
smooth surfaces in $\mathbb P^4$ which are not of general type. That is, the degree of such
surfaces in bounded from above. The best bound known so far is $52$, while
examples exist up to degree $15$.

For details, we refer to
- [DES93](@cite)
- [Pop93](@cite)
- [DS00](@cite)
and the references given there.

Below, we present functions which return one hard coded example for each family presented
in the first two papers above.  Based on these papers, the existence of further families has
been shown. Hard coded OSCAR examples for these surfaces are under construction.

!!! note
    To ease subsequent computations, all hard coded examples are defined over a finite prime field.


## Rational Surfaces

#### A Rational Surface with $d=3$, $\pi=0$

```@docs
cubic_scroll()
```

#### A Rational Surface with $d=4$, $\pi=0$

```@docs
veronese()
```

#### A Rational Surface with $d=5$, $\pi=2$

```@docs
castelnuovo()
```

#### A Rational Surface with $d=6$, $\pi=3$

```@docs
bordiga()
```

#### A Rational Surface with $d=7$, $\pi=4$

```@docs
rational_d7_pi4()
```

#### A Rational Surface with $d=8$, $\pi=5$

```@docs
rational_d8_pi5()
```

#### A Rational Surface with $d=8$, $\pi=6$

```@docs
rational_d8_pi6()
```

#### A Rational Surface with $d=9$, $\pi=6$

```@docs
rational_d9_pi6()
```

#### A Rational Surface with $d=9$, $\pi=7$

```@docs
rational_d9_pi7()
```

#### A Rational Surface with $d=10$, $\pi=8$

```@docs
rational_d10_pi8()
```

#### A Rational Surface with $d=10$, $\pi=9$ which is Contained in one Quartic

```@docs
rational_d10_pi9_quart_1()
```

#### A Rational Surface with $d=10$, $\pi=9$ which is Contained in a Pencil of Quartics

```@docs
rational_d10_pi9_quart_2()
```

#### A Rational Surface with $d=11$, $\pi=11$, and no 6-Secant

```@docs
rational_d11_pi11_ss_0()
```

#### A Rational Surface with $d=11$, $\pi=11$, and one 6-Secant

```@docs
rational_d11_pi11_ss_1()
```

#### A Rational Surface with $d=11$, $\pi=11$, and Infinitely Many 6-Secants

```@docs
rational_d11_pi11_ss_inf()
```

## Ruled Surfaces

#### A Ruled Surface with $d=5$, $\pi=1$

```@docs
quintic_elliptic_scroll()
```

## Enriques Surfaces

#### An Enriques Surface with $d=9$, $\pi=6$

```@docs
enriques_d9_pi6()
```

#### An Enriques Surface with $d=10$, $\pi=8$

```@docs
enriques_d10_pi8()
```

#### An Enriques Surface with $d=11$, $\pi=10$

```@docs
enriques_d11_pi10()
```

#### An Enriques Surface with $d=13$, $\pi=16$

```@docs
enriques_d13_pi16()
```

#### An Enriques Surface with $d=13$, $\pi=16$

```@docs
enriques_d13_pi16_two()
```

## K3 Surfaces

#### A K3 Surface with $d=7$, $\pi=5$

```@docs
k3_d7_pi5
```

#### A K3 Surface with $d=8$, $\pi=6$

```@docs
k3_d8_pi6
```

#### A K3 Surface with $d=9$, $\pi=8$

```@docs
k3_d9_pi8
```

#### A K3 Surface with $d=10$, $\pi=9$ which is Contained in one Quartic

```@docs
k3_d10_pi9_quart_1()
```

#### A K3 Surface with $d=10$, $\pi=9$ which is Contained in a Pencil of Quartics

```@docs
k3_d10_pi9_quart_2()
```

#### A K3 Surface with $d=11$, $\pi=11$ and no 6-Secant

```@docs
k3_d11_pi11_ss_0()
```

#### A K3 Surface with $d=11$, $\pi=11$ and one 6-Secant

```@docs
k3_d11_pi11_ss_1()
```

#### A K3 Surface with $d=11$, $\pi=11$ and two 6-Secants

```@docs
k3_d11_pi11_ss_2()
```

#### A K3 Surface with $d=11$, $\pi=11$ and three 6-Secants

```@docs
k3_d11_pi11_ss_3()
```

#### A K3 Surface with $d=11$, $\pi=12$

```@docs
k3_d11_pi12()
```


#### A K3 Surface with $d=12$, $\pi=14$

```@docs
k3_d12_pi14()
```

#### A K3 Surface with $d=13$, $\pi=16$

```@docs
k3_d13_pi16()
```

#### A K3 Surface with $d=14$, $\pi=19$

```@docs
k3_d14_pi19()
```

## Bielliptic Surfaces

#### A Bielliptic Surface with $d=10$, $\pi=6$

```@docs
bielliptic_d10_pi6()
```

#### A Bielliptic Surface with $d=15$, $\pi=21$

```@docs
bielliptic_d15_pi21()
```

## Abelian Surfaces

#### An Abelian Surface with $d=10$, $\pi=6$

```@docs
abelian_d10_pi6()
```

#### An Abelian Surface with $d=15$, $\pi=21$ which is Contained in a Net of Quintics

```@docs
abelian_d15_pi21_quintic_3()
```

#### An Abelian Surface with $d=15$, $\pi=21$ which is Contained in one Quintic

```@docs
abelian_d15_pi21_quintic_1()
```

## Elliptic Surfaces

#### An Elliptic Surface with $d=7$, $\pi=6$

```@docs
elliptic_d7_pi6()
```

#### An Elliptic Surface with $d=8$, $\pi=7$

```@docs
elliptic_d8_pi7()
```

#### An Elliptic Surface with $d=9$, $\pi=7$

```@docs
elliptic_d9_pi7()
```

#### An Elliptic Surface with $d=10$, $\pi=9$

```@docs
elliptic_d10_pi9()
```

#### An Elliptic Surface with $d=10$, $\pi=10$

```@docs
elliptic_d10_pi10()
```

#### An Elliptic Surface with $d=11$, $\pi=12$

```@docs
elliptic_d11_pi12()
```

#### An Elliptic Surface with $d=12$, $\pi=13$

```@docs
elliptic_d12_pi13()
```

#### An Elliptic Surface with $d=12$, $\pi=14$ and no 6-Secant

```@docs
elliptic_d12_pi14_ss_0()
```

#### An Elliptic Surface with $d=12$, $\pi=14$, and Infinitely Many 6-Secants

```@docs
elliptic_d12_pi14_ss_inf()
```

## Contact

Please direct questions about this part of OSCAR to the following people:
* [Wolfram Decker](https://math.rptu.de/en/wgs/agag/people/head/decker).

You can ask questions in the [OSCAR Slack](https://www.oscar-system.org/community/#slack).

Alternatively, you can [raise an issue on github](https://www.oscar-system.org/community/#how-to-report-issues).

