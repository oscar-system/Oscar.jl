```@meta
CurrentModule = Oscar
```

# Elliptic Surfaces
See [SS19](@cite) for the theory of elliptic surfaces.

```@docs
EllipticSurface
```

## Constructors
```@docs
elliptic_surface
kodaira_neron_model(E::EllipticCurve)
```

## Attributes
```@docs
generic_fiber(S::EllipticSurface)
zero_section(S::EllipticSurface)
euler_characteristic(X::EllipticSurface)
fibration(X::EllipticSurface)
weierstrass_chart(X::EllipticSurface)
weierstrass_chart_on_minimal_model(X::EllipticSurface)
weierstrass_model(X::EllipticSurface)
weierstrass_contraction(X::EllipticSurface)
fibration_on_weierstrass_model(X::EllipticSurface)
```

## Reducible fibers and the trivial lattice 
```@docs
trivial_lattice(X::EllipticSurface)
reducible_fibers(S::EllipticSurface)
fiber_cartier(S::EllipticSurface, P::Vector = ZZ.([0,1]))
fiber_components(S::EllipticSurface, P; algorithm=:exceptional_divisors)
```

## Mordell Weil group and related lattices
Since we require that ``\pi \colon X \to C`` has a section, the generic fibre of ``\pi`` is an elliptic curve ``E/k(C)``. The group ``E(k(C))`` of its rational points is called the Mordell-Weil group. It is a finitely generated abelian group. In general it is difficult to compute.

The height paring gives the free part of the Mordell-Weil group the structure of a quadratic lattice over the integers -- the so called Mordell-Weil lattice. 

```@docs
algebraic_lattice(X::EllipticSurface)
mordell_weil_sublattice(S::EllipticSurface)
mordell_weil_torsion(S::EllipticSurface)
section(X::EllipticSurface, P::EllipticCurvePoint)
basis_representation(X::EllipticSurface, D::AbsWeilDivisor)
```

### Updating the Mordell Weil group
Since the Mordell-Weil group is hard to compute, the user has to provide its generators, or at least generators of a subgroup
at construction of the surface. In some cases it may be convenient to enter the basis after creation, although it is in general not recommended. The following function are meant for this purpose.
```@docs
set_mordell_weil_basis!(X::EllipticSurface, mwl_basis::Vector{<:EllipticCurvePoint})
update_mwl_basis!(S::EllipticSurface, mwl_gens::Vector{<:EllipticCurvePoint})
algebraic_lattice_primitive_closure!(S::EllipticSurface)
algebraic_lattice_primitive_closure(S::EllipticSurface, p)
```

## Morphisms
```@docs
isomorphism_from_generic_fibers
isomorphisms(X::EllipticSurface, Y::EllipticSurface)
translation_morphism(X::EllipticSurface, P::EllipticCurvePoint)
```

## Fibration hopping
The methods in this section are available only for elliptically fibered K3 surfaces. 
A K3 surface ``X`` may admit several elliptic fibrations 
```math
\pi \colon X \to \mathbb{P}^1.
```
Fibration hopping is a way to compute them. 
See e.g. [BE23](@cite) and [BZ23](@cite).

A divisor ``F`` on ``X`` is called elliptic if it is primitive, isotropic and nef. 
The linear system ``|F|`` of an elliptic divisor ``F`` induces a genus one fibration on ``X``.
Conversely, the class of a fiber ``F`` of an elliptic fibration is an elliptic divisor.
Two elliptic fibrations on ``X`` with elliptic divisors ``F_1`` and ``F_2`` are called ``n``-neighbors if their intersection number is ``F_1.F_2=n``.
  
```@docs
linear_system(X::EllipticSurface, P::EllipticCurvePoint, k::Int)
two_neighbor_step(X::EllipticSurface, F::Vector{QQFieldElem})
elliptic_parameter(X::EllipticSurface, F::Vector{QQFieldElem})
pushforward_on_algebraic_lattices(f::MorphismFromRationalFunctions{<:EllipticSurface, <:EllipticSurface})
```
## Contact

Please direct questions about this part of OSCAR to the following people:
* [Simon Brandhorst](https://www.math.uni-sb.de/ag/brandhorst/index.php?lang=en).
* [Matthias Zach](https://math.rptu.de/en/wgs/agag/people/members),

You can ask questions in the [OSCAR Slack](https://www.oscar-system.org/community/#slack).

Alternatively, you can [raise an issue on github](https://www.oscar-system.org/community/#how-to-report-issues).
