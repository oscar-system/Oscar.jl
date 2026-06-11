```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Introduction

The central objetct here are $G$-modules: 
Let $G$ be a finite group and $M$ be any abelian group. A $G$-module
$C$ is given by specifying automorphisms of $M$ for each generator of $G$, 
thus as a $G$-action on $M$. Here $M$ can have a wide range of types, 
including
 - abelian groups
 - fp-modules (over fields/ the integers)
In fact, as long as the homomorphisms are sufficiently powerful, the $G$-module
code may work.


```jldoctest
julia> k, a = cyclotomic_field(11);

julia> u, mu = unit_group(maximal_order(k));

julia> C = gmodule(automorphism_group(PermGroup, k)[1], mu)
G-module for permutation group of degree 10 acting on u

julia> G = group(C)
Permutation group of degree 10

julia> C = natural_gmodule(G, QQ)
G-module for G acting on vector space of dimension 10 over QQ

```

The aim here is to model and work with
 - compute (Galois-) cohomology
 - group extensions
 - representations

# Creation
```@docs
gmodule
action
induce(C::GModule{GT, MT}, h::Map, D = nothing, mDC = nothing) where GT <: GAPGroup where MT
cohomology_group
extension
extension_with_abelian_kernel
is_central(NtoE::Map{<:Union{Group, FinGenAbGroup}, <:Group})
restriction_of_scalars
trivial_gmodule
regular_gmodule
natural_gmodule
gmodule_minimal_field
gmodule_over
factor_set
composition_factors_with_multiplicity(C::GModule{<:Any, <:AbstractAlgebra.FPModule{<:FinFieldElem}})
indecomposition
ghom
units_mod_ideal
is_coboundary
decomposition_group(K::AbsSimpleNumField, mK::Map, mG::Map = automorphism_group(K)[2], mGp::Map = automorphism_group(codomain(mK), absolute_base_field(codomain(mK))); _sub::Bool = false)
decomposition_group(K::AbsSimpleNumField, emb::Hecke.NumFieldEmb, mG::Map = automorphism_group(PermGroup, K)[2])
idele_class_gmodule
galois_group(A::ClassField, ::QQField)
relative_brauer_group
brauer_group
structure_constant_algebra(a::Oscar.GaloisCohomology_Mod.RelativeBrauerGroupElem)
structure_constant_algebra(CC::GrpCoh.CoChain{2, PermGroupElem, Oscar.GrpCoh.MultGrpElem{AbsSimpleNumFieldElem}}, mG::Map = automorphism_group(PermGroup, CC.C.M.data)[2], mkK::Union{<:Map, Nothing} = nothing)
```

## Tutorials

We encourage you to take a look at the tutorials on intersection theory in OSCAR,
which can be found [here](https://www.oscar-system.org/tutorials/IntersectionTheory/).


## Contact

Please direct questions about this part of OSCAR to the following
people:
* [Claus Fieker](https://math.rptu.de/en/wgs/agag/people/head/fieker),

You can ask questions in the [OSCAR Slack](https://www.oscar-system.org/community/#slack).

Alternatively, you can [raise an issue on github](https://www.oscar-system.org/community/#how-to-report-issues).
