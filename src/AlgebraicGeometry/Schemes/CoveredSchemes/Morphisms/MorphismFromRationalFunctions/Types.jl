########################################################################
# Morphisms from rational functions                                    #
########################################################################
@doc raw"""
    MorphismFromRationalFunctions{DomainType<:AbsCoveredScheme, CodomainType<:AbsCoveredScheme} 

A lazy type for a dominant morphism ``φ : X → Y`` of `AbsCoveredScheme`s which is given 
by a set of rational functions ``a₁,…,aₙ`` in the fraction field of the `base_ring`
of ``𝒪(U)`` for one of the dense open `affine_chart`s ``U`` of ``X``. 
The ``aᵢ`` represent the pullbacks of the coordinates (`gens`) of some 
`affine_chart` ``V`` of the codomain ``Y`` under this map. 
```jldoctest
julia> IP1 = covered_scheme(projective_space(QQ, [:s, :t]))
Scheme
  over rational field
with default covering
  described by patches
    1: affine 1-space
    2: affine 1-space
  in the coordinate(s)
    1: [(t//s)]
    2: [(s//t)]

julia> IP2 = projective_space(QQ, [:x, :y, :z]);

julia> S = homogeneous_coordinate_ring(IP2);

julia> x, y, z = gens(S);

julia> IPC, inc_IPC = sub(IP2, ideal(S, [x^2 - y*z]));

julia> C = covered_scheme(IPC);

julia> U = first(affine_charts(IP1))
Spectrum
  of multivariate polynomial ring in 1 variable (t//s)
    over rational field

julia> V = first(affine_charts(C))
Spectrum
  of quotient
    of multivariate polynomial ring in 2 variables (y//x), (z//x)
      over rational field
    by ideal (-(y//x)*(z//x) + 1)

julia> t = first(gens(OO(U)))
(t//s)

julia> Phi = MorphismFromRationalFunctions(IP1, C, U, V, [t//one(t), 1//t]);

julia> realizations = Oscar.realize_on_patch(Phi, U);

julia> realizations[3]
Affine scheme morphism
  from [(t//s)]          AA^1
  to   [(x//z), (y//z)]  scheme((x//z)^2 - (y//z))
given by the pullback function
  (x//z) -> (t//s)
  (y//z) -> (t//s)^2

```
"""
@attributes mutable struct MorphismFromRationalFunctions{DomainType<:AbsCoveredScheme, 
                                       CodomainType<:AbsCoveredScheme
                                      } <: AbsCoveredSchemeMorphism{DomainType, CodomainType, 
                                                                    MorphismFromRationalFunctions, Nothing}
  domain::DomainType
  codomain::CodomainType
  domain_covering::Covering
  codomain_covering::Covering
  domain_chart::AbsAffineScheme
  codomain_chart::AbsAffineScheme
  coord_imgs::Vector{<:FieldElem}
  run_internal_checks::Bool

  ### Various fields for caching
  patch_representatives::IdDict{<:AbsAffineScheme, <:Tuple{<:AbsAffineScheme, <:Vector{<:FieldElem}}}
  realizations::IdDict{<:AbsAffineScheme, <:Vector{<:AffineSchemeMor}}
  realization_previews::IdDict{<:Tuple{<:AbsAffineScheme, <:AbsAffineScheme}, <:Vector{<:FieldElem}}
  maximal_extensions::IdDict{<:Tuple{<:AbsAffineScheme, <:AbsAffineScheme}, <:Vector{<:AbsAffineSchemeMor}}
  cheap_realizations::IdDict{<:Tuple{<:AbsAffineScheme, <:AbsAffineScheme}, <:AbsAffineSchemeMor}
  full_realization::CoveredSchemeMorphism

  function MorphismFromRationalFunctions(
      X::AbsCoveredScheme, Y::AbsCoveredScheme, 
      U::AbsAffineScheme, V::AbsAffineScheme,
      a::Vector{<:FieldElem};
      check::Bool=true,
      domain_covering::Covering=default_covering(X),
      codomain_covering::Covering=default_covering(Y)
    )
    @check is_irreducible(X) "domain must be irreducible"
    @check is_irreducible(Y) "codomain must be irreducible"
    @check dim(Y) <= dim(X) "cannot be dominant"
    #_find_chart(U, default_covering(X)) !== nothing || error("patch not found in domain")
    #_find_chart(V, default_covering(Y)) !== nothing || error("patch not found in codomain")
    any(x->x===U, patches(default_covering(X))) || error("patch not found in domain")
    any(x->x===V, patches(default_covering(Y))) || error("patch not found in codomain")
    F = parent(first(a))
    R = base_ring(F)
    all(x->parent(x)===F, a) || error("coordinate images must be elements of the same field")
    R === ambient_coordinate_ring(U) || error("images of pullback of the coordinates do not live in the correct ring")
    patch_repr = IdDict{AbsAffineScheme, Tuple{AbsAffineScheme, Vector{FieldElem}}}()
    patch_repr[U] = (V, a)
    realizations = IdDict{AbsAffineScheme, Vector{AffineSchemeMor}}()
    realization_previews = IdDict{Tuple{AbsAffineScheme, AbsAffineScheme}, Vector{FieldElem}}()
    maximal_extensions = IdDict{Tuple{AbsAffineScheme, AbsAffineScheme}, Vector{AbsAffineSchemeMor}}()
    cheap_realizations = IdDict{Tuple{AbsAffineScheme, AbsAffineScheme}, AbsAffineSchemeMor}()
    return new{typeof(X), typeof(Y)}(X, Y, domain_covering, codomain_covering, 
                                     U, V, a, check, patch_repr, realizations, 
                                     realization_previews, maximal_extensions,
                                     cheap_realizations
                                    )
  end
end

