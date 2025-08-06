

abstract type AbsCoherentSheaf{
                               SpaceType, OpenType,
                               OutputType, RestrictionType
                              } <: AbsPreSheaf{
                                               SpaceType, OpenType,
                                               OutputType, RestrictionType
                                              }
end

### Interface for coherent sheaves

@doc raw"""
    scheme(F::AbsCoherentSheaf)

Return the scheme on which this sheaf is defined.
"""
scheme(F::AbsCoherentSheaf) = space(underlying_presheaf(F))

@doc raw"""
    sheaf_of_rings(F::AbsCoherentSheaf)

Return the sheaf of rings over which ``ℱ`` is defined.
"""
function sheaf_of_rings(F::AbsCoherentSheaf)
  error("method not implemented for coherent sheaves of type $(typeof(F))")
end

# Manage some left offsets so that the labels are aligned on the right - could
# have alignment issues in the case where we have more than 10 patches to
# describe the restrictions of the sheaf
function Base.show(io::IO, ::MIME"text/plain", M::AbsCoherentSheaf)
  io = pretty(io)
  X = scheme(M)
  cov = default_covering(X)
  print(io, "Coherent sheaf of modules")
  if has_attribute(M, :name)
    print(io, " ", get_attribute(M, :name))
  end
  println(io)
  print(io, Indent(), "on ", Lowercase())
  show(IOContext(io, :show_semi_compact => true, :covering => cov), X)
  if length(cov) > 0
    l = ndigits(length(cov))
    println(io)
    print(io, Dedent(), "with restriction")
    length(cov) > 1 && print(io, "s")
    print(io, Indent())
    for i in 1:length(cov)
      li = ndigits(i)
      U = cov[i]
      println(io)
      print(io, " "^(l-li)*"$i: ", Lowercase(), M(U))
    end
  end
  print(io, Dedent())
end

function Base.show(io::IO, M::AbsCoherentSheaf)
  io = pretty(io)
  if is_terse(io)
    print(io, "Coherent sheaf of modules")
  elseif has_attribute(M, :name)
    print(io, get_attribute(M, :name))
  else
    if is_unicode_allowed()
      print(io, "Coherent sheaf of $(sheaf_of_rings(M))-modules on ", Lowercase(), scheme(M))
    else
      print(io, "Coherent sheaf of modules on ", Lowercase(), scheme(M))
    end
  end
end


### The following provides a function for the internal checks that
# a given set U is open in and admissible for sheaves of modules on X.
#
# We allow the following cases:
#
#  * U::PrincipalOpenSubset in W===ambient_scheme(U) in the basic charts of X
#  * U::PrincipalOpenSubset ⊂ V::PrincipalOpenSubset with ambient_scheme(U) === ambient_scheme(V) in the basic charts of X
#  * U::PrincipalOpenSubset ⊂ V::PrincipalOpenSubset with ambient_scheme(U) != ambient_scheme(V) both in the basic charts of X
#    and U and V contained in the gluing domains of their ambient schemes
#  * U::AbsAffineScheme ⊂ U::AbsAffineScheme in the basic charts of X
#  * U::AbsAffineScheme ⊂ X for U in the basic charts
#  * U::PrincipalOpenSubset ⊂ X with ambient_scheme(U) in the basic charts of X
function _is_open_for_modules(X::AbsCoveredScheme)
  function is_open_func(U::PrincipalOpenSubset, V::PrincipalOpenSubset)
    C = default_covering(X)
    A = ambient_scheme(U)
    A in C || return false
    B = ambient_scheme(V)
    B in C || return false
    if A === B
      is_subset(U, V) || return false
    else
      G = C[A, B] # Get the gluing
      f, g = gluing_morphisms(G)
      is_subset(U, domain(f)) || return false
      gU = preimage(g, U)
      is_subset(gU, V) || return false
    end
    return true
  end
  function is_open_func(U::PrincipalOpenSubset, Y::AbsCoveredScheme)
    return Y === X && ambient_scheme(U) in affine_charts(X)
  end
  function is_open_func(U::AbsAffineScheme, Y::AbsCoveredScheme)
    return Y === X && U in affine_charts(X)
  end
  function is_open_func(Z::AbsCoveredScheme, Y::AbsCoveredScheme)
    return X === Y === Z
  end
  function is_open_func(U::AbsAffineScheme, V::AbsAffineScheme)
    any(x->x===U, affine_charts(X)) || return false
    any(x->x===V, affine_charts(X)) || return false
    G = default_covering(X)[U, V]
    return is_subscheme(U, gluing_domains(G)[1])
  end
  function is_open_func(
      U::AbsAffineScheme,
      V::Union{<:PrincipalOpenSubset, <:SimplifiedAffineScheme}
    )
    is_subscheme(U, V) && return true
    any(x->x===U, affine_charts(X)) || return false
    inc_V_flat = _flatten_open_subscheme(V, default_covering(X))
    A = ambient_scheme(codomain(inc_V_flat))
    Vdirect = codomain(inc_V_flat)
    W = ambient_scheme(Vdirect)
    haskey(gluings(default_covering(X)), (W, U)) || return false # In this case, they are not glued
    G = default_covering(X)[W, U]
    f, g = gluing_morphisms(G)
    pre_V = preimage(g, V)
    return is_subset(U, pre_V)
  end
  function is_open_func(
      U::Union{<:PrincipalOpenSubset, <:SimplifiedAffineScheme},
      V::AbsAffineScheme
    )
    any(x->x===V, affine_charts(X)) || return false
    inc_U_flat = _flatten_open_subscheme(U, default_covering(X))
    A = ambient_scheme(codomain(inc_U_flat))
    Udirect = codomain(inc_U_flat)
    W = ambient_scheme(Udirect)
    haskey(gluings(default_covering(X)), (W, V)) || return false # In this case, they are not glued
    G = default_covering(X)[W, V]
    return is_subset(Udirect, gluing_domains(G)[1])
  end
  return is_open_func
end


########################################################################
# Coherent sheaves of modules on covered schemes                       #
########################################################################
@doc raw"""
    SheafOfModules <: AbsPreSheaf

A sheaf of modules ``ℳ`` on an `AbsCoveredScheme` ``X``.

Note that due to technical reasons, the admissible open subsets are restricted
to the following:
 * `U::AbsAffineScheme` among the `basic_patches` of the `default_covering` of `X`;
 * `U::PrincipalOpenSubset` with `ambient_scheme(U)` in the `basic_patches` of the `default_covering` of `X`.

One can call the restriction maps of ``ℳ`` across charts implicitly using the
identifications given by the gluings in the `default_covering`.
"""
@attributes mutable struct SheafOfModules{SpaceType, OpenType,
                                          OutputType,
                                          RestrictionType
                                         } <: AbsCoherentSheaf{
                                                               SpaceType, OpenType,
                                                               OutputType, RestrictionType
                                                              }
  ID::IdDict{AbsAffineScheme, ModuleFP} # the modules on the basic patches of the default covering
  MG::IdDict{<:Tuple{<:AbsAffineScheme, <:AbsAffineScheme}, <:MatrixElem}; # A dictionary for pairs `(U, V)` of
  # `affine_charts` of `X` such that
  # A = MG[(U, V)] has entries aᵢⱼ with
  # gᵢ = ∑ⱼ aᵢⱼ ⋅ fⱼ on U ∩ V with gᵢ the
  # restrictions of the generators of M[U]
  # and fⱼ the restrictions of the generators
  # of MD[V]. The aᵢⱼ are elements of 𝒪(U ∩ V)
  # represented as a subset of V.
  OOX::StructureSheafOfRings # the structure sheaf on X
  M::PreSheafOnScheme # the underlying presheaf of modules for caching
  C::Covering # The covering of X on which the modules had first been described, a.k.a. the
              # `default_covering` of this sheaf ℱ.

  ### Sheaves of modules on covered schemes
  function SheafOfModules(X::AbsCoveredScheme,
      MD::IdDict{AbsAffineScheme, ModuleFP}, # A dictionary of modules on the `affine_charts` of `X`
      MG::IdDict{Tuple{AbsAffineScheme, AbsAffineScheme}, MatrixElem}; # A dictionary for pairs `(U, V)` of
                                                       # `affine_charts` of `X` such that
                                                       # A = MG[(U, V)] has entries aᵢⱼ with
                                                       # gᵢ = ∑ⱼ aᵢⱼ ⋅ fⱼ on U ∩ V with gᵢ the
                                                       # restrictions of the generators of M[U]
                                                       # and fⱼ the restrictions of the generators
                                                       # of MD[V]. The aᵢⱼ are elements of 𝒪(U ∩ V)
                                                       # represented as a subset of V.
      check::Bool=true,
      default_cov::Covering=begin                      # This will be the `default_covering` of the sheaf to be created.
        patch_list = collect(keys(MD))
        gluing_dict = IdDict{Tuple{AbsAffineScheme, AbsAffineScheme}, AbsGluing}()
        C = Covering(patch_list, gluing_dict)
        inherit_gluings!(C, default_covering(X))
        C
      end
    )
    OOX = OO(X)
    # Make sure that all patches and gluings of the
    # given `default_covering` of the sheaf ℱ to be created
    # are compatible with the data in the dictionaries.
    all(x->haskey(MD, x), patches(default_cov)) || error("all patches in the default covering must have a prescribed module")
    all(x->any(y->(x===y), patches(default_cov)), keys(MD)) || error("all prescribed modules must appear in the default covering")
    all(x->haskey(MG, x), keys(gluings(default_cov))) || error("all gluings in the default covering must have a prescribed transition")
    all(x->any(y->(x===y), keys(gluings(default_cov))), keys(MG)) || error("all prescribed transitions must correspond to gluings in the default covering")

    Mpre = PreSheafOnScheme(X,
                      OpenType=AbsAffineScheme, OutputType=ModuleFP,
                      RestrictionType=Map,
                      is_open_func=_is_open_func_for_schemes_without_affine_scheme_open_subscheme(X)
                      #is_open_func=_is_open_for_modules(X)
                     )
    M = new{typeof(X), AbsAffineScheme, ModuleFP, Map}(MD, MG, OOX, Mpre, default_cov)
    @check begin
      # Check that all sheaves of modules are compatible on the overlaps.
      # TODO: eventually replace by a check that on every basic
      # affine patch, the ideal sheaf can be inferred from what is
      # given on one dense open subset.
      true
    end
    return M
  end
end

### forwarding and implementing the required getters
underlying_presheaf(M::SheafOfModules) = M.M
sheaf_of_rings(M::SheafOfModules) = M.OOX

### Implementing the additional getters
default_covering(M::SheafOfModules) = M.C
restrictions_dict(M::SheafOfModules) = M.ID

@doc raw"""
    twisting_sheaf(IP::AbsProjectiveScheme{<:Field}, d::Int)

For a `ProjectiveScheme` ``ℙ`` return the ``d``-th twisting sheaf
``𝒪(d)`` as a `CoherentSheaf` on ``ℙ``.

# Examples
```jldoctest
julia> P = projective_space(QQ,3)
Projective space of dimension 3
  over rational field
with homogeneous coordinates [s0, s1, s2, s3]

julia> twisting_sheaf(P, 4)
Coherent sheaf of modules
  on scheme over QQ covered with 4 patches
    1: [(s1//s0), (s2//s0), (s3//s0)]   affine 3-space
    2: [(s0//s1), (s2//s1), (s3//s1)]   affine 3-space
    3: [(s0//s2), (s1//s2), (s3//s2)]   affine 3-space
    4: [(s0//s3), (s1//s3), (s2//s3)]   affine 3-space
with restrictions
  1: free module of rank 1 over multivariate polynomial ring in 3 variables over QQ
  2: free module of rank 1 over multivariate polynomial ring in 3 variables over QQ
  3: free module of rank 1 over multivariate polynomial ring in 3 variables over QQ
  4: free module of rank 1 over multivariate polynomial ring in 3 variables over QQ
```
"""
function twisting_sheaf(IP::AbsProjectiveScheme{<:Field}, d::Int)
  # First, look up whether this sheaf has already been computed:
  if !has_attribute(IP, :twisting_sheaves)
    set_attribute!(IP, :twisting_sheaves, Dict{Int, SheafOfModules}())
  end
  twisting_sheaves = get_attribute(IP, :twisting_sheaves)
  haskey(twisting_sheaves, d) && return twisting_sheaves[d]

  X = covered_scheme(IP)
  MD = IdDict{AbsAffineScheme, ModuleFP}()
  S = homogeneous_coordinate_ring(IP)
  n = ngens(S)-1
  for i in 1:n+1
    U = affine_charts(X)[i]
    MD[U] = FreeMod(OO(U), ["$(symbols(S)[i])^$d"])
  end

  MG = IdDict{Tuple{AbsAffineScheme, AbsAffineScheme}, MatrixElem}()
  C = default_covering(X)
  for G in values(gluings(C))
    (U, V) = patches(G)
    (UU, VV) = gluing_domains(G)
    h_U = complement_equation(UU)
    h_V = complement_equation(VV)
    MG[(U, V)] = diagonal_matrix((d>= 0 ? (x->OO(VV)(x, check=false))(h_V^d) : (inv((x->OO(VV)(x, check=false))(h_V))^(-d))), 1)
    MG[(V, U)] = diagonal_matrix((d>= 0 ? (x->OO(UU)(x, check=false))(h_U^d) : (inv((x->OO(UU)(x, check=false))(h_U))^(-d))), 1)
  end

  M = SheafOfModules(X, MD, MG)
  # Cache the result for the next usage
  twisting_sheaves[d] = M
  return M
end


@doc raw"""
    tautological_bundle(IP::AbsProjectiveScheme{<:Field})

For a `ProjectiveScheme` ``ℙ`` return the sheaf ``𝒪(-1)`` as a `CoherentSheaf` on ``ℙ``.

# Examples
```jldoctest
julia> P = projective_space(QQ,3)
Projective space of dimension 3
  over rational field
with homogeneous coordinates [s0, s1, s2, s3]

julia> tautological_bundle(P)
Coherent sheaf of modules
  on scheme over QQ covered with 4 patches
    1: [(s1//s0), (s2//s0), (s3//s0)]   affine 3-space
    2: [(s0//s1), (s2//s1), (s3//s1)]   affine 3-space
    3: [(s0//s2), (s1//s2), (s3//s2)]   affine 3-space
    4: [(s0//s3), (s1//s3), (s2//s3)]   affine 3-space
with restrictions
  1: free module of rank 1 over multivariate polynomial ring in 3 variables over QQ
  2: free module of rank 1 over multivariate polynomial ring in 3 variables over QQ
  3: free module of rank 1 over multivariate polynomial ring in 3 variables over QQ
  4: free module of rank 1 over multivariate polynomial ring in 3 variables over QQ
```
"""
function tautological_bundle(IP::AbsProjectiveScheme{<:Field})
    return twisting_sheaf(IP, -1)
end

@doc raw"""
    cotangent_sheaf(X::AbsCoveredScheme)

For an `AbsCoveredScheme` ``X``, return the sheaf ``Ω¹(X)`` of Kaehler-differentials
on ``X`` as a `CoherentSheaf`.

# Examples
```jldoctest
julia> IP1 = projective_space(QQ, 1);

julia> X = covered_scheme(IP1)
Scheme
  over rational field
with default covering
  described by patches
    1: affine 1-space
    2: affine 1-space
  in the coordinate(s)
    1: [(s1//s0)]
    2: [(s0//s1)]

julia> Omega = cotangent_sheaf(X);

julia> U, V = affine_charts(X);

julia> UV, VU = gluing_domains(default_covering(X)[U, V]);

julia> dx = Omega(U)[1]
d(s1//s0)

julia> Omega(V)
Free module of rank 1 over multivariate polynomial ring in 1 variable over QQ

julia> Omega(U, VU)(dx)
-1/(s0//s1)^2*d(s0//s1)

```
"""
@attr SheafOfModules function cotangent_sheaf(X::AbsCoveredScheme)
  MD = IdDict{AbsAffineScheme, ModuleFP}()
  for U in affine_charts(X)
    MD[U] = cotangent_module(U)
  end
  MG = IdDict{Tuple{AbsAffineScheme, AbsAffineScheme}, MatrixElem}()
  C = default_covering(X)
  for G in values(gluings(C))
    (U, V) = patches(G)
    (UU, VV) = gluing_domains(G)
    (f, g) = gluing_morphisms(G)
    MG[(U, V)] = transpose(jacobian_matrix(pullback(g).(gens(OO(UU)))))
    MG[(V, U)] = transpose(jacobian_matrix(pullback(f).(gens(OO(VV)))))
  end


  M = SheafOfModules(X, MD, MG)
  return M
end

@doc raw"""
    cotangent_module(X::AbsAffineScheme)

Return the ``𝒪(X)``-module ``Ω¹(X)`` of Kaehler-differentials on ``X``.
"""
function cotangent_module(X::AbsAffineScheme)
  error("method not implemented for this type of ring")
end

@attr ModuleFP function cotangent_module(X::AbsAffineScheme{<:Field, <:MPolyRing})
  R = OO(X)
  F = FreeMod(R, ["d$(x)" for x in symbols(R)])
  return F
end

@attr ModuleFP function cotangent_module(X::AbsAffineScheme{<:Field, <:MPolyLocRing})
  R = OO(X)
  P = base_ring(R)
  F = FreeMod(R, ["d$(x)" for x in symbols(P)])
  return F
end

@attr ModuleFP function cotangent_module(X::AbsAffineScheme{<:Field, <:MPolyQuoRing})
  R = OO(X)
  P = base_ring(R)
  F = FreeMod(R, ["d$(x)" for x in symbols(P)])
  rels, _ = sub(F, transpose(change_base_ring(R, jacobian_matrix(base_ring(modulus(R)), gens(modulus(R))))))
  M, _ = quo(F, rels)
  return M
end

@attr ModuleFP function cotangent_module(X::AbsAffineScheme{<:Field, <:MPolyQuoLocRing})
  R = OO(X)
  P = base_ring(R)
  F = FreeMod(R, ["d$(x)" for x in symbols(P)])
  rels, _ = sub(F, transpose(change_base_ring(R, jacobian_matrix(base_ring(modulus(R)), gens(modulus(R))))))
  M, _ = quo(F, rels)
  return M
end

########################################################################
# Hom-Sheaves                                                          #
########################################################################
#=
  Hom sheaves ℋom(ℱ, 𝒢) are special.

  First of all, they can be made completely lazy, as their modules
  on U ⊂ X can be created from ℱ(U) and 𝒢(U) on the fly and the same
  holds for their transition- and restriction functions.

  Second, Hom sheaves come with a domain, a codomain, and an
  interpretation mapping of their sections as homomorphisms.

  We realize hom sheaves in this way, taking only ℱ and 𝒢 as input
  in the constructor.
=#

@doc raw"""
    HomSheaf
    
For two `AbsCoherentSheaf`s `F` and `G` on an `AbsCoveredScheme` `X` 
this computes the sheaf associated to `U -> Hom(F(U), G(U))`.
# Examples
```jldoctest
julia> IP1 = projective_space(GF(7), [:x, :y])
Projective space of dimension 1
  over prime field of characteristic 7
with homogeneous coordinates [x, y]

julia> Y = covered_scheme(IP1);

julia> Omega = cotangent_sheaf(Y)
Coherent sheaf of modules
  on scheme over GF(7) covered with 2 patches
    1: [(y//x)]   affine 1-space
    2: [(x//y)]   affine 1-space
with restrictions
  1: free module of rank 1 over multivariate polynomial ring in 1 variable over GF(7)
  2: free module of rank 1 over multivariate polynomial ring in 1 variable over GF(7)

julia> F = free_module(OO(Y), 1)
Coherent sheaf of modules
  on scheme over GF(7) covered with 2 patches
    1: [(y//x)]   affine 1-space
    2: [(x//y)]   affine 1-space
with restrictions
  1: free module of rank 1 over multivariate polynomial ring in 1 variable over GF(7)
  2: free module of rank 1 over multivariate polynomial ring in 1 variable over GF(7)

julia> T = Oscar.HomSheaf(Omega, F)
Coherent sheaf of modules
  on scheme over GF(7) covered with 2 patches
    1: [(y//x)]   affine 1-space
    2: [(x//y)]   affine 1-space
with restrictions
  1: hom of (Multivariate polynomial ring in 1 variable over GF(7)^1, Multivariate polynomial ring in 1 variable over GF(7)^1)
  2: hom of (Multivariate polynomial ring in 1 variable over GF(7)^1, Multivariate polynomial ring in 1 variable over GF(7)^1)

julia> typeof(T)
Oscar.HomSheaf{CoveredScheme{FqField}, AbsAffineScheme, ModuleFP, Map}


```
"""
@attributes mutable struct HomSheaf{SpaceType, OpenType, OutputType,
                                    RestrictionType
                                   } <: AbsCoherentSheaf{
                                                         SpaceType, OpenType,
                                                         OutputType, RestrictionType
                                                        }
  domain::AbsCoherentSheaf{SpaceType, OpenType, OutputType, RestrictionType}
  codomain::AbsCoherentSheaf{SpaceType, OpenType, OutputType, RestrictionType}
  OOX::StructureSheafOfRings
  M::PreSheafOnScheme

  function HomSheaf(F::AbsCoherentSheaf, G::AbsCoherentSheaf)
    X = scheme(F)
    X === scheme(G) || error("sheaves must be defined over the same scheme")
    OOX = sheaf_of_rings(F)
    OOX === sheaf_of_rings(G) || error("sheaves must be defined over the same sheaves of rings")

    Mpre = PreSheafOnScheme(X,
                      OpenType=AbsAffineScheme, OutputType=ModuleFP,
                      RestrictionType=Map,
                      is_open_func=_is_open_func_for_schemes_without_affine_scheme_open_subscheme(X)
                     )
    M = new{typeof(X), AbsAffineScheme, ModuleFP, Map}(F, G, OOX, Mpre)

    return M
  end
end

### forwarding and implementation of the essential getters
underlying_presheaf(M::HomSheaf) = M.M
sheaf_of_rings(M::HomSheaf) = M.OOX
domain(M::HomSheaf) = M.domain
codomain(M::HomSheaf) = M.codomain
#default_covering(M::HomSheaf) = default_covering(domain(M)) # TODO: This is only a temporary fix!

########################################################################
# Sheaves of direct sums                                               #
########################################################################
@doc raw"""
    DirectSumSheaf

Given two or more `AbsCoherentSheaf`s `F` and `G` on an `AbsCoveredScheme` `X`, 
this holds the sheaf associated to the direct sum of `F` and `G`

# Examples
```jldoctest
julia> IP1 = projective_space(GF(7), [:x, :y])
Projective space of dimension 1
  over prime field of characteristic 7
with homogeneous coordinates [x, y]

julia> Y = covered_scheme(IP1);

julia> Omega = cotangent_sheaf(Y)
Coherent sheaf of modules
  on scheme over GF(7) covered with 2 patches
    1: [(y//x)]   affine 1-space
    2: [(x//y)]   affine 1-space
with restrictions
  1: free module of rank 1 over multivariate polynomial ring in 1 variable over GF(7)
  2: free module of rank 1 over multivariate polynomial ring in 1 variable over GF(7)

julia> F = free_module(OO(Y), 1)
Coherent sheaf of modules
  on scheme over GF(7) covered with 2 patches
    1: [(y//x)]   affine 1-space
    2: [(x//y)]   affine 1-space
with restrictions
  1: free module of rank 1 over multivariate polynomial ring in 1 variable over GF(7)
  2: free module of rank 1 over multivariate polynomial ring in 1 variable over GF(7)

julia> W = Oscar.DirectSumSheaf(Y, [Omega, F])
Coherent sheaf of modules
  on scheme over GF(7) covered with 2 patches
    1: [(y//x)]   affine 1-space
    2: [(x//y)]   affine 1-space
with restrictions
  1: direct sum of FreeMod{FqMPolyRingElem}[Multivariate polynomial ring in 1 variable over GF(7)^1, Multivariate polynomial ring in 1 variable over GF(7)^1]
  2: direct sum of FreeMod{FqMPolyRingElem}[Multivariate polynomial ring in 1 variable over GF(7)^1, Multivariate polynomial ring in 1 variable over GF(7)^1]

julia> typeof(W)
DirectSumSheaf{CoveredScheme{FqField}, AbsAffineScheme, ModuleFP, Map}

```
"""
@attributes mutable struct DirectSumSheaf{SpaceType, OpenType, OutputType,
                                          RestrictionType
                                         } <: AbsCoherentSheaf{
                                                               SpaceType, OpenType,
                                                               OutputType, RestrictionType
                                                              }
  summands::Vector{AbsCoherentSheaf{SpaceType, OpenType, OutputType, RestrictionType}}
  #projections::Vector #TODO: Realize as sections in HomSheafs
  #inclusions::Vector # TODO: same.
  OOX::StructureSheafOfRings
  M::PreSheafOnScheme

  function DirectSumSheaf(X::AbsCoveredScheme, summands::Vector{<:AbsCoherentSheaf})
    all(x->(scheme(x)===X), summands) || error("summands must be defined over the same scheme")
    OOX = OO(X)
    all(x->(sheaf_of_rings(x)===OOX), summands) || error("summands must be defined over the same sheaves of rings")

    Mpre = PreSheafOnScheme(X, 
                      OpenType=AbsAffineScheme, OutputType=ModuleFP,
                      RestrictionType=Map,
                      is_open_func=_is_open_func_for_schemes_without_affine_scheme_open_subscheme(X)
                     )
    M = new{typeof(X), AbsAffineScheme, ModuleFP, Map}(summands, OOX, Mpre)

    return M
  end
end

### forwarding and implementation of the essential getters
underlying_presheaf(M::DirectSumSheaf) = M.M
sheaf_of_rings(M::DirectSumSheaf) = M.OOX
summands(M::DirectSumSheaf) = M.summands

### user facing constructors
function direct_sum(summands::Vector{<:AbsCoherentSheaf})
  length(summands) == 0 && error("list of summands must not be empty")
  X = scheme(first(summands))
  return DirectSumSheaf(X, summands)
end

function Base.show(io::IO, M::DirectSumSheaf)
  if is_terse(io)
    print(io, "Direct sum of sheaves")
  else
    s = summands(M)
    if is_unicode_allowed() && length(s) > 0
      for i in 1:length(M) - 1
        print(io, "$(s[i]) ⊕ ")
      end
      print(io, "$(s[end])")
    else
      print(io, "Direct sum of sheaves of modules on covered scheme")
    end
  end
end

@doc raw"""
    free_module(R::StructureSheafOfRings, n::Int)

Return the sheaf of free ``𝒪``-modules ``𝒪ⁿ`` for a structure
sheaf of rings ``𝒪 = R``.
"""
function free_module(R::StructureSheafOfRings, n::Int)
  return free_module(R, ["e_$i" for i in 1:n])
end

function free_module(R::StructureSheafOfRings, gen_names::Vector{String})
  return free_module(R, Symbol.(gen_names))
end

function free_module(R::StructureSheafOfRings, gen_names::Vector{Symbol})
  X = space(R)
  n = length(gen_names)
  MD = IdDict{AbsAffineScheme, ModuleFP}()
  for U in affine_charts(X)
    MD[U] = FreeMod(OO(U), gen_names)
  end

  MG = IdDict{Tuple{AbsAffineScheme, AbsAffineScheme}, MatrixElem}()
  C = default_covering(X)
  for G in values(gluings(C))
    (U, V) = patches(G)
    (UU, VV) = gluing_domains(G)
    MG[(U, V)] = identity_matrix(OO(VV), n)
    MG[(V, U)] = identity_matrix(OO(UU), n)
  end

  M = SheafOfModules(X, MD, MG)
  return M
end

@doc raw"""
    dual(M::SheafOfModules)

For a `SheafOfModules` ``ℳ`` on an `AbsCoveredScheme` ``X``, return
the ``𝒪_X``-dual ``ℋ om_{𝒪_X}(ℳ , 𝒪_X)`` of ``ℳ``.
"""
@attr AbsCoherentSheaf function dual(M::SheafOfModules)
  OOX = sheaf_of_rings(M)
  F = free_module(OOX, ["1"])
  return HomSheaf(M, F)
end

@doc raw"""
    tangent_sheaf(X::AbsCoveredScheme)

Return the tangent sheaf ``T_X`` of `X`, constructed as ``ℋ om_{𝒪_X}(Ω¹_X, 𝒪_X)``.
"""
@attr HomSheaf function tangent_sheaf(X::AbsCoveredScheme)
  return dual(cotangent_sheaf(X))
end

########################################################################
# Pushforwards of sheaves for closed embeddings                        #
########################################################################

#=
# Let f : X ↪ Y be a closed embedding with ideal sheaf ℐ on Y and ℳ
# a sheaf of modules on X. For an open set U ⊂ Y we have
# f_* ℳ (U) to be the 𝒪_X(f⁻¹(U))-module ℳ (f⁻¹(U)), but seen as
# an 𝒪_Y(U)-module via the natural restriction of functions.
#
# Mathematically, this is almost an implicit operation. But since we
# do not have natural bi-module structures, we need to set up a new
# sheaf of modules on Y, together with canonical identifications
# with the modules on X.
#
# It is clear that this can and should be made lazy.
#                                                                     =#

@doc raw"""
    PushforwardSheaf
    
For a `CoveredClosedEmbedding` `i : X -> Y` and an `AbsCoherentSheaf` `F`
on `X` this computes the coherent sheaf `i_* F` on `Y`.

# Examples
```jldoctest
julia> IP2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> S = cox_ring(IP2)
Multivariate polynomial ring in 3 variables over QQ graded by
  x1 -> [1]
  x2 -> [1]
  x3 -> [1]

julia> x, y, z = gens(S);

julia> I = ideal(S, x^3 - y*z^2);

julia> II = ideal_sheaf(IP2, I);

julia> X, inc_X = sub(II);

julia> F = cotangent_sheaf(X)
Coherent sheaf of modules
  on scheme over QQ covered with 3 patches
    1: [x_1_1, x_2_1]   scheme(x_1_1^3 - x_2_1)
    2: [x_1_2, x_2_2]   scheme(x_1_2^2*x_2_2 - 1)
    3: [x_1_3, x_2_3]   scheme(x_1_3^3 - x_2_3^2)
with restrictions
  1: subquotient of submodule with 2 generators
    1: dx_1_1
    2: dx_2_1
  by submodule with 1 generator
    1: 3*x_1_1^2*dx_1_1 - dx_2_1
  2: subquotient of submodule with 2 generators
    1: dx_1_2
    2: dx_2_2
  by submodule with 1 generator
    1: 2*x_1_2*x_2_2*dx_1_2 + x_1_2^2*dx_2_2
  3: subquotient of submodule with 2 generators
    1: dx_1_3
    2: dx_2_3
  by submodule with 1 generator
    1: 3*x_1_3^2*dx_1_3 - 2*x_2_3*dx_2_3

julia> inc_F = Oscar.PushforwardSheaf(inc_X, F)
Coherent sheaf of modules
  on normal toric variety
with restrictions
  1: subquotient of submodule with 2 generators
    1: dx_1_1
    2: dx_2_1
  by submodule with 5 generators
    1: (x_1_1^3 - x_2_1)*dx_1_1
    2: (x_1_1^3 - x_2_1)*dx_2_1
    3: 3*x_1_1^2*dx_1_1 - dx_2_1
    4: (x_1_1^3 - x_2_1)*dx_1_1
    5: (x_1_1^3 - x_2_1)*dx_2_1
  2: subquotient of submodule with 2 generators
    1: dx_1_2
    2: dx_2_2
  by submodule with 5 generators
    1: (x_1_2^2*x_2_2 - 1)*dx_1_2
    2: (x_1_2^2*x_2_2 - 1)*dx_2_2
    3: 2*x_1_2*x_2_2*dx_1_2 + x_1_2^2*dx_2_2
    4: (x_1_2^2*x_2_2 - 1)*dx_1_2
    5: (x_1_2^2*x_2_2 - 1)*dx_2_2
  3: subquotient of submodule with 2 generators
    1: dx_1_3
    2: dx_2_3
  by submodule with 5 generators
    1: (x_1_3^3 - x_2_3^2)*dx_1_3
    2: (x_1_3^3 - x_2_3^2)*dx_2_3
    3: 3*x_1_3^2*dx_1_3 - 2*x_2_3*dx_2_3
    4: (x_1_3^3 - x_2_3^2)*dx_1_3
    5: (x_1_3^3 - x_2_3^2)*dx_2_3

julia> typeof(inc_F)
PushforwardSheaf{NormalToricVariety, AbsAffineScheme, ModuleFP, Map}

```
"""
@attributes mutable struct PushforwardSheaf{SpaceType, OpenType, OutputType,
                                            RestrictionType
                                           } <: AbsCoherentSheaf{
                                                                 SpaceType, OpenType,
                                                                 OutputType, RestrictionType
                                                                }
  inc::CoveredClosedEmbedding
  OOX::StructureSheafOfRings
  OOY::StructureSheafOfRings
  M::AbsCoherentSheaf
  ident::IdDict{AbsAffineScheme, Union{Map, Nothing}} # a dictionary caching the identifications
  F::PreSheafOnScheme

  function PushforwardSheaf(inc::CoveredClosedEmbedding, M::AbsCoherentSheaf)
    X = domain(inc)
    X === scheme(M) || error("sheaf must be defined over the domain of the embedding")
    OOX = sheaf_of_rings(M)
    Y = codomain(inc)
    OOY = OO(Y)

    ident = IdDict{AbsAffineScheme, Union{Map, Nothing}}()

    Blubber = PreSheafOnScheme(Y, 
                      OpenType=AbsAffineScheme, OutputType=ModuleFP,
                      RestrictionType=Map,
                      is_open_func=_is_open_func_for_schemes_without_affine_scheme_open_subscheme(Y)
                      #is_open_func=_is_open_for_modules(Y)
                     )
    MY = new{typeof(Y), AbsAffineScheme, ModuleFP, Map}(inc, OOX, OOY, M, ident, Blubber)
    return MY
  end
end

### forwarding and implementing the required getters
underlying_presheaf(M::PushforwardSheaf) = M.F
sheaf_of_rings(M::PushforwardSheaf) = M.OOY
original_sheaf(M::PushforwardSheaf) = M.M
morphism(M::PushforwardSheaf) = M.inc
identification_dict(M::PushforwardSheaf) = M.ident

function Base.show(io::IO, M::PushforwardSheaf)
  print(io, "pushforward of $(original_sheaf(M)) along $(morphism(M))")
end

########################################################################
# Pullback of sheaves along morphisms                                  #
########################################################################

#=
# Let f : X → Y be a morphism and ℳ
# a sheaf of modules on Y. For an open set U ⊂ X we have
# f^* ℳ (U) to be the 𝒪_X(U)-module
#
#   𝒪_X ⊗_{f⁻¹𝒪_Y} f⁻¹ℳ
#
# where f⁻¹(ℱ) denotes the sheaf associated to U ↦ lim_{V ⊃ f(U)} ℱ(V).
# On the algebraic side, this merely means carrying out a change of bases
# for the module ℳ (V) where V is some affine open containing f(U).
# To find the latter might be a subtle task for general morphisms of
# schemes. In fact, f will in general only be given with respect to
# fixed coverings CX of X and CY of Y. Since the pullback of sheaves
# is a local question on X, we need to restrict to sufficiently small
# neighborhoods such that
#
#   fᵢ : Uᵢ → Vᵢ
#
# is a local affine representative of the map f. But then the Uᵢ might
# not be `affine_charts` of X, anymore. Thus, we can a priori only
# construct the modules locally on X and `produce_object` has
# to take care of extending them to the `affine_charts` if necessary.
#
# Again, it is clear that this can and should be made lazy.
#                                                                     =#
#
# TODO: PullbackSheaf is not yet fully functional.
#
# Missing parts:
#
#  - If ℳ  is given only on the patches of a refinement Vⱼ, j ∈ J of
#    the `default_covering` of X, then there is no method implemented
#    to create a module for ℳ (U) when U ⊂ X is an `affine_chart` of X.
#    The user is hence forced to work in the refinement only.


@doc raw"""
    PullbackSheaf
    
For a morphism `f : X -> Y` of `AbsCoveredScheme`s and a coherent 
sheaf `F` on `Y` this computes the pullback `f^* F` on `X`.
# Examples
```jldoctest
julia> IP2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> S = cox_ring(IP2)
Multivariate polynomial ring in 3 variables over QQ graded by
  x1 -> [1]
  x2 -> [1]
  x3 -> [1]

julia> x, y, z = gens(S);

julia> I = ideal(S, x^3 - y*z^2);

julia> II = ideal_sheaf(IP2, I);

julia> X, inc_X = sub(II);

julia> F = cotangent_sheaf(codomain(inc_X))
Coherent sheaf of modules
  on normal toric variety
with restrictions
  1: submodule with 2 generators
    1: dx_1_1
    2: dx_2_1
  represented as subquotient with no relations
  2: submodule with 2 generators
    1: dx_1_2
    2: dx_2_2
  represented as subquotient with no relations
  3: submodule with 2 generators
    1: dx_1_3
    2: dx_2_3
  represented as subquotient with no relations

julia> inc_F = Oscar.PullbackSheaf(inc_X, F)
Coherent sheaf of modules
  on scheme over QQ covered with 3 patches
    1: [x_1_1, x_2_1]   scheme(x_1_1^3 - x_2_1)
    2: [x_1_2, x_2_2]   scheme(x_1_2^2*x_2_2 - 1)
    3: [x_1_3, x_2_3]   scheme(x_1_3^3 - x_2_3^2)
with restrictions
  1: submodule with 2 generators
    1: dx_1_1
    2: dx_2_1
  represented as subquotient with no relations
  2: submodule with 2 generators
    1: dx_1_2
    2: dx_2_2
  represented as subquotient with no relations
  3: submodule with 2 generators
    1: dx_1_3
    2: dx_2_3
  represented as subquotient with no relations

julia> typeof(inc_F)
PullbackSheaf{CoveredScheme{QQField}, AbsAffineScheme, ModuleFP, Map}

```
"""
@attributes mutable struct PullbackSheaf{SpaceType, OpenType, OutputType,
                                         RestrictionType
                                        } <: AbsCoherentSheaf{
                                                              SpaceType, OpenType,
                                                              OutputType, RestrictionType
                                                             }
  f::AbsCoveredSchemeMorphism
  OOX::StructureSheafOfRings # the sheaf of rings in the domain
  OOY::StructureSheafOfRings # the sheaf of rings in the codomain
  M::AbsCoherentSheaf        # the sheaf of modules on Y
  pullback_of_sections::IdDict{AbsAffineScheme, Union{Map, Nothing}} # a dictionary caching the natural
                                                                   # pullback maps along the maps in the `covering_morphism` of f
  F::PreSheafOnScheme        # the internal caching instance doing the bookkeeping

  function PullbackSheaf(f::AbsCoveredSchemeMorphism, M::AbsCoherentSheaf)
    X = domain(f)
    Y = codomain(f)
    Y === scheme(M) || error("sheaf must be defined over the domain of the embedding")
    OOY = sheaf_of_rings(M)
    OOX = OO(X)
    fcov = covering_morphism(f)::CoveringMorphism
    CX = domain(fcov)::Covering
    CY = codomain(fcov)::Covering
    pullbacks = IdDict{AbsAffineScheme, Map}()

    ident = IdDict{AbsAffineScheme, Union{Map, Nothing}}()

    Blubber = PreSheafOnScheme(X, 
                      OpenType=AbsAffineScheme, OutputType=ModuleFP,
                      RestrictionType=Map,
                      is_open_func=_is_open_func_for_schemes_without_affine_scheme_open_subscheme(X)
                     )
    MY = new{typeof(X), AbsAffineScheme, ModuleFP, Map}(f, OOX, OOY, M, pullbacks, Blubber)
    return MY
  end
end

underlying_presheaf(M::PullbackSheaf) = M.F
sheaf_of_rings(M::PullbackSheaf) = M.OOX
original_sheaf(M::PullbackSheaf) = M.M
morphism(M::PullbackSheaf) = M.f
pullbacks_on_patches(M::PullbackSheaf) = M.pullback_of_sections

function Base.show(io::IO, M::PullbackSheaf)
  print(io, "pullback of $(original_sheaf(M)) along $(morphism(M))")
end


########################################################################
# pushforward of modules                                               #
########################################################################
#
# It is assumed that f : R → S is a map of rings such that S ≅ R/I and
# M is an S-module. We transform M into an R-module by adding the
# necessary relations. The return value is that new module M', together
# with its identification map M' → M. Note that we can not give the
# inverse of this map, since there is no well-defined underlying ring
# homomorphism.
function _pushforward(f::Map{<:Ring, <:Ring}, I::Ideal, M::FreeMod)
  R = domain(f)
  S = codomain(f)
  base_ring(I) === R || error("ideal is not defined over the correct ring")
  base_ring(M) === S || error("module is not defined over the correct ring")
  FR = FreeMod(R, M.S) # M.S are the symbols of M
  #FRtoM = hom(FR, M, gens(M), f)
  MR, res = quo(FR, (I*FR)[1])
  ident = hom(MR, M, gens(M), f)
  return MR, ident
end

function _pushforward(f::Map{<:Ring, <:Ring}, I::Ideal, M::SubquoModule)
  R = domain(f)
  S = codomain(f)
  base_ring(I) === R || error("ideal is not defined over the correct ring")
  base_ring(M) === S || error("module is not defined over the correct ring")
  FS = ambient_free_module(M)
  Mgens = ambient_representatives_generators(M)
  Mrels = relations(M)
  FR, identF = _pushforward(f, I, FS)
  MRgens = [preimage(identF, v) for v in ambient_representatives_generators(M)]
  MRrels = elem_type(FR)[preimage(identF, w) for w in relations(M)]
  MRrels_ext = vcat(MRrels, elem_type(FR)[g*e for g in gens(I) for e in gens(FR)])
  MR = quo(sub(FR, MRgens)[1], sub(FR, MRrels_ext)[1])[1]
  ident = hom(MR, M, gens(M), f)
  return MR, ident
end

@attr Bool function is_locally_free(M::AbsCoherentSheaf)
  return all(U->is_projective(M(U))[1], affine_charts(scheme(M)))
end

#@attr Covering function trivializing_covering(M::AbsCoherentSheaf)
@attr Covering function trivializing_covering(M::AbsCoherentSheaf)
  X = scheme(M)
  OOX = OO(X)
  patch_list = Vector{AbsAffineScheme}()
  for U in affine_charts(X)
    patch_list = vcat(patch_list, _trivializing_covering(M, U))
  end
  C = Covering(patch_list)
  inherit_gluings!(C, default_covering(X))
  if has_decomposition_info(default_covering(X))
    for U in patches(C)
      V = __find_chart(U, default_covering(X))
      phi = OOX(V, U)
      new_info = Vector{elem_type(OO(U))}(phi.(decomposition_info(default_covering(X))[V]))
      set_decomposition_info!(C, U, new_info)
    end
  end

  push!(coverings(X), C)
  return C
end

@attr Covering function trivializing_covering(M::HomSheaf)
  X = scheme(M)
  OOX = OO(X)
  # The problem is that every module of a HomSheaf must know that it is
  # a hom-module. Hence, the way to go is to pass through a common
  # refinement of domain and codomain and recreate all the hom modules
  # as free modules on this covering.
  #
  # But for this, it is not yet clear where to locate the patches of these
  # refinements in the tree and how to deal with the restriction maps in a
  # clean way. Say M = Hom(F, G) where F is trivialized on {Uᵢ} and G on
  # {Vⱼ}. Then W = Uᵢ∩ Vⱼ would have to be a PrincipalOpenSubset of both
  # Uᵢ and Vⱼ for the restrictions of F and G to induce the proper job
  # on restrictions to M(W) automatically.
  # Hence, we need to manually prescribe how to trivialize and restrict
  # M on the Ws.
  dom_triv = trivializing_covering(domain(M))
  cod_triv = trivializing_covering(codomain(M))
  patch_list = AbsAffineScheme[]
  for U in patches(dom_triv)
    for V in patches(cod_triv)
      success, W = _have_common_ancestor(U, V)
      if success
        incU = _flatten_open_subscheme(U, W)
        incV = _flatten_open_subscheme(V, W)
        UV = intersect(codomain(incU), codomain(incV))::PrincipalOpenSubset
        push!(patch_list, UV)

        dom_UV, dom_res = change_base_ring(OOX(U, UV), domain(M)(U))
        add_incoming_restriction!(domain(M), U, dom_UV, dom_res)
        add_incoming_restriction!(domain(M), W, dom_UV,
                                  compose(domain(M)(W, U), dom_res))
        object_cache(domain(M))[UV] = dom_UV

        cod_UV, cod_res = change_base_ring(OOX(V, UV), codomain(M)(V))
        add_incoming_restriction!(codomain(M), V, cod_UV, cod_res)
        add_incoming_restriction!(codomain(M), W, cod_UV,
                                  compose(codomain(M)(W, V), cod_res))
        object_cache(codomain(M))[UV] = cod_UV

        MUV = M(UV) # This will be a free module; we need to prescribe the restrictions!
        MW = M(W)
        img_gens = elem_type(MUV)[]
        # every generator g of MW is a homomorphism. It takes an element
        # v ∈ domain(M)(W) to w = ϕ_{g}(v) ∈ codomain(M)(W).
        # Where does g map to when restricting to MUV?
        #
        for g in gens(MW)
          phi = element_to_homomorphism(g)
          img_gens_phi = cod_res.(codomain(M)(W, V).(phi.(gens(domain(M)(W)))))
          sub_dom, inc_dom = sub(domain(M)(UV),
                                 domain(M)(W, UV).(gens(domain(M)(W))))
                                 #dom_res.(domain(M)(W, U).(gens(domain(M)(W)))))
          img_gens_psi = elem_type(codomain(M)(UV))[]
          for v in gens(domain(M)(UV))
            w = preimage(inc_dom, v)
            c = coordinates(w) # These are the coordinates in the original set
                               # of generators in the domain
            # We use this to compute the image of v
            phi_v = sum([c[i]*img_gens_phi[i] for i in 1:length(img_gens_phi)], init=zero(codomain(M)(UV)))
            # and push it to the list.
            push!(img_gens_psi, phi_v)
          end
          # From that list, we can assemble what the restriction of phi
          # looks like as a homomorphism
          psi = hom(domain(M)(UV), codomain(M)(UV), img_gens_psi)
          # and convert it to a module element.
          img_g = homomorphism_to_element(MUV, psi)
          push!(img_gens, img_g)
        end

        # Finally, this allows us to assemble the restriction map
        res = hom(MW, MUV, img_gens, OOX(W, UV))
        add_incoming_restriction!(M, W, MUV, res)
      end
    end
  end
  C = Covering(patch_list)
  inherit_gluings!(C, default_covering(scheme(M)))
  push!(coverings(X), C)
  return C
end

function _trivializing_covering(M::AbsCoherentSheaf, U::AbsAffineScheme)
  X = scheme(M)
  OOX = OO(X)
  MU = M(U)
  MU isa FreeMod && return [U]
  MU::SubquoModule
  A = _presentation_matrix(MU)
  if iszero(A)
    # Trivial shortcut in the recursion.
    # We nevertheless need to recreate U as a PrincipalOpenSubset of itself
    # as we are not allowed to alter the values of the sheaf M on U directly.
    V = PrincipalOpenSubset(U, one(OO(U)))
    F = FreeMod(OO(V), ncols(A))
    res = hom(MU, F, gens(F), OOX(U, V))
    add_incoming_restriction!(M, U, F, res)
    object_cache(M)[V] = F
    return [V]
  end

  # We do not need to go through all entries of A, but only those
  # necessary to generate the unit ideal.
  I = ideal(OOX(U), [A[i, j] for i in 1:nrows(A) for j in 1:ncols(A)])
  if !(one(OOX(U)) in I)
    # Now two things could be happening.
    # 1. The sheaf is not locally trivial.
    # 2. We might have different disjoint components on which
    #    the sheaf has different ranks.
    Y = connected_components(U)
    length(Y) == 1 && error("sheaf is not locally free")

    return_patches = AbsAffineScheme[]
    for V in Y
      # Test for being locally free on V
      rho = OOX(U, V)
      MV, res = change_base_ring(rho, MU)
      add_incoming_restriction!(M, U, MV, res)
      object_cache(M)[V] = MV
      return_patches = vcat(return_patches, _trivializing_covering(M, V))
    end
    return return_patches
  end

  # The non-zero coordinates provide us with a list of entries which
  # are sufficient to do so. This set can not assumed to be minimal, though.
  a = coordinates(one(OOX(U)), I)
  nonzero_entries = [ i for i in 1:ngens(I) if !iszero(a[i])]
  return_patches = AbsAffineScheme[]

  for t in nonzero_entries
    i = div(t-1, ncols(A)) + 1
    j = mod(t-1, ncols(A)) + 1 # The matrix coordinates of the nonzero entry
    # We invert the (i,j)-th entry of A.
    # Then we can reduce the presentation matrix so that we can throw away one
    # of the generators of the module.
    V = PrincipalOpenSubset(U, A[i, j])
    Ares = map_entries(OOX(U, V), A)
    uinv = inv(Ares[i, j])
    multiply_row!(Ares, uinv, i)
    for k in 1:i-1
      #multiply_row!(Ares, u, k)
      add_row!(Ares, -A[k, j], i, k)
    end
    for k in i+1:nrows(Ares)
      #multiply_row!(Ares, u, k)
      add_row!(Ares, -A[k, j], i, k)
    end
    Asub = Ares[[k for k in 1:nrows(Ares) if k != i], [k for k in 1:ncols(Ares) if k !=j]]

    # Assemble the restriction map from the parent node
    if iszero(Asub)
      # End of recursion.
      # Create a free module and the corresponding restriction morphism.
      F = FreeMod(OO(V), ncols(Asub))
      img_gens = elem_type(F)[]
      for k in 1:j-1
        push!(img_gens, F[k])
      end
      push!(img_gens,
            -sum([Ares[i, k]*F[(k>j ? k-1 : k)] for k in 1:ncols(Ares) if k!=j],
                 init=zero(F))
           )
      for k in j+1:ncols(Ares)
        push!(img_gens, F[k-1])
      end
      res = hom(MU, F, img_gens, OOX(U, V))

      # Since we are messing with the internals of the sheaf, we need
      # to leave everything clean. This includes manual caching.
      add_incoming_restriction!(M, U, F, res)
      object_cache(M)[V] = F
      set_attribute!(F, :_presentation_matrix, Asub)
      push!(return_patches, V)
    else
      # Intermediate recursion step.
      # Recreate the restriction of the module to the open subset but with one generator
      # less and construct the restriction map.
      F, amb_res = change_base_ring(OOX(U, V), ambient_free_module(MU))
      v = ambient_representatives_generators(MU)
      M_gens = amb_res.(v)
      rest_gens = [M_gens[k] for k in 1:length(M_gens) if k!=j]
      rels = [amb_res(w) for w in relations(MU)]
      MV = SubquoModule(F, rest_gens, rels)
      img_gens = elem_type(F)[]
      for k in 1:j-1
        push!(img_gens, M_gens[k])
      end
      push!(img_gens,
            -sum([Ares[i, k]*M_gens[k] for k in 1:length(M_gens) if k!=j],
                 init=zero(F))
           )
      for k in j+1:length(M_gens)
        push!(img_gens, M_gens[k])
      end
      res = hom(MU, MV, MV.(img_gens), OOX(U, V))
      add_incoming_restriction!(M, U, MV, res)
      object_cache(M)[V] = MV
      set_attribute!(MV, :_presentation_matrix, Asub)
      return_patches = vcat(return_patches, _trivializing_covering(M, V))
    end
  end
  return return_patches
end

@attr MatrixElem function _presentation_matrix(M::ModuleFP)
  return matrix(map(presentation(M), 1))
end

########################################################################
# Projectivization of vector bundles                                   #
#
# To any vector bundle E → X one can associate its projectivization
# ℙ(E) → X. In general, E will be given as a locally free
# AbsCoherentSheaf.
########################################################################

@doc raw"""
    projectivization(E::AbsCoherentSheaf;
        var_names::Vector{String}=Vector{String}(),
        check::Bool=true
      )

For a locally free sheaf ``E`` on an `AbsCoveredScheme` ``X`` this produces
the associated projectivization ``ℙ (E) → X`` as a `CoveredProjectiveScheme`.

A list of names for the variables of the relative homogeneous coordinate
rings can be provided with `var_names`.

!!! note: The sheaf ``E`` needs to be locally free so that a `trivializing_covering`
can be computed. The check for this can be turned off by setting `check=false`.
"""
function projectivization(E::AbsCoherentSheaf;
    var_names::Vector{<:VarName}=Vector{Symbol}(),
    check::Bool=true
  )
  X = scheme(E)
  check && (is_locally_free(E) || error("coherent sheaf must be locally free"))
  C = trivializing_covering(E)
  algebras = IdDict{AbsAffineScheme, Union{MPolyQuoRing, MPolyRing}}()
  on_patches = IdDict{AbsAffineScheme, AbsProjectiveScheme}()

  # Fill in the names of the variables in case there are none provided.
  if length(var_names) == 0
    # determine the global bound on the rank
    r = 0
    for U in patches(C)
      F = E(U)
      F isa FreeMod || error("modules must locally be free")
      r = (rank(F) > r ? rank(F) : r)
    end
    var_names = [Symbol(:s, i) for i in 0:r-1]
  end

  for U in patches(C)
    F = E(U)
    F isa FreeMod || error("modules must locally be free")
    r = rank(F)
    length(var_names) >= r || error("number of names for the variables must greater or equal to the local rank of the module")
    RU = rees_algebra(E(U); var_names=var_names[1:r])
    algebras[U] = RU
    SU, _ = grade(RU)
    PU = proj(SU)
    set_base_scheme!(PU, U)
    on_patches[U] = PU
  end
  projective_gluings = IdDict{Tuple{AbsAffineScheme, AbsAffineScheme}, AbsProjectiveGluing}()
  OX = StructureSheafOfRings(X)

  # prepare for the projective gluings
  for (U, V) in keys(gluings(C))
    P = on_patches[U]
    SP = homogeneous_coordinate_ring(P)
    Q = on_patches[V]
    SQ = homogeneous_coordinate_ring(Q)
    G = C[U, V]
    UV, VU = gluing_domains(G)
    f, g = gluing_morphisms(G)

    PUV, PUVtoP = fiber_product(OX(U, UV), P)
    QVU, QVUtoQ = fiber_product(OX(V, VU), Q)

    # to construct the identifications of PUV with QVU we need to
    # express the generators of I(U) in terms of the generators of I(V)
    # on the overlap U ∩ V.
    !(G isa Gluing) || error("method not implemented for this type of gluing")

    # The problem is that on a AffineSchemeOpenSubscheme U ∩ V
    # despite I(U)|U ∩ V == I(V)|U ∩ V, we
    # have no method to find coefficients aᵢⱼ such that fᵢ = ∑ⱼaᵢⱼ⋅gⱼ
    # for the generators fᵢ of I(U) and gⱼ of I(V): Even though
    # we can do this locally on the patches of a AffineSchemeOpenSubscheme, the result
    # is not guaranteed to glue to global functions on the overlap.
    # Abstractly, we know that the intersection of affine charts
    # in a separated scheme must be affine, but we do not have a
    # model of this overlap as an affine scheme and hence no computational
    # backup.

    # fᵢ the generators of I(U)
    # gⱼ the generators of I(V)
    # aᵢⱼ the coefficients for fᵢ = ∑ⱼ aᵢⱼ⋅gⱼ in VU
    # bⱼᵢ the coefficients for gⱼ = ∑ᵢ bⱼᵢ⋅fᵢ in UV
    # sᵢ the variables for the homogeneous ring over U
    # tⱼ the variables for the homogenesous ring over V
    A = [coordinates(E(U, VU)(v)) for v in gens(E(U))]
    B = [coordinates(E(V, UV)(w)) for w in gens(E(V))]
    #A = [coordinates(OX(U, VU)(f), I(VU)) for f in gens(I(U))] # A[i][j] = aᵢⱼ
    #B = [coordinates(OX(V, UV)(g), I(UV)) for g in gens(I(V))] # B[j][i] = bⱼᵢ
    SQVU = homogeneous_coordinate_ring(QVU)
    SPUV = homogeneous_coordinate_ring(PUV)
    # the induced map is ℙ(UV) → ℙ(VU), tⱼ ↦ ∑ᵢ bⱼᵢ ⋅ sᵢ
    # and ℙ(VU) → ℙ(UV), sᵢ ↦ ∑ⱼ aᵢⱼ ⋅ tⱼ
    fup = ProjectiveSchemeMor(PUV, QVU,
                              hom(SQVU, SPUV, pullback(f),
                                  [sum([B[j][i]*SPUV[i] for i in 1:ngens(SPUV)]) for j in 1:length(B)],
                                  check=false),
                              check=false
                            )
    gup = ProjectiveSchemeMor(QVU, PUV,
                              hom(SPUV, SQVU, pullback(g),
                                  [sum([A[i][j]*SQVU[j] for j in 1:ngens(SQVU)]) for i in 1:length(A)],
                                  check=false
                                 ),
                              check=false
                             )

    projective_gluings[U, V] = ProjectiveGluing(G, PUVtoP, QVUtoQ, fup, gup, check=false)
  end
  return CoveredProjectiveScheme(X, C, on_patches, projective_gluings, check=false)
end

function Base.hash(X::AbsCoherentSheaf, u::UInt)
  return u
end

