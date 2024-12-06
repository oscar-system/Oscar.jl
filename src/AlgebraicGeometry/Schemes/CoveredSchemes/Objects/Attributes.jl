
########################################################################
# Attributes of AbsCoveredScheme                                       #
########################################################################

########################################################################
# Type getters                                                         #
########################################################################
base_ring_type(::Type{T}) where {BRT, T<:AbsCoveredScheme{BRT}} = BRT

########################################################################
# Basic getters                                                        #
########################################################################
base_ring(X::AbsCoveredScheme) = base_ring(underlying_scheme(X))

@doc raw"""
    coverings(X::AbsCoveredScheme)

Return the list of internally stored `Covering`s of ``X``.

# Examples
```jldoctest
julia> P = projective_space(QQ, 2);

julia> Pcov = covered_scheme(P)
Scheme
  over rational field
with default covering
  described by patches
    1: affine 2-space
    2: affine 2-space
    3: affine 2-space
  in the coordinate(s)
    1: [(s1//s0), (s2//s0)]
    2: [(s0//s1), (s2//s1)]
    3: [(s0//s2), (s1//s2)]

julia> coverings(Pcov)
1-element Vector{Covering{QQField}}:
 Covering with 3 patches
```
"""
function coverings(X::AbsCoveredScheme) ::Vector{<:Covering}
  return coverings(underlying_scheme(X))
end

@doc raw"""
    default_covering(X::AbsCoveredScheme)

Return the default covering for ``X``.

# Examples
```jldoctest
julia> P = projective_space(QQ, 2);

julia> S = homogeneous_coordinate_ring(P);

julia> I = ideal(S, [S[1]*S[2]-S[3]^2]);

julia> X = subscheme(P, I)
Projective scheme
  over rational field
defined by ideal (s0*s1 - s2^2)

julia> Xcov = covered_scheme(X)
Scheme
  over rational field
with default covering
  described by patches
    1: scheme((s1//s0) - (s2//s0)^2)
    2: scheme((s0//s1) - (s2//s1)^2)
    3: scheme((s0//s2)*(s1//s2) - 1)
  in the coordinate(s)
    1: [(s1//s0), (s2//s0)]
    2: [(s0//s1), (s2//s1)]
    3: [(s0//s2), (s1//s2)]

julia> default_covering(Xcov)
Covering
  described by patches
    1: scheme((s1//s0) - (s2//s0)^2)
    2: scheme((s0//s1) - (s2//s1)^2)
    3: scheme((s0//s2)*(s1//s2) - 1)
  in the coordinate(s)
    1: [(s1//s0), (s2//s0)]
    2: [(s0//s1), (s2//s1)]
    3: [(s0//s2), (s1//s2)]

```
"""
function default_covering(X::AbsCoveredScheme)
  return default_covering(underlying_scheme(X))::Covering
end

@doc raw"""
    patches(X::AbsCoveredScheme) = patches(default_covering(X))

Return the affine patches in the `default_covering` of ``X``.
"""
patches(X::AbsCoveredScheme) = patches(default_covering(X))

@doc raw"""
    affine_charts(X::AbsCoveredScheme)

Return the affine charts in the `default_covering` of ``X``.

# Examples
```jldoctest
julia> P = projective_space(QQ, 2);

julia> S = homogeneous_coordinate_ring(P);

julia> I = ideal(S, [S[1]*S[2]-S[3]^2]);

julia> X = subscheme(P, I)
Projective scheme
  over rational field
defined by ideal (s0*s1 - s2^2)

julia> Xcov = covered_scheme(X)
Scheme
  over rational field
with default covering
  described by patches
    1: scheme((s1//s0) - (s2//s0)^2)
    2: scheme((s0//s1) - (s2//s1)^2)
    3: scheme((s0//s2)*(s1//s2) - 1)
  in the coordinate(s)
    1: [(s1//s0), (s2//s0)]
    2: [(s0//s1), (s2//s1)]
    3: [(s0//s2), (s1//s2)]

julia> affine_charts(Xcov)
3-element Vector{AffineScheme{QQField, MPolyQuoRing{QQMPolyRingElem}}}:
 scheme((s1//s0) - (s2//s0)^2)
 scheme((s0//s1) - (s2//s1)^2)
 scheme((s0//s2)*(s1//s2) - 1)

```
"""
affine_charts(X::AbsCoveredScheme) = basic_patches(default_covering(X))


########################################################################
# Attributes of CoveredScheme                                          #
########################################################################

########################################################################
# Basic getters                                                        #
########################################################################
coverings(X::CoveredScheme) = X.coverings
default_covering(X::CoveredScheme) = X.default_covering
patches(X::CoveredScheme) = patches(default_covering(X))
gluings(X::CoveredScheme) = gluings(default_covering(X))
base_ring(X::CoveredScheme) = X.kk


########################################################################
# Names of CoveredSchemes                                              #
########################################################################
set_name!(X::AbsCoveredScheme, name::String) = set_attribute!(X, :name, name)
name(X::AbsCoveredScheme) = get_attribute(X, :name)::String
has_name(X::AbsCoveredScheme) = has_attribute(X, :name)

########################################################################
# Auxiliary attributes                                                 #
########################################################################
function dim(X::AbsCoveredScheme)
  if !has_attribute(X, :dim)
    d = -1
    is_equidimensional=true
    for U in patches(default_covering(X))
      e = dim(U)
      if e > d
        d == -1 || (is_equidimensional=false)
        d = e
      end
    end
    set_attribute!(X, :dim, d)
    if !is_equidimensional
      # the above is not an honest check for equidimensionality,
      # because in each chart the output of `dim` is only the
      # supremum of all components. Thus we can only infer
      # non-equidimensionality in case this is already visible
      # from comparing the different charts
      set_attribute!(X, :is_equidimensional, false)
    end
  end
  return get_attribute(X, :dim)::Int
end

@attr Any function singular_locus_reduced(X::AbsCoveredScheme)
  D = IdDict{AbsAffineScheme, Ideal}()
  for U in affine_charts(X)
    _, inc_sing = singular_locus_reduced(U)
    D[U] = image_ideal(inc_sing)
  end
  Ising = IdealSheaf(X, D)
  inc = CoveredClosedEmbedding(X, Ising)
  return domain(inc), inc
end

@doc raw"""
    singular_locus(X::AbsCoveredScheme) -> AbsCoveredScheme, AbsCoveredSchemeMorphism

Return the singular locus of `X` as a covered scheme.

For the singular locus of the reduced scheme induced by `X`, please use
`singular_locus_reduced`.

# Examples
```jldoctest
julia> P, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z]);

julia> Y = variety(ideal([x^3-y^2*z]))
Projective variety
  in projective 2-space over QQ with coordinates [x, y, z]
defined by ideal (x^3 - y^2*z)

julia> Ycov = covered_scheme(Y)
Scheme
  over rational field
with default covering
  described by patches
    1: scheme(-(y//x)^2*(z//x) + 1)
    2: scheme((x//y)^3 - (z//y))
    3: scheme((x//z)^3 - (y//z)^2)
  in the coordinate(s)
    1: [(y//x), (z//x)]
    2: [(x//y), (z//y)]
    3: [(x//z), (y//z)]

julia> I, s = singular_locus(Ycov)
(Scheme over QQ covered with 1 patch, Hom: scheme over QQ covered with 1 patch -> scheme over QQ covered with 3 patches)

julia> I # singular locus actually lives in the patch {z != 0}
Scheme
  over rational field
with default covering
  described by patches
    1: scheme((x//z)^3 - (y//z)^2, (y//z), (x//z))
  in the coordinate(s)
    1: [(x//z), (y//z)]

julia> s
Covered scheme morphism
  from scheme over QQ covered with 1 patch
    1a: [(x//z), (y//z)]   scheme((x//z)^3 - (y//z)^2, (y//z), (x//z))
  to scheme over QQ covered with 3 patches
    1b: [(y//x), (z//x)]   scheme(-(y//x)^2*(z//x) + 1)
    2b: [(x//y), (z//y)]   scheme((x//y)^3 - (z//y))
    3b: [(x//z), (y//z)]   scheme((x//z)^3 - (y//z)^2)
given by the pullback function
  1a -> 3b
    (x//z) -> 0
    (y//z) -> 0
```
"""
@attr Any function singular_locus(
    X::AbsCoveredScheme;
  )
  D = IdDict{AbsAffineScheme, Ideal}()
  covering = get_attribute(X, :simplified_covering, default_covering(X))
  for U in covering
    _, inc_sing = singular_locus(U)
    D[U] = image_ideal(inc_sing)
  end
  Ising = IdealSheaf(X, D)
  inc = CoveredClosedEmbedding(X, Ising)
  return domain(inc), inc
end

function ideal_sheaf_of_singular_locus(
    X::AbsCoveredScheme;
    focus=zero_ideal_sheaf(X) # This should really be an AbsIdealSheaf, but the inclusion order forbids mentioning this here.
  )
  return get_attribute!(X, :ideal_sheaf_of_singular_locus) do
    SingularLocusIdealSheaf(X; focus)
  end::SingularLocusIdealSheaf
  D = IdDict{AbsAffineScheme, Ideal}()
  covering = get_attribute(X, :simplified_covering, default_covering(X))
  for U in covering
    _, inc_sing = singular_locus(U)
    D[U] = radical(image_ideal(inc_sing))
  end
  Ising = IdealSheaf(X, D, check=false)
  return Ising
end

function simplified_covering(X::AbsCoveredScheme)
  if !has_attribute(X, :simplified_covering)
    simplify!(X)
  end
  return get_attribute(X, :simplified_covering)::Covering
end
