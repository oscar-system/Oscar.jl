
function ==(f::ProjectiveSchemeMor, g::ProjectiveSchemeMor) 
  domain(f) === domain(g) || return false
  codomain(f) === codomain(g) || return false
  for s in gens(homogeneous_coordinate_ring(codomain(f)))
    iszero(pullback(f)(s) - pullback(g)(s)) || return false
  end
  return true
end

function ==(f::ProjectiveSchemeMor{<:AbsProjectiveScheme{<:Union{<:MPolyRing, <:MPolyQuoRing, <:MPolyLocRing, <:MPolyQuoLocRing}}},
            g::ProjectiveSchemeMor{<:AbsProjectiveScheme{<:Union{<:MPolyRing, <:MPolyQuoRing, <:MPolyLocRing, <:MPolyQuoLocRing}}})
  domain(f) === domain(g) || return false
  codomain(f) === codomain(g) || return false
  return map_on_affine_cones(f) == map_on_affine_cones(g)
end

@doc raw"""
    covered_scheme_morphism(f::ProjectiveSchemeMor)

Given a morphism of `ProjectiveScheme`s ``f : X → Y``, construct and 
return the same morphism as a `CoveredSchemeMorphism` of the `covered_scheme`s 
of ``X`` and ``Y``, respectively.

# Examples
```jldoctest
julia> P, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z]);

julia> I = ideal([x^3-y^2*z]);

julia> Y = projective_scheme(P, I);

julia> f = identity_map(Y)
Morphism
  from projective scheme in IP^2 over QQ
  to   projective scheme in IP^2 over QQ

julia> covered_scheme_morphism(f)
Morphism
  from scheme over QQ covered with 9 patches
    1a: [(y//x), (z//x)]   spec of localized quotient of multivariate polynomial ring
    2a: [(x//y), (z//y)]   spec of localized quotient of multivariate polynomial ring
    3a: [(y//x), (z//x)]   spec of localized quotient of multivariate polynomial ring
    4a: [(x//y), (z//y)]   spec of localized quotient of multivariate polynomial ring
    5a: [(x//y), (z//y)]   spec of localized quotient of multivariate polynomial ring
    6a: [(y//x), (z//x)]   spec of localized quotient of multivariate polynomial ring
    7a: [(x//z), (y//z)]   spec of localized quotient of multivariate polynomial ring
    8a: [(x//z), (y//z)]   spec of localized quotient of multivariate polynomial ring
    9a: [(x//z), (y//z)]   spec of localized quotient of multivariate polynomial ring
  to   scheme over QQ covered with 3 patches
    1b: [(y//x), (z//x)]   spec of quotient of multivariate polynomial ring
    2b: [(x//y), (z//y)]   spec of quotient of multivariate polynomial ring
    3b: [(x//z), (y//z)]   spec of quotient of multivariate polynomial ring
given by
  1a -> 2b
    (x//y) -> 1/(y//x)
    (z//y) -> (z//x)/(y//x)
    ----------------------------------------
  2a -> 2b
    (x//y) -> (x//y)
    (z//y) -> (z//y)
    ----------------------------------------
  3a -> 3b
    (x//z) -> 1/(z//x)
    (y//z) -> (y//x)/(z//x)
    ----------------------------------------
  4a -> 3b
    (x//z) -> (x//y)/(z//y)
    (y//z) -> 1/(z//y)
    ----------------------------------------
  5a -> 1b
    (y//x) -> 1/(x//y)
    (z//x) -> (z//y)/(x//y)
    ----------------------------------------
  6a -> 1b
    (y//x) -> (y//x)
    (z//x) -> (z//x)
    ----------------------------------------
  7a -> 1b
    (y//x) -> (y//z)/(x//z)
    (z//x) -> 1/(x//z)
    ----------------------------------------
  8a -> 3b
    (x//z) -> (x//z)
    (y//z) -> (y//z)
    ----------------------------------------
  9a -> 2b
    (x//y) -> (x//z)/(y//z)
    (z//y) -> 1/(y//z)
```
"""
@attr function covered_scheme_morphism(f::ProjectiveSchemeMor)
  PX = domain(f)
  PY = codomain(f)
  SX = ambient_coordinate_ring(PX)
  SY = ambient_coordinate_ring(PY)
  pbf = pullback(f) # The pullback on the free polynomial rings, not the quotients

  X = covered_scheme(PX)
  Y = covered_scheme(PY)

  mor_dict = IdDict{AbsSpec, AbsSpecMor}()
  U = affine_charts(X)
  for i in 1:ngens(SX)
    U_i = U[i]
    dehom = dehomogenization_map(PX, U_i) # the dehomogenization map SX → 𝒪(Uᵢ)
    for j in 1:ngens(SY)
      y = gen(SY, j)
      denom = dehom(pbf(y))
      V_j = affine_charts(Y)[j]
      U_ij = PrincipalOpenSubset(U_i, denom)
      u = inv(OO(U_ij)(denom))
      mor_dict[U_ij] = SpecMor(U_ij, V_j, 
                               hom(OO(V_j), OO(U_ij), 
                                   [OO(U_ij)(dehom(pbf(gen(SY, k))))*u for k in 1:ngens(SY) if k != j]
                                  )
                              )
    end
  end
  # We skip the glueings for the time being.
  # Eventually, they should be made lazy.
  CC = Covering(collect(keys(mor_dict)), IdDict{Tuple{AbsSpec, AbsSpec}, AbsGlueing}())
  inherit_glueings!(CC, default_covering(X))
  phi = CoveringMorphism(CC, default_covering(Y), mor_dict)
  push!(coverings(X), CC)

  ff = CoveredSchemeMorphism(X, Y, phi)
  return ff
end

###############################################################################
#
#  Printing
#
###############################################################################

function Base.show(io::IO, f::ProjectiveSchemeMor)
  if get(io, :supercompact, false)
    print(io, "Morphism")
  else
    io = pretty(io)
    print(io, "Morphism: ", Lowercase(), domain(f), " -> ", Lowercase(), codomain(f))
  end
end

function Base.show(io::IO, ::MIME"text/plain", f::ProjectiveSchemeMor)
  io = pretty(io)
  X = domain(f)
  Y = codomain(f)
  println(io, "Morphism")
  print(io, Indent(), "from ")
  if typeof(X) <: AbsProjectiveScheme{<:Field, <:MPolyAnyRing} # X is not a V(bla)
    print(io, Lowercase())
  end
  println(io, X)
  print(io, "to   ")
  if typeof(Y) <: AbsProjectiveScheme{<:Field, <:MPolyAnyRing} # same as above
    print(io, Lowercase())
  end
  print(io, Y)
  print(io, Dedent())
  if has_attribute(f, :covered_scheme_morphism)
    println(io)
    println(io, "defined by the map")
    print(io, Indent(), Lowercase())
    show(io, MIME"text/plain"(), covered_scheme_morphism(f))
    print(io, Dedent())
  end
end

function _show_semi_compact(io::IO, f::ProjectiveSchemeMor)
  io = pretty(io)
  print(io, "Morphism of projective schemes")
  if has_attribute(f, :covered_scheme_morphism)
    println(io, " given by")
    print(io, Indent(), Lowercase())
    show(io, MIME"text/plain"(), covered_scheme_morphism(f))
    print(io, Dedent())
  end
end

