
function ==(f::AbsProjectiveSchemeMorphism, g::AbsProjectiveSchemeMorphism) 
  domain(f) === domain(g) || return false
  codomain(f) === codomain(g) || return false
  for s in gens(homogeneous_coordinate_ring(codomain(f)))
    iszero(pullback(f)(s) - pullback(g)(s)) || return false
  end
  return true
end

function ==(f::AbsProjectiveSchemeMorphism{<:AbsProjectiveScheme{<:Union{<:MPolyRing, <:MPolyQuoRing, <:MPolyLocRing, <:MPolyQuoLocRing}}},
            g::AbsProjectiveSchemeMorphism{<:AbsProjectiveScheme{<:Union{<:MPolyRing, <:MPolyQuoRing, <:MPolyLocRing, <:MPolyQuoLocRing}}})
  domain(f) === domain(g) || return false
  codomain(f) === codomain(g) || return false
  return map_on_affine_cones(f) == map_on_affine_cones(g)
end

@attr Any function covered_scheme_morphism(
    f::AbsProjectiveSchemeMorphism{<:Any, <:Any, <:Any, Nothing} # No map on base rings
  )
  PX = domain(f)
  PY = codomain(f)
  SX = ambient_coordinate_ring(PX)
  SY = ambient_coordinate_ring(PY)
  pbf = pullback(f) # The pullback on the free polynomial rings, not the quotients

  X = covered_scheme(PX)
  Y = covered_scheme(PY)

  mor_dict = IdDict{AbsAffineScheme, AbsAffineSchemeMor}()
  U = affine_charts(X)

  if ngens(SY) == ngens(SX) && all(k->pbf(SY[k]) == SX[k], 1:ngens(SY))
    for i in 1:ngens(SX)
      U_i = U[i]
      dehom = dehomogenization_map(PX, U_i) # the dehomogenization map SX → 𝒪(Uᵢ)
      y = gen(SY, i)
      denom = dehom(pbf(y))
      V_i = affine_charts(Y)[i]
      mor_dict[U_i] = AffineSchemeMor(U_i, V_i, 
                               hom(OO(V_i), OO(U_i), 
                                   [OO(U_i)(dehom(pbf(gen(SY, k)))) for k in 1:ngens(SY) if k != i];
                                   check=false
                                  )
                              )
    end
    phi = CoveringMorphism(default_covering(X), default_covering(Y), mor_dict, check=false)
    ff = CoveredSchemeMorphism(X, Y, phi; check=false)
    return ff
  end

  # the default case
  for i in 1:ngens(SX)
    U_i = U[i]
    dehom = dehomogenization_map(PX, U_i) # the dehomogenization map SX → 𝒪(Uᵢ)
    for j in 1:ngens(SY)
      y = gen(SY, j)
      denom = dehom(pbf(y))
      V_j = affine_charts(Y)[j]
      U_ij = PrincipalOpenSubset(U_i, denom)
      u = inv(OO(U_ij)(denom))
      mor_dict[U_ij] = morphism(U_ij, V_j, 
                               hom(OO(V_j), OO(U_ij), 
                                   [OO(U_ij)(dehom(pbf(gen(SY, k))))*u for k in 1:ngens(SY) if k != j];
                                   check=false
                                  )
                              )
    end
  end
  # We skip the gluings for the time being.
  # Eventually, they should be made lazy.
  CC = Covering(collect(keys(mor_dict)), IdDict{Tuple{AbsAffineScheme, AbsAffineScheme}, AbsGluing}())
  inherit_gluings!(CC, default_covering(X))
  phi = CoveringMorphism(CC, default_covering(Y), mor_dict, check=false)
  push!(coverings(X), CC)

  ff = CoveredSchemeMorphism(X, Y, phi; check=false)
  return ff
end

@doc raw"""
    covered_scheme_morphism(f::AbsProjectiveSchemeMorphism)

Given a morphism of `ProjectiveScheme`s ``f : X → Y``, construct and 
return the same morphism as a `CoveredSchemeMorphism` of the `covered_scheme`s 
of ``X`` and ``Y``, respectively.

# Examples
```jldoctest
julia> P, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z]);

julia> I = ideal([x^3-y^2*z]);

julia> Y = proj(P, I);

julia> f = id_hom(Y)
Projective scheme morphism
  from projective scheme in IP^2 over QQ
  to projective scheme in IP^2 over QQ

julia> fcov = covered_scheme_morphism(f);

julia> codomain(fcov)
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
```
"""
@attr Any function covered_scheme_morphism(
    f::AbsProjectiveSchemeMorphism
  ) # with map on the base rings
  PX = domain(f)
  PY = codomain(f)
  SX = ambient_coordinate_ring(PX)
  SY = ambient_coordinate_ring(PY)
  pbf = pullback(f) # The pullback on the free polynomial rings, not the quotients
  coeff_map = coefficient_map(pbf)

  X = covered_scheme(PX)
  Y = covered_scheme(PY)

  mor_dict = IdDict{AbsAffineScheme, AbsAffineSchemeMor}()
  U = affine_charts(X)
  if ngens(SY) == ngens(SX) && all(k->pullback(f)(SY[k]) == SX[k], 1:ngens(SY))
    for i in 1:ngens(SX)
      U_i = U[i]
      dehom = dehomogenization_map(PX, U_i) # the dehomogenization map SX → 𝒪(Uᵢ)
      y = gen(SY, i)
      denom = dehom(pbf(y))
      V_i = affine_charts(Y)[i]
      mor_dict[U_i] = AffineSchemeMor(U_i, V_i, 
                               hom(OO(V_i), OO(U_i), coeff_map,
                                   [OO(U_i)(dehom(pbf(gen(SY, k)))) for k in 1:ngens(SY) if k != i];
                                   check=false
                                  );
                               check=false
                              )
    end
    phi = CoveringMorphism(default_covering(X), default_covering(Y), mor_dict, check=false)
    ff = CoveredSchemeMorphism(X, Y, phi; check=false)
    return ff
  end

  # the default case
  for i in 1:ngens(SX)
    U_i = U[i]
    dehom = dehomogenization_map(PX, U_i) # the dehomogenization map SX → 𝒪(Uᵢ)
    for j in 1:ngens(SY)
      y = gen(SY, j)
      denom = dehom(pbf(y))
      V_j = affine_charts(Y)[j]
      U_ij = PrincipalOpenSubset(U_i, denom)
      u = inv(OO(U_ij)(denom))
      mor_dict[U_ij] = AffineSchemeMor(U_ij, V_j, 
                               hom(OO(V_j), OO(U_ij), coeff_map,
                                   [OO(U_ij)(dehom(pbf(gen(SY, k))))*u for k in 1:ngens(SY) if k != j];
                                   check=false
                                  );
                               check=false
                              )
      #@assert _has_coefficient_map(pullback(mor_dict[U_ij]))
    end
  end
  # We skip the gluings for the time being.
  # Eventually, they should be made lazy.
  CC = Covering(collect(keys(mor_dict)), IdDict{Tuple{AbsAffineScheme, AbsAffineSchemeMor}, AbsGluing}())
  inherit_gluings!(CC, default_covering(X))
  phi = CoveringMorphism(CC, default_covering(Y), mor_dict, check=false)
  push!(coverings(X), CC)

  ff = CoveredSchemeMorphism(X, Y, phi; check=false)
  return ff
end

@attr CoveredClosedEmbedding function covered_scheme_morphism(f::ProjectiveClosedEmbedding)
  PX = domain(f)
  PY = codomain(f)
  SX = homogeneous_coordinate_ring(PX)
  SY = homogeneous_coordinate_ring(PY)
  pbf = pullback(f)

  X = covered_scheme(PX)
  Y = covered_scheme(PY)

  mor_dict = IdDict{AbsAffineScheme, ClosedEmbedding}()
  U = affine_charts(X)
  # TODO: The code below does not run when we have empty affine charts. Adjust accordingly when cleaning up!
  II = IdealSheaf(PY, image_ideal(f))
  for i in 1:ngens(SX)
    U_i = U[i]
    V_i = affine_charts(Y)[i]
    mor_dict[U_i] = ClosedEmbedding(morphism(U_i, V_i, gens(OO(U_i)), check=false), II(V_i))
  end
  f_cov = CoveringMorphism(default_covering(X), default_covering(Y), mor_dict, check=false)

  result = CoveredClosedEmbedding(X, Y, f_cov, check=false)
  return result
end

function pushforward(inc::ProjectiveClosedEmbedding, M::FreeMod)
  f = pullback(inc)
  S = codomain(f)
  T = domain(f)
  S === base_ring(M) || error("rings do not match")
  FT = graded_free_module(T, [_degree_fast(a) for a in gens(M)])
  I = image_ideal(inc)
  IFT, inc_IFT = I*FT
  MT = cokernel(inc_IFT)
  id = hom(MT, M, gens(M), f; check=false)
  return MT, id
end

function pushforward(inc::ProjectiveClosedEmbedding, M::SubquoModule)
  F = ambient_free_module(M)
  f = pullback(inc)
  S = codomain(f)
  T = domain(f)
  FT, id_F = pushforward(inc, F)
  g = ambient_representatives_generators(M)
  gT = elem_type(FT)[sum(lift(x[i])*FT[i] for i in 1:length(g); init=zero(FT)) for x in coordinates.(g)]
  rel = relations(M)
  relT = elem_type(FT)[sum(lift(x[i])*FT[i] for i in 1:length(g); init=zero(FT)) for x in coordinates.(rel)]
  G, inc_G = sub(FT, vcat(gT, relT))
  Q, inc_Q = sub(G, gens(G)[length(gT)+1:end])
  MT = cokernel(inc_Q)
  id = hom(MT, M, vcat(gens(M), elem_type(M)[zero(M) for i in 1:length(relT)]), S; check=false)
  return MT, id
end

function SubquoModule(F::FreeMod, g::Vector{T}, rel::Vector{T}) where {T<:FreeModElem}
  G = SubModuleOfFreeModule(F, g)
  Rel = SubModuleOfFreeModule(F, rel)
  return SubquoModule(F, G, Rel)
end

###############################################################################
#
#  Printing
#
###############################################################################

function Base.show(io::IO, f::AbsProjectiveSchemeMorphism)
  if get(io, :show_semi_compact, false)
    _show_semi_compact(io, f)
  elseif is_terse(io)
    print(io, "Projective scheme morphism")
  else
    io = pretty(io)
    print(io, "Hom: ", Lowercase(), domain(f), " -> ", Lowercase(), codomain(f))
  end
end

function Base.show(io::IO, ::MIME"text/plain", f::AbsProjectiveSchemeMorphism)
  io = pretty(io)
  println(io, "Projective scheme morphism")
  print(io, Indent())
  println(io, "from ", Lowercase(), domain(f))
  print(io, "to ", Lowercase(), codomain(f))
  print(io, Dedent())
  if has_attribute(f, :covered_scheme_morphism)
    println(io)
    println(io, "defined by the map")
    print(io, Indent(), Lowercase())
    show(io, MIME"text/plain"(), covered_scheme_morphism(f))
    print(io, Dedent())
  end
end

function _show_semi_compact(io::IO, f::AbsProjectiveSchemeMorphism)
  io = pretty(io)
  print(io, "Projective scheme morphism")
  if has_attribute(f, :covered_scheme_morphism)
    println(io, " defined by the map")
    print(io, Indent(), Lowercase())
    show(io, MIME"text/plain"(), covered_scheme_morphism(f))
    print(io, Dedent())
  end
end

