
########################################################################
# Dummy constructor for documentation only                             #
########################################################################
@doc raw"""
    Gluing(X::AbsAffineScheme, Y::AbsAffineScheme, f::SchemeMor, g::SchemeMor)

Glue two affine schemes ``X`` and ``Y`` along mutual isomorphisms
``f`` and ``g`` of open subsets ``U`` of ``X`` and ``V`` of ``Y``.

# Examples
```jldoctest
julia> P1, (x,y) = QQ[:x, :y]; P2, (u,v) = QQ[:u, :v];

julia> U1 = spec(P1); U2 = spec(P2);

julia> V1 = PrincipalOpenSubset(U1, x); # Preparations for gluing

julia> V2 = PrincipalOpenSubset(U2, u);

julia> f = morphism(V1, V2, [1//x, y//x]); # The gluing isomorphism

julia> g = morphism(V2, V1, [1//u, v//u]); # and its inverse

julia> G = Gluing(U1, U2, f, g) # Construct the gluing
Gluing
  of affine 2-space
  and affine 2-space
along the open subsets
  [x, y]   AA^2 \ scheme(x)
  [u, v]   AA^2 \ scheme(u)
given by the pullback function
  u -> 1/x
  v -> y/x

julia> G isa SimpleGluing # Since the gluing domains were `PrincipalOpenSubsets`, this defaults to a `SimpleGluing`
true

julia> # Alternative using AffineSchemeOpenSubschemes as gluing domains:

julia> W1 = AffineSchemeOpenSubscheme(U1, [x]); W2 = AffineSchemeOpenSubscheme(U2, [u]);

julia> h1 = AffineSchemeOpenSubschemeMor(W1, W2, [1//x, y//x]);

julia> h2 = AffineSchemeOpenSubschemeMor(W2, W1, [1//u, v//u]);

julia> H = Gluing(U1, U2, h1, h2)
Gluing
  of affine 2-space
  and affine 2-space
along the open subsets
  [x, y]   complement to V(x) in affine scheme with coordinates [x, y]
  [u, v]   complement to V(u) in affine scheme with coordinates [u, v]
defined by the map
  affine scheme morphism
    from [x, y]  AA^2 \ scheme(x)
    to   [u, v]  affine 2-space
  given by the pullback function
    u -> 1/x
    v -> y/x

julia> H isa Gluing
true
```
"""
function Gluing(X::AbsAffineScheme, Y::AbsAffineScheme, f::SchemeMor, g::SchemeMor)
  error("constructor for `Gluing` not implemented for morphisms of type $(typeof(f)) and $(typeof(g))")
end

########################################################################
# Conversion from SimpleGluings to Gluings                           #
########################################################################
function Gluing(G::SimpleGluing)
  X, Y = patches(G)
  U, V = gluing_domains(G)
  f, g = gluing_morphisms(G)
  incY = inclusion_morphism(V, Y)
  incX = inclusion_morphism(U, X)
  Uo = AffineSchemeOpenSubscheme(U)
  Vo = AffineSchemeOpenSubscheme(V)
  fo = AffineSchemeOpenSubschemeMor(Uo, Vo, [compose(f, incY)], check=false)
  go = AffineSchemeOpenSubschemeMor(Vo, Uo, [compose(g, incX)], check=false)
  return Gluing(X, Y, fo, go, check=false)
end

### For compatibility, return the gluing itself whenever it is of the right type
function Gluing(G::Gluing)
  return G
end

########################################################################
# Simplified constructor for common special case                       #
########################################################################

function Gluing(
    X::AbsAffineScheme, Y::AbsAffineScheme,
    f::AbsAffineSchemeMor{<:PrincipalOpenSubset},
    g::AbsAffineSchemeMor{<:PrincipalOpenSubset};
    check::Bool=true)
  return SimpleGluing(X, Y, f, g, check=check)
end

########################################################################
# Restrictions of Gluings to closed subschemes                        #
########################################################################
function restrict(G::AbsGluing, X::AbsAffineScheme, Y::AbsAffineScheme; check::Bool=true)
  return restrict(underlying_gluing(G), X, Y, check=check)
end

function restrict(G::AbsGluing, X::AbsAffineScheme, Y::AbsAffineScheme,
    Ures::Scheme, Vres::Scheme;
    check::Bool=true
  )
  return restrict(underlying_gluing(G), X, Y, Ures, Vres, check=check)
end


function restrict(G::Gluing, X::AbsAffineScheme, Y::AbsAffineScheme; check::Bool=true)
  U, V = gluing_domains(G)
  f, g = gluing_morphisms(G)
  @check is_closed_embedding(intersect(X, ambient_scheme(U)), ambient_scheme(U)) "the scheme is not a closed in the ambient scheme of the open set"
  @check is_closed_embedding(intersect(Y, ambient_scheme(V)), ambient_scheme(V)) "the scheme is not a closed in the ambient scheme of the open set"
  Ures = intersect(X, U, check=false)
  Vres = intersect(Y, V, check=false)
  return Gluing(X, Y, restrict(f, Ures, Vres, check=check), restrict(g, Vres, Ures, check=check), check=check)
end

function restrict(G::SimpleGluing, X::AbsAffineScheme, Y::AbsAffineScheme; check::Bool=true)
  U, V = gluing_domains(G)
  f, g = gluing_morphisms(G)
  @check is_closed_embedding(intersect(X, ambient_scheme(U)), ambient_scheme(U)) "the scheme is not a closed in the ambient scheme of the open set"
  @check is_closed_embedding(intersect(Y, ambient_scheme(V)), ambient_scheme(V)) "the scheme is not a closed in the ambient scheme of the open set"
  UX = PrincipalOpenSubset(X, OO(X)(lifted_numerator(complement_equation(U))))
  VY = PrincipalOpenSubset(Y, OO(Y)(lifted_numerator(complement_equation(V))))
  f_res = restrict(f, UX, VY, check=check)
  g_res = restrict(g, VY, UX, check=check)
  return SimpleGluing(X, Y, f_res, g_res, check=check)
end

function restrict(G::Gluing, X::AbsAffineScheme, Y::AbsAffineScheme,
    Ures::AffineSchemeOpenSubscheme, Vres::AffineSchemeOpenSubscheme;
    check::Bool=true
  )
  U, V = gluing_domains(G)
  f, g = gluing_morphisms(G)
  @check is_closed_embedding(intersect(X, ambient_scheme(U)), ambient_scheme(U)) "the scheme is not a closed in the ambient scheme of the open set"
  @check is_closed_embedding(intersect(Y, ambient_scheme(V)), ambient_scheme(V)) "the scheme is not a closed in the ambient scheme of the open set"
  return Gluing(X, Y, restrict(f, Ures, Vres, check=check), restrict(g, Vres, Ures, check=check), check=check)
end

function restrict(G::SimpleGluing, X::AbsAffineScheme, Y::AbsAffineScheme,
    UX::PrincipalOpenSubset, VY::PrincipalOpenSubset;
    check::Bool=true
  )
  U, V = gluing_domains(G)
  f, g = gluing_morphisms(G)
  @check is_closed_embedding(intersect(X, ambient_scheme(U)), ambient_scheme(U)) "the scheme is not a closed in the ambient scheme of the open set"
  @check is_closed_embedding(intersect(Y, ambient_scheme(V)), ambient_scheme(V)) "the scheme is not a closed in the ambient scheme of the open set"
  f_res = restrict(f, UX, VY, check=check)
  g_res = restrict(g, VY, UX, check=check)
  return SimpleGluing(X, Y, f_res, g_res, check=check)
end

########################################################################
# Identification of gluings under isomorphisms of patches             #
########################################################################
@doc raw"""
    restrict(G::AbsGluing, f::AbsAffineSchemeMor, g::AbsAffineSchemeMor; check::Bool=true)

Given a gluing ``X ↩ U ≅ V ↪ Y`` and isomorphisms ``f : X → X'`` and
``g: Y → Y'``, return the induced gluing of ``X'`` and ``Y'``.
"""
function restrict(G::AbsGluing, f::AbsAffineSchemeMor, g::AbsAffineSchemeMor; check::Bool=true)
  (X1, Y1) = patches(G)
  X1 === domain(f) || error("maps not compatible")
  X2 = codomain(f)
  finv = inverse(f)

  Y1 === domain(g) || error("maps not compatible")
  Y2 = codomain(g)
  ginv = inverse(g)

  (h1, h2) = gluing_morphisms(G)

  U2 = preimage(finv, domain(h1), check=check)
  V2 = preimage(ginv, domain(h2), check=check)

  return Gluing(X2, Y2,
                 compose(restrict(finv, U2, domain(h1), check=check),
                         compose(h1, restrict(g, domain(h2), V2, check=check))
                        ),
                 compose(restrict(ginv, V2, domain(h2), check=check),
                         compose(h2, restrict(f, domain(h1), U2, check=check))
                        ),
                 check=check
                )
end

function restrict(G::AbsGluing, f::AbsAffineSchemeMor, g::AbsAffineSchemeMor,
    f_res::SchemeMor, g_res::SchemeMor;
    check::Bool=true
  )
  (X1, Y1) = patches(G)
  X1 === domain(f) || error("maps not compatible")
  X2 = codomain(f)
  finv = inverse(f)

  Y1 === domain(g) || error("maps not compatible")
  Y2 = codomain(g)
  ginv = inverse(g)

  (h1, h2) = gluing_morphisms(G)

  U2 = codomain(f_res)
  V2 = codomain(f_res)

  return Gluing(X2, Y2,
                 compose(inverse(f_res),
                         compose(h1, g_res)
                        ),
                 compose(inverse(g_res),
                         compose(h2, f_res)
                        ),
                 check=check
                )
end

########################################################################
# Maximal extensions from SimpleGluings                               #
########################################################################

function maximal_extension(G::SimpleGluing)
  # We first check, whether this gluing already is maximal.
  # If it is not, we pass it on to the more complicated methods.
  U, V = gluing_domains(G)
  f, g = gluing_morphisms(G)
  X, Y = patches(G)
  x = gens(OO(X))
  y = gens(OO(Y))
  pbx = pullback(g).(x)
  pby = pullback(f).(y)
  U2, ext_y = maximal_extension(X, fraction.(pby))
  V2, ext_x = maximal_extension(Y, fraction.(pbx))
  is_subscheme(U2, U) && is_subscheme(V2, V) && return G
  return maximal_extension(Gluing(G))
end

