export compose, maximal_extension, restrict

########################################################################
# Dummy constructor for documentation only                             #
########################################################################
@Markdown.doc """
    Glueing(X::AbsSpec, Y::AbsSpec, f::SchemeMor, g::SchemeMor)

Glue two affine schemes ``X`` and ``Y`` along mutual isomorphisms 
``f`` and ``g`` of open subsets ``U`` of ``X`` and ``V`` of ``Y``.

# Examples
```jldoctest
julia> P1, (x,y) = QQ["x", "y"]; P2, (u,v) = QQ["u", "v"];

julia> U1 = Spec(P1); U2 = Spec(P2);

julia> V1 = PrincipalOpenSubset(U1, x); # Preparations for glueing

julia> V2 = PrincipalOpenSubset(U2, u);

julia> f = SpecMor(V1, V2, [1//x, y//x]); # The glueing isomorphism

julia> g = SpecMor(V2, V1, [1//u, v//u]); # and its inverse

julia> G = Glueing(U1, U2, f, g) # Construct the glueing
Glueing of Spec of Multivariate Polynomial Ring in x, y over Rational Field and Spec of Multivariate Polynomial Ring in u, v over Rational Field along the map morphism from

	Spec of localization of Multivariate Polynomial Ring in x, y over Rational Field at the powers of fmpq_mpoly[x]

to

	Spec of localization of Multivariate Polynomial Ring in u, v over Rational Field at the powers of fmpq_mpoly[u]

with coordinates

	1//x, y//x

julia> typeof(G)<:SimpleGlueing # Since the glueing domains were `PrincipalOpenSubsets`, this defaults to a `SimpleGlueing`
true

julia> # Alternative using SpecOpens as glueing domains:

julia> W1 = SpecOpen(U1, [x]); W2 = SpecOpen(U2, [u]);

julia> h1 = SpecOpenMor(W1, W2, [1//x, y//x]);

julia> h2 = SpecOpenMor(W2, W1, [1//u, v//u]);

julia> H = Glueing(U1, U2, h1, h2)
Glueing of Spec of Multivariate Polynomial Ring in x, y over Rational Field and Spec of Multivariate Polynomial Ring in u, v over Rational Field along the map Morphism from complement of zero locus of fmpq_mpoly[x] in Spec of Multivariate Polynomial Ring in x, y over Rational Field to complement of zero locus of fmpq_mpoly[u] in Spec of Multivariate Polynomial Ring in u, v over Rational Field

julia> typeof(H)<:Glueing
true
```
"""
function Glueing(X::AbsSpec, Y::AbsSpec, f::SchemeMor, g::SchemeMor)
  error("constructor for `Glueing` not implemented for morphisms of type $(typeof(f)) and $(typeof(g))")
end

########################################################################
# Conversion from SimpleGlueings to Glueings                           #
########################################################################
function Glueing(G::SimpleGlueing)
  X, Y = patches(G)
  U, V = glueing_domains(G)
  f, g = glueing_morphisms(G)
  incY = inclusion_morphism(V, Y)
  incX = inclusion_morphism(U, X)
  Uo = SpecOpen(U)
  Vo = SpecOpen(V)
  fo = SpecOpenMor(Uo, Vo, [compose(f, incY)], check=false)
  go = SpecOpenMor(Vo, Uo, [compose(g, incX)], check=false)
  return Glueing(X, Y, fo, go, check=false)
end

### For compatibility, return the glueing itself whenever it is of the right type
function Glueing(G::Glueing)
  return G
end

########################################################################
# Simplified constructor for common special case                       #
########################################################################

function Glueing(
    X::AbsSpec, Y::AbsSpec, 
    f::AbsSpecMor{<:PrincipalOpenSubset}, 
    g::AbsSpecMor{<:PrincipalOpenSubset};
    check::Bool=true)
  return SimpleGlueing(X, Y, f, g, check=check)
end

########################################################################
# Restrictions of Glueings to closed subschemes                        #
########################################################################
function restrict(G::Glueing, X::AbsSpec, Y::AbsSpec; check::Bool=true)
  U, V = glueing_domains(G)
  f, g = glueing_morphisms(G)
  if check
    is_closed_embedding(intersect(X, ambient_scheme(U)), ambient_scheme(U)) || error("the scheme is not a closed in the ambient scheme of the open set")
    is_closed_embedding(intersect(Y, ambient_scheme(V)), ambient_scheme(V)) || error("the scheme is not a closed in the ambient scheme of the open set")
  end
  Ures = intersect(X, U)
  Vres = intersect(Y, V)
  return Glueing(X, Y, restrict(f, Ures, Vres, check=check), restrict(g, Vres, Ures, check=check), check=check)
end

function restrict(G::SimpleGlueing, X::AbsSpec, Y::AbsSpec; check::Bool=true)
  U, V = glueing_domains(G)
  f, g = glueing_morphisms(G)
  if check
    is_closed_embedding(intersect(X, ambient_scheme(U)), ambient_scheme(U)) || error("the scheme is not a closed in the ambient scheme of the open set")
    is_closed_embedding(intersect(Y, ambient_scheme(V)), ambient_scheme(V)) || error("the scheme is not a closed in the ambient scheme of the open set")
  end
  UX = PrincipalOpenSubset(X, OO(X)(lifted_numerator(complement_equation(U))))
  VY = PrincipalOpenSubset(Y, OO(Y)(lifted_numerator(complement_equation(V))))
  f_res = restrict(f, UX, VY, check=check)
  g_res = restrict(g, VY, UX, check=check)
  return SimpleGlueing(X, Y, f_res, g_res)
end

########################################################################
# Identification of glueings under isomorphisms of patches             #
########################################################################
@Markdown.doc """
    restrict(G::AbsGlueing, f::AbsSpecMor, g::AbsSpecMor; check::Bool=true)

Given a glueing ``X ↩ U ≅ V ↪ Y`` and isomorphisms ``f : X → X'`` and 
``g: Y → Y'``, return the induced glueing of ``X'`` and ``Y'``.
"""
function restrict(G::AbsGlueing, f::AbsSpecMor, g::AbsSpecMor; check::Bool=true)
  (X1, Y1) = patches(G)
  X1 === domain(f) || error("maps not compatible")
  X2 = codomain(f)
  finv = inverse(f)

  Y1 === domain(g) || error("maps not compatible")
  Y2 = codomain(g)
  ginv = inverse(g)

  (h1, h2) = glueing_morphisms(G)

  U2 = preimage(finv, domain(h1), check=check)
  V2 = preimage(ginv, domain(h2), check=check)

  return Glueing(X2, Y2, 
                 compose(restrict(finv, U2, domain(h1), check=check), 
                         compose(h1, restrict(g, domain(h2), V2, check=check), check=check),
                         check=check
                        ),
                 compose(restrict(ginv, V2, domain(h2), check=check), 
                         compose(h2, restrict(f, domain(h1), U2), check=check), 
                         check=check
                        ),
                 check=check
                )
end
