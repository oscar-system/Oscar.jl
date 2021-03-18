export AlgebraHomomorphism, codomain, compose, domain, hom,
       IdentityAlgebraHomomorphism, kernel, preimage
        

###############################################################################
#
#   IdentityAlgebraHomomorphism
#
###############################################################################

struct IdAlgHom{T} <: AbstractAlgebra.Map{Ring, Ring,
         AbstractAlgebra.IdentityMap, IdAlgHom} where T <: Union{AbstractAlgebra.Ring, AbstractAlgebra.Field}

   domain::Union{MPolyRing, MPolyQuo}
   image::Vector{U} where U <: Union{MPolyElem, MPolyQuoElem}
   salghom::Singular.SIdAlgHom
   kernel::Union{MPolyIdeal, MPolyQuoIdeal}

   function IdAlgHom{T}(R::U) where U <: Union{MPolyRing{T}, MPolyQuo{T}} where T
      V = gens(R)
      Sx = Oscar.singular_ring(R)
      ty = typeof(base_ring(Sx))
      z = new(R, V, Singular.IdentityAlgebraHomomorphism(Sx), ideal(R, [zero(R)]))
      return z
   end
end

function IdentityAlgebraHomomorphism(R::U) where U <: Union{MPolyRing{T}, MPolyQuo{T}} where T
   return IdAlgHom{T}(R)
end

###############################################################################
#
#   I/O for Identity Algebra Homomorphisms
#
###############################################################################

function show(io::IO, M::Map(IdAlgHom))
   println(io, "Identity algebra homomorphism with")
   println(io, "")
   println(io, "domain: ", domain(M))
   println(io, "")
   println(io, "defining equations: ", M.image)
end

###############################################################################
#
#   Basic Operations with Identity Algebra Homomorphisms
#
###############################################################################

function map_poly(f::Map(IdAlgHom), p::U) where U <: Union{MPolyElem, MPolyQuoElem}
   @assert parent(p) == domain(f)
   return p
end

function (f::IdAlgHom)(p::U) where U <: Union{MPolyElem, MPolyQuoElem}
   return map_poly(f, p)
end

###############################################################################
#
#   Preimage and Kernel for Identity Algebra Homomorphisms
#
###############################################################################

function preimage(f::Map(IdAlgHom), I::U) where U <: Union{MPolyIdeal, MPolyQuoIdeal}
   @assert base_ring(I) == domain(f)
   return I
end

function kernel(f::Map(IdAlgHom))
   return f.kernel
end

###############################################################################
#
#   AlgebraHomomorphism
#
###############################################################################

# this terrible function is used to get polynomials in and out of Singular and
# to change: polys over the quotient ring <=> polys over the non-quotient ring
function _badpolymap(f, R::MPolyRing)
  parent(f) == R && return f
  @assert ngens(parent(f)) == ngens(R)
  B = base_ring(R)
  g = MPolyBuildCtx(R)
  for (c, e) = zip(Nemo.coeffs(f), Nemo.exponent_vectors(f))
    push_term!(g, B(c), e)
  end
  return finish(g)
end

mutable struct AlgHom{T} <: AbstractAlgebra.Map{Ring, Ring,
         AbstractAlgebra.SetMap, AlgHom} where T <: Union{AbstractAlgebra.Ring, AbstractAlgebra.Field}
   domain::Union{MPolyRing, MPolyQuo}
   codomain::Union{MPolyRing, MPolyQuo}
   image::Vector{U} where U <: Union{MPolyElem, MPolyQuoElem}
   salghom::Singular.SAlgHom
   kernel::Union{MPolyIdeal, MPolyQuoIdeal}

   function AlgHom{T}(R::U, S::W, V::Vector{X}) where {T, U, W, X}
      Rx = singular_ring(R)
      Sx = singular_ring(S)
      if isdefined(S, :I) ## Check if S is a quotient ring
         Vx = map(p -> _badpolymap(S(p).f, Sx), V)
      else
         Vx = map(p -> _badpolymap(S(p), Sx), V)
      end

      z = new(R, S, V, Singular.AlgebraHomomorphism(Rx, Sx, Vx))
      return z
   end

   function AlgHom{T}(R::U, S::W, V::Vector{X}, phi::Singular.SAlgHom) where {T, U, W, X}
      return new(R, S, V, phi)
   end
end

###############################################################################
#
#   I/O for Algebra Homomorphisms
#
###############################################################################

function show(io::IO, M::Map(AlgHom))
   println(io, "Algebra homomorphism with")
   println(io, "")
   println(io, "domain: ", domain(M))
   println(io, "")
   println(io, "codomain: ", codomain(M))
   println(io, "")
   println(io, "defining equations: ", M.image)
end

###############################################################################
#
#   Algebra Homomorphism constructor
#
###############################################################################

#Compute the type of the base_ring of the underlying poly ring
function _type_helper(R)
   if isdefined(R, :I)
      return typeof(base_ring(R.R))
   else
      return typeof(base_ring(R))
   end
end

@doc Markdown.doc"""
    AlgebraHomomorphism(D::U, C::W, V::Vector{X}) where 
    {T, U <: Union{MPolyRing{T}, MPolyQuo}, 
    W <: Union{MPolyRing{T}, MPolyQuo}, 
    X <: Union{MPolyElem{T}, MPolyQuoElem}}

Creates the algebra homomorphism $D \rightarrow C$ defined by sending the $i$th generator of $D$ to the $i$th element of $V$. 
Allows types `MPolyRing` and `MPolyQuo` for $C$ and $D$ as well as entries of type `MPolyElem` and `MPolyQuoElem` for `X`.
Alternatively, use `hom(D::U, C::W, V::Vector{X})`.
"""
function AlgebraHomomorphism(D::U, C::W, V::Vector{X}) where 
    {T, U <: Union{MPolyRing{T}, MPolyQuo}, 
    W <: Union{MPolyRing{T}, MPolyQuo}, 
    X <: Union{MPolyElem{T}, MPolyQuoElem}}
   n = length(V)
   @assert n == ngens(D)
   ty = _type_helper(D)
   return AlgHom{ty}(D, C, V)
end

hom(D::U, C::W, V::Vector{X}) where {T,
   U <: Union{MPolyRing{T}, MPolyQuo}, W <: Union{MPolyRing{T}, MPolyQuo},
   X <: Union{MPolyElem{T}, MPolyQuoElem}} = AlgebraHomomorphism(D, C, V)

###############################################################################
#
#   Basic Operations with Algebra Homomorphisms
#
###############################################################################

function map_poly(F::Map(AlgHom), p::U) where U <: Union{MPolyElem, MPolyQuoElem}
   @assert parent(p) == domain(F)
   D = domain(F)
   Dx = domain(F.salghom)
   C = codomain(F)
   Cx = codomain(F.salghom)
   # TODO: _badpolymap really has to be replaced somehow ...
   if isdefined(D, :I) ## Check if D is a quotient ring
         px = _badpolymap(p.f, Dx)
   else
         px = _badpolymap(p, Dx)
   end

   if isdefined(C, :R) ## Check if C is a quotient ring
      return C(_badpolymap(F.salghom(px), C.R))
   else
      return _badpolymap(F.salghom(px), C)
   end
end

function (F::AlgHom)(p::U) where U <: Union{MPolyElem, MPolyQuoElem}
   return map_poly(F, p)
end

@doc Markdown.doc"""
    function domain(F::AlgHom)

Returns the domain of `F`.
"""
function domain(F::AlgHom)
   return F.domain
end

@doc Markdown.doc"""
    function codomain(F::AlgHom)

Returns the codomain of `F`.
"""
function codomain(F::AlgHom)
   return F.codomain
end

###############################################################################
#
#   Composition of Algebra Homomorphisms
#
###############################################################################

@doc Markdown.doc"""
    compose(F::AlgHom, G::AlgHom)

Returns the algebra homomorphism $H = G\circ F: domain(F) --> codomain(G)$.
"""
function compose(F::AlgHom, G::AlgHom)
   check_composable(F, G)
   phi = Singular.compose(F.salghom, G.salghom)
   D = domain(F)
   ty = _type_helper(D)
   C = codomain(G)
   V = C.(phi.image)
   return AlgHom{ty}(D, C, V, phi)
end

###############################################################################
#
#   Preimage and Kernel for Algebra Homomorphisms
#
###############################################################################

@doc Markdown.doc"""
    preimage(F::AlgHom, I::U) where U <: Union{MPolyIdeal, MPolyQuoIdeal}

Returns the preimage of the ideal $I$ under the algebra homomorphism $F$.
"""
function preimage(F::AlgHom, I::U) where U <: Union{MPolyIdeal, MPolyQuoIdeal}

   @assert base_ring(I) == codomain(F)
   D = domain(F)
   C = codomain(F)
   Cx = codomain(F.salghom)
   V = gens(I)
   if isdefined(C, :I) ## Check if C is a quotient ring
         Vx = map(p -> _badpolymap(C(p).f, Cx), V)
      else
         Vx = map(p -> _badpolymap(C(p), Cx), V)
      end
   Ix = Singular.Ideal(Cx, Vx)
   prIx = Singular.preimage(F.salghom, Ix)
   return ideal(D, D.(gens(prIx)))
end

@doc Markdown.doc"""
    kernel(F::AlgHom)

Returns the kernel of the algebra homomorphism $F$.
"""
function kernel(F::AlgHom)
   isdefined(F, :kernel) && return F.kernel
   C = codomain(F)
   F.kernel = preimage(F, ideal(C, [zero(C)]))
   return F.kernel
end