
import AbstractAlgebra.Ring, Oscar.AlgHom, Oscar.compose
import Base: ∘

export AffineScheme, PrincipalSubScheme,Spec, SpecPrincipalOpen, affine_space
export AffSchMorphism, base_ring,ambient_ring,defining_ideal
export localize

abstract type Scheme end
abstract type AffineScheme end
abstract type SchemeMorphism end

mutable struct Spec{S <: Ring, T <: MPolyRing, U <: MPolyElem} <: AffineScheme
  # the basic fields 
  k::S			# the base ring (usually a field) of definition for the scheme
  R::T		  	# the ambient polynomial ring to model this affine scheme
  I::MPolyIdeal{U}	# The ideal in R defining the scheme

  function Spec(k::S, R::T, I::MPolyIdeal{U} ) where{
			S <: Ring, T <:MPolyRing , U <: MPolyElem}
    if k != coefficient_ring(R)
      error( "Base ring of the affine scheme does not coincide with the base ring of the associated algebra." )
    end
    # TODO: Implement further plausibility checks to be performed at runtime.
    return new{S, T, U}(k, R, I )
  end
end

function base_ring(A::Spec)
  return A.k
end

function ambient_ring(A::Spec)
  return A.R
end

function defining_ideal(A::Spec)
  return A.I
end

mutable struct SpecPrincipalOpen{S <: Ring, T<:MPolyRing, U<:MPolyElem} <: AffineScheme
  parent::Union{Spec{S,T,U},SpecPrincipalOpen{S,T,U}}
  denom::U  # element of the "ambient" polynomial ring of the root
  # these fields are only initialized if needed and used to comply to the AffineScheme interface
  k::S  
  R::T
  I::MPolyIdeal{U} 
  
  function SpecPrincipalOpen(parent::Union{Spec{S,T,U},SpecPrincipalOpen{S,T,U}}, denom::U) where {S <: Ring, T<:MPolyRing, U<:MPolyElem}
    x = new{S,T,U}() 
    x.parent = parent
    x.denom = denom
    return x 
  end 
end

function PrincipalSubScheme(parent, denom)
  return SpecPrincipalOpen(parent, denom)
end

function localize(parent, denom)
  return SpecPrincipalOpen(parent, denom)
end

function parent(S::AffineScheme)
  return S
end

function parent(D::SpecPrincipalOpen)
  return D.parent
end 

function root(D::SpecPrincipalOpen)
  return parent(parent(D))
end

function base_ring(D::SpecPrincipalOpen)
  if isdefined(D,Symbol("k"))
    return D.k
  end
  k = base_ring(root(D))
  D.k = k
  return k
end


function denoms(D::SpecPrincipalOpen{S,T,U}) where {S <: Ring, T <: MPolyRing, U<:MPolyElem}
  result = U[]
  P = D
  while typeof(P) <: SpecPrincipalOpen
    push!(result, P.denom)
    P = parent(P)
  end
  return result
end 

function ambient_ring(D::SpecPrincipalOpen)
  if isdefined(D,Symbol("R"))
    return D.R
  end
  names = ["t$i" for i in 1:length(denoms(D))]
  R = ambient_ring(root(D))
  names = vcat( names, String.( symbols( R)) )
  S, x = PolynomialRing( base_ring(D), names)
  D.R = S
  return S
end 

function get_ring_hom(D::SpecPrincipalOpen)
  R = ambient_ring(root(D))
  S = ambient_ring(D)
  d = length(denoms(D))
  n = length(gens(R))
  return AlgebraHomomorphism(R,S, gens(S)[d+1:n+d])
end

function defining_ideal(D::SpecPrincipalOpen)
  phi = get_ring_hom(D)
  I = defining_ideal(root(D))
  den = denoms(D)
  R = ambient_ring(root(D))
  S = ambient_ring(D)
  gen = gens(S)
  tmp = [gen[i]*phi(den[i]) - 1 for i in 1:length(den)]
  for g in gens(R)
    push!(tmp,phi(g))
  end
  return ideal(tmp)
end


# outer constructors
function Spec( k::S, R::T ) where{S <: Ring, T <:MPolyRing}
  I = ideal(R, zero(R))
  return Spec(k, R, I )
end

@doc Markdown.doc"""
    Spec(R::MPolyRing) -> Spec

Return the affine scheme corresponding to the ring $R/(0)$.
"""
function Spec( R::T ) where{T <: MPolyRing}
  I = ideal(R, zero(R))
  k = coefficient_ring( R )
  return Spec(k, R, I )
end

function Spec( R::T, I::MPolyIdeal{U} ) where{ T <: MPolyRing, U <: MPolyElem }
  k = coefficient_ring(R)
  return Spec(k, R, I )
end

# Construct affine n-space over the ring k.
function affine_space( k::Ring, n::Int, name::String="x" )
  R, x = PolynomialRing( k, name => (1:n))
  return Spec( R )
end

function Base.show( io::Base.IO, X::AffineScheme )
  Base.print( io, "Affine scheme over " )
  Base.print( io, base_ring(X) )
  Base.print( io, "\n" ) 
  Base.print( io, "given as the quotient of the polynomial ring \n" )
  Base.print( io, ambient_ring(X) )
  Base.print( io, "\nby the ideal \n" )
  Base.print( io, defining_ideal(X) )
end


function Base.show( io::Base.IO, X::SpecPrincipalOpen)
  Base.print( io, "Principal open subscheme of \n" )
  Base.print( io, root(X) )
  Base.print( io, "\n" )
  Base.print( io, "defined by\n" )
  Base.print( io,  denoms(X))
end


mutable struct AffSchMorphism{
    Sdom <: Ring, Tdom <: MPolyRing, Udom <: MPolyElem, 
    Scod <: Ring, Tcod <: MPolyRing, Ucod <: MPolyElem 
  } <: SchemeMorphism
  domain::Spec{ Sdom, Tdom, Udom }
  codomain::Spec{ Scod, Tcod, Ucod }
  pullback::AlgHom

  function AffSchMorphism( domain::AffineScheme, codomain::AffineScheme, pullback::AlgHom )
    if coefficient_ring( domain.R ) != coefficient_ring( codomain.R )
      error( "the base rings of the domain and the codomain do not coincide!" )
    end
    #if domain(pullback) != codomain.R || codomain(pullback) != domain.R
      #error( "the domain and codomain of the given ring homomorphism is not compatible with the affine schemes." )
    #end 
    return new{ 
      typeof(coefficient_ring(domain.R)), typeof(domain.R), elem_type(domain.R),
      typeof(coefficient_ring(codomain.R)), typeof(codomain.R), elem_type(codomain.R)
    }( domain, codomain, pullback )
    # return new( domain, codomain, pullback )
  end
end


@doc Markdown.doc"""
    struct Glueing( domain::AffineScheme, codomain::AffineScheme, pullback::AlgHom )

    Maintains the data for the glueing of two affine varieties.
    In practice, domain and codomain will be distinguished open subsets
    of different affine algebras and then pullback will specify
    a concrete isomorphism between them.

    Compared to AffSchMorphism, this structure has an additional field
    in which the inverse morphism can be stored.
"""
mutable struct Isomorphism{
    Sdom <: Ring, Tdom <: MPolyRing, Udom <: MPolyElem,
    Scod <: Ring, Tcod <: MPolyRing, Ucod <: MPolyElem
  } <: SchemeMorphism
  domain::Spec{ Sdom, Tdom, Udom }
  codomain::Spec{ Scod, Tcod, Ucod }
  pullback::AlgHom
  inverse::AlgHom

  function Isomorphism( domain::AffineScheme, codomain::AffineScheme, pullback::AlgHom, inverse::AlgHom )
    if coefficient_ring( domain.R ) != coefficient_ring( codomain.R )
      error( "the base rings of the domain and the codomain do not coincide!" )
    end
    #if domain(pullback) != codomain.R || codomain(pullback) != domain.R
      #error( "the domain and codomain of the given ring homomorphism is not compatible with the affine schemes." )
    #end
    return new{
      typeof(coefficient_ring(domain.R)), typeof(domain.R), elem_type(domain.R),
      typeof(coefficient_ring(codomain.R)), typeof(codomain.R), elem_type(codomain.R)
    }( domain, codomain, pullback, inverse )
  end

  function Isomorphism( domain::AffineScheme, codomain::AffineScheme, pullback::AlgHom )
    if coefficient_ring( domain.R ) != coefficient_ring( codomain.R )
      error( "the base rings of the domain and the codomain do not coincide!" )
    end
    inv_list = [ preimage( pullback, g )[1] for g in gens( codomain( pullback ))]
    inverse = AlgebraHomomorphism( pullback.codomain, pullback.domain, inv_list )
    return new{
      typeof(coefficient_ring(domain.R)), typeof(domain.R), elem_type(domain.R),
      typeof(coefficient_ring(codomain.R)), typeof(codomain.R), elem_type(codomain.R)
    }( domain, codomain, pullback, inverse )
  end
end

# The inclusion of open subsets. 
# Note that this struct can only instantiate open inclusions which 
# are determined by an explicit ring homomorphism. Not all open 
# inclusions (not even of affine schemes) can be described this way!
mutable struct OpenInclusion{
    Sdom <: Ring, Tdom <: MPolyRing, Udom <: MPolyElem, 
    Scod <: Ring, Tcod <: MPolyRing, Ucod <: MPolyElem 
  } <: SchemeMorphism
  domain::Spec{ Sdom, Tdom, Udom }
  codomain::Spec{ Scod, Tcod, Ucod }
  pullback::AlgHom
end

function Oscar.compose( phi::AffSchMorphism, psi::AffSchMorphism )
  if psi.codomain != phi.domain 
    error( "Morphisms of schemes can not be composed." )
  end
  return AffSchMorphism( phi.domain, psi.codomain, compose( psi.pullback, phi.pullback ))
end

function Oscar.compose( phi::Isomorphism, psi::AffSchMorphism )
  if psi.codomain != phi.domain 
    error( "Morphisms of schemes can not be composed." )
  end
  return AffSchMorphism( phi.domain, psi.codomain, compose( psi.pullback, phi.pullback ))
end

function Oscar.compose( phi::AffSchMorphism, psi::Isomorphism )
  if psi.codomain != phi.domain 
    error( "Morphisms of schemes can not be composed." )
  end
  return AffSchMorphism( phi.domain, psi.codomain, compose( psi.pullback, phi.pullback ))
end

function Oscar.compose( phi::Isomorphism, psi::Isomorphism )
  if psi.codomain != phi.domain 
    error( "Morphisms of schemes can not be composed." )
  end
  return Isomorphism( phi.domain, psi.codomain, 
		 compose( psi.pullback, phi.pullback ),
		 compose( phi.inverse, psi.inverse) )
end

function ∘( phi::SchemeMorphism, psi::SchemeMorphism )
  return compose( phi, psi )
end

# A struct for glueing of two affine schemes X = Spec A and Y = Spec B along a common 
# principal open subset U. Then U is given as both, Spec A_f and Spec B_g for 
# elements f ∈  A and g ∈  B, and there is an isomorphism ϕ : A_f →  B_g. The 
# Glueing maintains this information, together with the two inclusions of 
# open sets i : U ↪  X and j : U ↪  Y given by A →  A_f, and B →  B_f, respectively.
mutable struct Glueing
  X::AffineScheme
  Y::AffineScheme
  i::SchemeMorphism
  j::SchemeMorphism
  ϕ::Isomorphism

end

# A struct maintaining the information for a covering 
# of a scheme. A scheme itself can have multiple coverings 
# with a partial order given by refinement. Any object relying 
# on a covering such as, for instance, a vector bundle given by 
# transition functions, can refer to such a covering. 
#
# TODO: Simon argued that such a struct should not be mutable, 
# because the objects depending on it can not keep track 
# of such changes. Technically, however, that means all data 
# for the covering has to be computed before the actual instantiation 
# of the covering. That will involve a lot of temporary variables 
# resembling the final struct and seems a bit tedious. Therefore, 
# I make it mutable for now. 
mutable struct Covering 
  charts::Vector{AffineScheme} 
  glueings::Vector{Glueing}
end

mutable struct CoveredScheme <: Scheme
  coverings::Vector{Covering}
end

abstract type ChowCycle end

# This is the data structure for a cycle in an affine 
# scheme given by a list of coefficients and prime ideals
# TODO: Parametrize by types.
mutable struct AffineCycle{ CoefficientType <: Ring } <: ChowCycle
  parent::AffineScheme
  summands::Tuple{CoefficientType, AbstractAlgebra.Ideal}
end

