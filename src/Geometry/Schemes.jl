module Schemes

using Oscar
using Markdown

import AbstractAlgebra.Ring, Oscar.AlgHom, Oscar.compose
import Base: ∘

export AffineScheme, affine_space, AffSchMorphism

export localize, base_ring, defining_ideal, ambient_ring

abstract type Scheme end
abstract type SchemeMorphism end

mutable struct AffineScheme{S <: Ring, T <: MPolyRing, U <: MPolyElem} <: Scheme # this is the mutuble struct 
  # the basic fields 
  k::S			# the base ring (usually a field) of definition for the scheme
  R::T		  	# the ambient polynomial ring to model this affine scheme
  I::MPolyIdeal{U}	# The ideal in R defining the scheme

  function AffineScheme(k::S, R::T, I::MPolyIdeal{U} ) where{
			S <: Ring, T <:MPolyRing , U <: MPolyElem}
    if k != base_ring(R)
      error( "Base ring of the affine scheme does not coincide with the base ring of the associated algebra." )
    end
    # TODO: Implement further plausibility checks to be performed at runtime.
    return new{S, T, U}(k, R, I )
  end
end

function base_ring(A::AffineScheme)
  return A.k
end

function ambient_ring(A::AffineScheme)
  return A.R
end

function defining_ideal(A::AffineScheme)
  return A.I
end

mutable struct PrincipalOpenAffSubSch{S <: Ring, T<:MPolyRing, U<:MPolyElem}
  parent::Union{AffineScheme{S,T,U},PrincipalOpenAffSubSch{S,T,U}}
  denom::MPolyElem{U}  # element of the "ambient" polynomial ring of the root
  # these fields are only initialized if needed and used to comply to the AffineScheme interface
  k::S  
  R::T
  I::MPolyIdeal{U} 
  
  function PrincipalOpenAffSubScheme(parent, denom)
    x = new{S,T,U}() 
    x.parent = parent
    x.denom = denom
  end 
end
  
# outer constructors
function AffineScheme( k::S, R::T ) where{S <: Ring, T <:MPolyRing}
  I = ideal(R, zero(R))
  return AffineScheme(k, R, I )
end

@doc Markdown.doc"""
    AffineScheme(R::MPolyRing) -> AffineScheme

Return the affine scheme corresponding to the ring $R$.
"""
function AffineScheme( R::T ) where{T <: MPolyRing}
  I = ideal(R, zero(R))
  k = base_ring( R )
  return AffineScheme(k, R, I )
end

function AffineScheme( R::T, I::MPolyIdeal{U} ) where{ T <: MPolyRing, U <: MPolyElem }
  k = base_ring(R)
  return AffineScheme(k, R, I )
end

# Construct affine n-space over the ring k.
function affine_space( k::Ring, n::Int, name::String="x" )
  R, x = PolynomialRing( k, name => (1:n))
  return AffineScheme( R )
end

function Base.show( io::Base.IO, X::AffineScheme )
  Base.print( io, "Affine scheme over " )
  Base.print( io, X.k )
  Base.print( io, "\n" ) 
  Base.print( io, "given as the quotient of the polynomial ring \n" )
  Base.print( io, X.R )
  Base.print( io, "\nby the ideal \n" )
  Base.print( io, X.I )
end

mutable struct AffSchMorphism{
    Sdom <: Ring, Tdom <: MPolyRing, Udom <: MPolyElem, 
    Scod <: Ring, Tcod <: MPolyRing, Ucod <: MPolyElem 
  } <: SchemeMorphism
  domain::AffineScheme{ Sdom, Tdom, Udom } 
  codomain::AffineScheme{ Scod, Tcod, Ucod } 
  pullback::AlgHom
  
  function AffSchMorphism( domain::AffineScheme, codomain::AffineScheme, pullback::AlgHom )
    if base_ring( domain.R ) != base_ring( codomain.R )
      error( "the base rings of the domain and the codomain do not coincide!" )
    end
    #if domain(pullback) != codomain.R || codomain(pullback) != domain.R
      #error( "the domain and codomain of the given ring homomorphism is not compatible with the affine schemes." )
    #end 
    return new{ 
      typeof(base_ring(domain.R)), typeof(domain.R), elem_type(domain.R), 
      typeof(base_ring(codomain.R)), typeof(codomain.R), elem_type(codomain.R) 
    }( domain, codomain, pullback )
    # return new( domain, codomain, pullback )
  end
end

# construct the affine scheme given by localization 
# and return that new scheme together with the algebra 
# morphism for the open inclusion
function localize( X, f::MPolyElem, name::String="denom" ) 
  # Check whether there already exists a variable named 'denom'
  var_symb = symbols( X.R )
  if Symbol( name ) in var_symb 
    # If yes, recycle it.
    # We implicitly assume that the variable called "denom" is 
    # coming from a previous localization and the defining ideal 
    # I in X.R is of the form I = I' + ⟨1-denom⋅f⟩ for some f 
    # and some ideal I' not involving the variable denom.
    y = gens( X.R )
    unit_index = indexin( [Symbol(name)], var_symb )[1]
    unit = y[unit_index]
    var_images = gens( X.R )
    var_images[unit_index] = unit*f
    pullback = AlgebraHomomorphism( X.R, X.R, var_images )
    I = ideal( pullback.( gens(X.I) ))
    Y = AffineScheme( X.R, I )
    phi = AffSchMorphism( X, Y, pullback )
  else
    # If not, then localize by introducing a new denominator
    # This will be the first variable in the new ring, since we 
    # expect an elimination ordering for 'denom' to be imposed 
    # on this ring for saturation purposes.
    B, y = PolynomialRing( base_ring(X.R), vcat( [name], String.( symbols( X.R )) ))
    n = length( gens( X.R ))
    pullback = AlgebraHomomorphism( X.R, B, y[2:n+1] )
    I = ideal( B, [ pullback( g ) for g in gens(X.I) ] )
    I = I + ideal( B, one(B) - y[1]*pullback(f) )
    Y = AffineScheme( B, I )
    phi = AffSchMorphism( X, Y, pullback )
  end
  return Y, phi
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
  domain::AffineScheme{ Sdom, Tdom, Udom } 
  codomain::AffineScheme{ Scod, Tcod, Ucod } 
  pullback::AlgHom
  inverse::AlgHom
  
  function Isomorphism( domain::AffineScheme, codomain::AffineScheme, pullback::AlgHom, inverse::AlgHom )
    if base_ring( domain.R ) != base_ring( codomain.R )
      error( "the base rings of the domain and the codomain do not coincide!" )
    end
    #if domain(pullback) != codomain.R || codomain(pullback) != domain.R
      #error( "the domain and codomain of the given ring homomorphism is not compatible with the affine schemes." )
    #end 
    return new{ 
      typeof(base_ring(domain.R)), typeof(domain.R), elem_type(domain.R), 
      typeof(base_ring(codomain.R)), typeof(codomain.R), elem_type(codomain.R) 
    }( domain, codomain, pullback, inverse )
  end

  function Isomorphism( domain::AffineScheme, codomain::AffineScheme, pullback::AlgHom )
    if base_ring( domain.R ) != base_ring( codomain.R )
      error( "the base rings of the domain and the codomain do not coincide!" )
    end
    inv_list = [ preimage( pullback, g )[1] for g in gens( codomain( pullback ))]
    inverse = AlgebraHomomorphism( pullback.codomain, pullback.domain, inv_list )
    return new{ 
      typeof(base_ring(domain.R)), typeof(domain.R), elem_type(domain.R), 
      typeof(base_ring(codomain.R)), typeof(codomain.R), elem_type(codomain.R) 
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
  domain::AffineScheme{ Sdom, Tdom, Udom } 
  codomain::AffineScheme{ Scod, Tcod, Ucod } 
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

end
