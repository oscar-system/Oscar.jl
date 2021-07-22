module Schemes

using Oscar

import AbstractAlgebra.Ring, Oscar.AlgHom

export AffineScheme, affine_space, localize, MorphismAffSch

abstract type Scheme end

mutable struct AffineScheme{S <: Ring, T <: MPolyRing, U <: MPolyElem} <: Scheme
  k::S			# the base ring (usually a field) of definition for the scheme
  R::T		  	# the ambient polynomial ring to model this affine scheme
  I::MPolyIdeal{U}	# The ideal in R defining the scheme
#mutable struct AffineScheme{S <: Ring, T <: MPolyRing, U <: MPolyElem} <: Scheme
# k::S			# the base ring (usually a field) of definition for the scheme
# R::T		  	# the ambient polynomial ring to model this affine scheme
# I::MPolyIdeal{U}	# The ideal in R defining the scheme

  function AffineScheme(k::S, R::T, I::MPolyIdeal{U}) where{S <: Ring, T <:MPolyRing , U <: MPolyElem}
    if k != base_ring(R)
      error( "Base ring of the affine scheme does not coincide with the base ring of the associated algebra." )
    end
    # TODO: Implement further plausibility checks to be performed at runtime.
    return new{S, T, U}(k, R, I)
  end

  function AffineScheme( k::S, R::T ) where{S <: Ring, T <:MPolyRing}
    I = ideal(R, zero(R))
    return new{S, T, elem_type(T)}(k, R, I)
  end

  function AffineScheme( R::T ) where{T <: MPolyRing}
    I = ideal(R, zero(R))
    k = base_ring( R )
    #@show typeof(k), T, elem_type(T)
    return new{typeof(k), T, elem_type(T)}(k, R, I)
  end

  function AffineScheme( R::T, I::MPolyIdeal{U} ) where { T <: MPolyRing, U <: MPolyElem }
    k = base_ring(R)
    return new{typeof(k), T, elem_type(T)}(k, R, I)
  end

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

mutable struct MorphismAffSch{ 
    Sdom <: Ring, Tdom <: MPolyRing, Udom <: MPolyElem, 
    Scod <: Ring, Tcod <: MPolyRing, Ucod <: MPolyElem 
  } 
  domain::AffineScheme{ Sdom, Tdom, Udom } 
  codomain::AffineScheme{ Scod, Tcod, Ucod } 
  pullback::AlgHom
  
  function MorphismAffSch( domain::AffineScheme, codomain::AffineScheme, pullback::AlgHom )
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
    phi = MorphismAffSch( X, Y, pullback )
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
    phi = MorphismAffSch( X, Y, pullback )
  end
  return Y, phi
end

end
