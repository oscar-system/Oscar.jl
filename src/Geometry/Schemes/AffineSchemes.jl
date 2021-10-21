########################################################################
# Spectra of polynomial algebras                                       #
########################################################################

@doc Markdown.doc"""
    Spec{BaseRingType, RingType, RingElemType} <: AffineScheme{BaseRingType, RingType, RingElemType}

An affine scheme ``X = Spec R/I`` with ``R = k[x₁,…,xₙ]`` a free 
polynomial algebra of type `RingType` over a base ring ``k`` of type 
`BaseRingType` and ``I ⊂ R`` a finitely generated ideal 
with elements of type ``RingElemType``.
"""
mutable struct Spec{BaseRingType, RingType, RingElemType} <: AffineScheme{BaseRingType, RingType, RingElemType}
  # the basic fields 
  k::BaseRingType		# the base ring (usually a field) of definition for the scheme
  R::RingType 			# the ambient polynomial ring to model this affine scheme
  I::MPolyIdeal{RingElemType}	# The ideal in R defining the scheme

  # fields for caching
  name::String # the name of this scheme for printing

  function Spec(k::BaseRingType, R::RingType, I::MPolyIdeal{RingElemType} ) where{
			BaseRingType <: Ring, RingType <:MPolyRing , RingElemType <: MPolyElem}
    k == coefficient_ring(R) || error("Base ring of the affine scheme does not coincide with the base ring of the associated algebra")
    base_ring(I) == R || error("Ideal does not belong to the ring")
    return new{BaseRingType, RingType, RingElemType}(k, R, I )
  end
end


### Getter functions ###################################################

@Markdown.doc """
    base_ring(A::Spec{BaseRingType, RingType, RingElemType}) where {BaseRingType, RingType, RingElemType}

Returns the base ring ``k`` over which the affine scheme ``A`` is defined.
"""
function base_ring(A::Spec{BaseRingType, RingType, RingElemType}) where {BaseRingType, RingType, RingElemType}
  return A.k
end


@Markdown.doc """
    ambient_ring(A::Spec{BaseRingType, RingType, RingElemType}) where {BaseRingType, RingType, RingElemType}

Returns the "ambient ring" ``R`` of the affine scheme ``A``, where ``A = Spec R/I``.
"""
function ambient_ring(A::Spec{BaseRingType, RingType, RingElemType}) where {BaseRingType, RingType, RingElemType}
  return A.R
end

@Markdown.doc """
     defining_ideal(A::Spec{BaseRingType, RingType, RingElemType}) where {BaseRingType, RingType, RingElemType}
Returns the "defining ideal" ``I`` of the affine scheme ``A``, where ``A = Spec R/I``.
"""
function defining_ideal(A::Spec{BaseRingType, RingType, RingElemType}) where {BaseRingType, RingType, RingElemType}
  return A.I
end


@Markdown.doc """
    denoms(A::Spec{BaseRingType, RingType, RingElemType}) where {BaseRingType, RingType, RingElemType}

This function simply extends the `denoms(D::SpecPrincipalOpen)` so that it 
can conveniently be used on all of `AffineScheme`. It returns an empty list.
"""
function denoms(A::Spec{BaseRingType, RingType, RingElemType}) where {BaseRingType, RingType, RingElemType}
  return typeof(A.R)[]
end


### outer constructors ################################################

@Markdown.doc """
    Spec( X::Spec{BaseRingType, RingType, RingElemType} ) where {BaseRingType, RingType, RingElemType}

Returns a copy of the affine scheme ``X``.
"""
Spec(X::Spec{BaseRingType, RingType, RingElemType}) where {BaseRingType, RingType, RingElemType} = Spec(base_ring(X), ambient_ring(X), deepcopy(defining_ideal(X)))

# the routine for deepcopies internally used by julia. 
# Note that rings don't need to be deepcopied. Ideals on 
# the other hand, posess information that can be altered. 
Base.deepcopy_internal(X::Spec{BaseRingType, RingType, RingElemType}, dict::IdDict) where {BaseRingType, RingType, RingElemType} = Spec(X)


@Markdown.doc """
    Spec( k::BaseRingType, R::RingType ) where{BaseRingType <: Ring, RingType <:MPolyRing}

Returns the scheme X = Spec R with R an instance of MPolyRing over k.
"""
function Spec( k::BaseRingType, R::RingType ) where{BaseRingType <: Ring, RingType <:MPolyRing}
  I = ideal(R, zero(R))
  return Spec(k, R, I )
end


@doc Markdown.doc"""
    Spec(R::RingType) where{RingType <: MPolyRing}

Return the affine scheme corresponding to the ring ``R/⟨0⟩``.
"""
function Spec(R::RingType) where{RingType <: MPolyRing}
  I = ideal(R, zero(R))
  k = coefficient_ring(R)
  return Spec(k, R, I)
end


@doc Markdown.doc"""
    Spec(R::RingType, I::MPolyIdeal{RingElemType}) where{ RingType <: MPolyRing, RingElemType <: MPolyElem }

Return the affine scheme corresponding to the ring ``R/I``.
"""
function Spec(R::RingType, I::MPolyIdeal{RingElemType}) where{ RingType <: MPolyRing, RingElemType <: MPolyElem }
  k = coefficient_ring(R)
  return Spec(k, R, I)
end


function Spec(I::MPolyIdeal{RingElemType}) where {RingElemType <: MPolyElem}
  R = base_ring(I)
  k = coefficient_ring(R)
  return Spec(k, R, I)
end


function Spec(Q::MPolyQuo{RingElemType}) where {RingElemType <: MPolyElem}
  R = base_ring(Q)
  k = coefficient_ring(R)
  I = modulus(Q)
  return Spec(k, R, I)
end


### other constructors ################################################

@Markdown.doc """
    function affine_space( k::Ring, n::Int; var_name::String="x" )

Return the Affine scheme 𝔸ⁿ over the base ring k. An optional argument 
`var_name` can be 
passed on to determine how the variables will be printed.
"""
function affine_space( k::BaseRingType, n::Int; var_name::String="x" ) where {BaseRingType <: Ring}
  R, x = PolynomialRing( k, var_name => (1:n))
  return Spec(R)
end


@Markdown.doc """
    subscheme(X::Spec{BaseRingType, RingType, RingElemType}, J::MPolyIdeal{RingElemType}) where {BaseRingType, RingType, RingElemType}

Returns the subscheme of ``X = Spec R/I`` defined by the ideal ``I + J`` in the ring ``R``.
"""
function subscheme(X::Spec{BaseRingType, RingType, RingElemType}, J::MPolyIdeal{RingElemType}) where {BaseRingType, RingType, RingElemType}
  base_ring(J) == ambient_ring(X) || error( "Ideal does not belong to the ambient ring of the variety" )
  return Spec(base_ring(X), ambient_ring(X), defining_ideal(X) + J)
end


@Markdown.doc """
    subscheme(X::Spec{BaseRingType, RingType, RingElemType}, f::RingElemType) where {BaseRingType, RingType, RingElemType}

Returns the subscheme of ``X = Spec R/I`` defined by the ideal`` I + ⟨f⟩`` in the ring ``R``.
"""
function subscheme(X::Spec{BaseRingType, RingType, RingElemType}, f::RingElemType) where {BaseRingType, RingType, RingElemType}
  return subscheme(X, ideal(ambient_ring(X), f))
end


### routines for printing #############################################

@Markdown.doc """
    Base.show(io::Base.IO, X::AffineScheme)

This prints the information of X to the stream io. Whenever 
the field X.name is set, the output will be abbreviated to 
just the given name of X.
"""
function Base.show( io::Base.IO, X::AffineScheme )
  if isdefined( X, :name )
    Base.print( io, X.name )
    return
  end
  Base.print( io, "Affine scheme over " )
  Base.print( io, base_ring(X) )
  Base.print( io, "\n" ) 
  Base.print( io, "given as the quotient of the polynomial ring \n" )
  Base.print( io, ambient_ring(X) )
  Base.print( io, "\nby the ideal \n" )
  Base.print( io, defining_ideal(X) )
end

@Markdown.doc """
    set_name!(X::AffineScheme, name::String)

Each affine scheme `X` can be assigned an optional `name`. 
This will then be used for printing, once the field is set.
"""
function set_name!( X::AffineScheme, name::String )
  X.name = name
end


#######################################################################
# Morphisms of Spectra                                                #
#######################################################################


