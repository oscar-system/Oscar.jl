import AbstractAlgebra.Ring, Oscar.AlgHom, Oscar.compose, AbstractAlgebra.Generic.Frac
import Base: ∘
import Oscar: base_ring
import AbstractAlgebra: FPModule, FPModuleElem
import Base.copy
import AbstractAlgebra.Generic.MatSpaceElem

include( "./Multiindices.jl" )
using Oscar.Multiindices

include( "./Misc.jl" )
using Oscar.Misc

export AffineScheme, Spec, SpecPrincipalOpen
export base_ring, ambient_ring, defining_ideal, imgs_frac, pullback, pullback_from_parent, pullback_from_root, inclusion_in_parent, inclusion_in_root, set_name!, inverted_element, identity_morphism, denoms, inverses, divide_by_units
export affine_space, localize, subscheme

export AffSchMorphism
export domain, codomain, pullback

export Glueing
export first_patch, second_patch, overlap, first_inclusion, second_inclusion, glueing_isomorphism

export Covering
export patches, glueings, _add_patch!

export Refinement

export CoveredScheme
export coverings

export CovSchMorphism

export AbstractCoherentSheaf

export IdealSheaf
export parent, covering, ideals

export CoherentSheaf
export parent, covering, modules

export LineBundle, VectorBundle
export parent, covering, transitions, rank, OO, direct_sum

export VectorBundleSection, LineBundleSection
export tautological_sections

export projectivize
export cartesian_product

@Markdown.doc """
    abstract type Scheme{ S <: Ring }

The type parameter S is the type of the base ring k over which the Scheme is defined.
"""
abstract type Scheme{ S <: Ring }end
@doc Markdown.doc"""
    abstract type AffineScheme{S, T <: MPolyRing, U <: MPolyElem} <: Scheme{S}

The type parameters are: 
 * T: The type of the ambient polynomial ring R
 * U: The type of the elements in R
"""
abstract type AffineScheme{S, T <: MPolyRing, U <: MPolyElem} <: Scheme{S} end
abstract type SchemeMorphism end

abstract type SchematicPoint end
abstract type Sheaf end
abstract type AbstractCoherentSheaf <: Sheaf end

@doc Markdown.doc"""
    mutable struct Spec{S,T,U} <: AffineScheme{S,T,U}

Models an affine scheme X = Spec R/I with R = k[x₁,…,xₙ] a free 
polynomial algebra of type T over a base ring k of type S and I ⊂ R a finitely generated ideal 
with elements of type U.
"""
mutable struct Spec{S,T,U} <: AffineScheme{S,T,U}
  # the basic fields 
  k::S			# the base ring (usually a field) of definition for the scheme
  R::T		  	# the ambient polynomial ring to model this affine scheme
  I::MPolyIdeal{U}	# The ideal in R defining the scheme

  # fields for caching
  name::String # the name of this scheme for printing

  function Spec(k::S, R::T, I::MPolyIdeal{U} ) where{
			S <: Ring, T <:MPolyRing , U <: MPolyElem}
    if k != coefficient_ring(R)
      error("Base ring of the affine scheme does not coincide with the base ring of the associated algebra")
    end
    # TODO: Implement further plausibility checks to be performed at runtime.
    return new{S, T, U}(k, R, I )
  end
end

# The copy constructors.
# In this case, it's easy, since non of the used data structures 
# can be modified from the outside with "Fernwirkung" affecting 
# our instance of Spec. 
#
# To explain what I mean by "Fernwirkung": Say, your struct stores 
# a vector v as a member variable and v is passed on to the constructor. 
# Simply storing v is wrong, because it is an object that could potentially
# also be accessed from other parts of the program. If v is altered, then 
# also the member variable of our Spec is altered! Thus in this case, 
# we would really need to copy v before storing it in the Spec.
#
# Such behaviour is not expected for instances of rings or ideals. 
# As far as I understand, these objects are supposed to behave like 
# static ones to the user. But we should watch out for exceptions and 
# bugs here.
@Markdown.doc """
    Spec( X::Spec ) -> Spec

Returns a copy of the affine scheme X.
"""
Spec( X::Spec ) = Spec( X.k, X.R, X.I )

Base.copy( X::Spec ) = Spec(X)

###################################################################################
# Getter functions
#
# No caching is needed in this case, since all these variables need to be assigned 
# at instantiation
#
@Markdown.doc """
    base_ring(A::Spec)

Returns the base ring k over which the affine scheme A is defined.
"""
function base_ring(A::Spec)
  return A.k
end

@Markdown.doc """
    ambient_ring(A::Spec)

Returns the "ambient ring" R of the affine scheme A, where A = Spec R/I.
"""
function ambient_ring(A::Spec)
  return A.R
end

@Markdown.doc """
     defining_ideal(A::Spec)
Returns the "defining ideal" I of the affine scheme A, where A = Spec R/I.
"""
function defining_ideal(A::Spec)
  return A.I
end

@Markdown.doc """
    denoms(A::Spec)

This function simply extends the denoms(D::SpecPrincipalOpen) so that it 
can conveniently be used on all of AffineScheme. It returns an empty list.
"""
function denoms(A::Spec)
  return typeof(A.R)[]
end

@Markdown.doc """
    subscheme( X::Spec, J::MPolyIdeal ) 

Returns the subscheme of X = Spec R/I defined by the ideal I + J in the ring R.
"""
function subscheme( X::Spec, J::MPolyIdeal ) 
  if base_ring(J) != X.R 
    error( "Ideal does not belong to the ambient ring of the variety" )
  end
  Y = Spec( X.k, X.R, X.I )
  Y.I = Y.I + J
  return Y
end

@Markdown.doc """
    subscheme( X::Spec, f::MPolyElem )

Returns the subscheme of X = Spec R/I defined by the ideal I + ⟨f⟩ in the ring R.
"""
function subscheme( X::Spec, f::MPolyElem )
  return subscheme( X, ideal( ambient_ring(X), f ))
end

@Markdown.doc """
    mutable struct SpecPrincipalOpen{S,T,U} <: AffineScheme{S,T,U}

A principal open subset D ⊂ X of an affine scheme X = Spec R/I (with R of type T) over a base 
ring k (of type S), arising as the complement of the zero locus of a polynomial in R (of type U). 

The principal open subsets are implemented in a recursive way: Each instance 
D has a parent D' which is either again of type SpecPrincipalOpen, or of type Spec. 
This ancestry relation forms a tree structure, the so-called *localization tree*,
with a unique instance X = Spec R/I of type Spec at 
its root. Then D is obtained from D' as the complement of the zero locus of some 
function f in the ambient ring R of the root X. At each such 'localization'-step, 
this function f is stored in the field D.denom. 

The implementation is lazy in the sense that a priori, each instance D of SpecPrincipalOpen 
merely stores the information about its parent D' and the element f ∈ R that has been inverted 
on D compared to its parent D'. Whenever the user asks for it, an explicit algebra 
D = R[t]/(I + ⟨1 - t⋅∏ᵢ fᵢ⟩) 
can be constructed, where the fᵢ∈ R are the elements which have been inverted along the 
path from D to the root X in the localization tree. 
"""
mutable struct SpecPrincipalOpen{S,T,U} <: AffineScheme{S,T,U}
  parent::AffineScheme{S,T,U}
  denom::U  # element of the "ambient" polynomial ring of the root
  	    # which has been localized in the parent to get this scheme.
  d::U # the element that really had to be inverted given the previous localizations.

  # Fields for caching. These provide the data as a Spec for this 
  # affine scheme
  k::S  # the ground field
  R::T	# the ambient polynomial ring
  I::MPolyIdeal{U} # the defining ideal
  pullbackFromParent::AlgHom
  pullbackFromRoot::AlgHom
  u::U # The inverse of denom in R
  name::String
  loc_var_name::String
  #inclusionInParent::AffSchMorphism
  #inclusionInRoot::AffSchMorphism

  
  function SpecPrincipalOpen(orig::Union{Spec{S,T,U},SpecPrincipalOpen{S,T,U}}, denom::U; loc_var_name="t" ) where {S <: Ring, T<:MPolyRing, U<:MPolyElem}
    x = new{S,T,U}() 
    x.parent = orig
    # TODO: Implement a plausibility check: Does denom belong to the coordinate ring of the root?
    parent(denom) == ambient_ring(root(x)) || error( "denominator has the wrong ring as parent" )
    x.denom = denom
    x.loc_var_name=loc_var_name
    return x 
  end 

  function SpecPrincipalOpen( k::S, R::T, I::MPolyIdeal{U}; loc_var_name="t" ) where {S <: Ring, T<:MPolyRing, U<:MPolyElem}
    if k != base_ring(R) 
      error( "basering not compatible with ambient ring" )
    end
    if base_ring(I) != R
      error( "ideal does not live in the given ambient ring" )
    end
    X = new{S,T,U}(k,R,I,loc_var_name)
    return X
  end
end

# Copy constructor
@Markdown.doc """
    SpecPrincipalOpen( X::SpecPrincipalOpen )

Returns a copy of X
"""
function SpecPrincipalOpen( X::SpecPrincipalOpen )
  Y = SpecPrincipalOpen( X.k, X.R, X.I )
  if isdefined( X, :denom )
    Y.denom = X.denom # this is safe to copy. It's excempt from "Fernwirkung" (tested).
  end
  # The following two variables are being cached. So we probably don't want 
  # to copy them, because the need to be build on the new instance and not 
  # accidentally refer to the member-variables of the original instance. 
  # We should talk about this. 
#  if isdefined( X, :pullbackFromParent )
#    Y.pullbackFromParent = X.pullbackFromParent 
#  end
#  if isdefined( X, :pullbackFromRoot )
#    Y.pullbackFromRoot = X.pullbackFromRoot
#  end
  if isdefined( X, :u )
    Y.u = X.u
  end
  return Y
end

Base.copy( X::SpecPrincipalOpen ) = SpecPrincipalOpen( X )

@Markdown.doc """
    Spec( D::SpecPrincipalOpen ) 
    
Returns an instance of Spec isomorphic to D. This uses 
the Rabinowitsch version for the localizations and forgets 
about the ancestry in the localization tree.
"""
function Spec( D::SpecPrincipalOpen ) 
  return Spec( base_ring(D), ambient_ring(D), defining_ideal(D) )
end
###############################################################################
# Getter functions

@Markdown.doc """
    localize(X::AffineScheme{S,T,U}, f::U; loc_var_name="t" ) where {S,T,U}

Returns an instance D of SpecPrincipalOpen where D = X ∖ Z(f).
The optional argument `loc_var_name` can be used to specify a name 
for printing the additional variable used in Rabinowitsch's trick 
for localizations.
"""
function localize(X::AffineScheme{S,T,U}, f::U; loc_var_name="t" ) where {S,T,U}
  return SpecPrincipalOpen(X, f; loc_var_name )
end

@Markdown.doc """
    localize(X::AffineScheme{S,T,U}, v::Vector{U}; loc_var_name="t" ) where {S,T,U}

Apply `localize(X,f)` recursively to the list of functions f in v.
"""
function localize(X::AffineScheme{S,T,U}, v::Vector{U}; loc_var_name="t" ) where {S,T,U}
  D = X
  for f in v 
    D = localize(D, f, loc_var_name)
  end
  return D
end

@Markdown.doc """
    function parent(S::AffineScheme)

Return the parent of S in the tree structure for Spec and SpecPrincipalOpen.
"""
function parent(S::AffineScheme)
  return S
end

# Check whether X is a parent of U and return the inclusion morphism 
# if true. Otherwise, return nothing.
@Markdown.doc """
    function is_parent_of( X::AffineScheme, U::SpecPrincipalOpen )

Returns true whenever X occurs as a parent of U in the tree structure 
for Spec and SpecPrincipalOpen.
"""
function is_parent_of( X::AffineScheme, U::SpecPrincipalOpen )
  i = identity_morphism(U)
  if X == U
    return i
  end
  V = U
  while typeof(V)<:SpecPrincipalOpen
    W = parent(V)
    i = compose( inclusion_in_parent(V), i )
    if W == X 
      return i
    end
    V = W
  end
  return nothing
end
    

@Markdown.doc """
    function parent(D::SpecPrincipalOpen)

Returns the parent of D in the localization tree.
"""
function parent(D::SpecPrincipalOpen)
  return D.parent
end 

@Markdown.doc """
    function root(D::SpecPrincipalOpen)

Returns the root of D in the localization tree.
"""
function root(D::SpecPrincipalOpen)
  return root(parent(D))
end

@Markdown.doc """
    root(A::Spec)

Returns A itself since this is necessarily the root in the localization tree.
"""
root(A::Spec) = A

@Markdown.doc """
    function base_ring(D::SpecPrincipalOpen)

Returns the base ring k of D where D = Spec (R/I)[(f₁)⁻¹,…,(fᵣ)⁻¹] with 
R a free polynomial ring over k.
"""
function base_ring(D::SpecPrincipalOpen)
  if isdefined(D,:k)
    return D.k
  end
  k = base_ring(root(D))
  D.k = k
  return k
end

@Markdown.doc """
    function inverted_element(D::SpecPrincipalOpen)

Returns a representative u of the inverse 1/f in the 
ambient ring R[t] of D where f is the element of R that 
has last been inverted in the localization tree.
"""
function inverted_element(D::SpecPrincipalOpen)
  if isdefined( D, :u )
    return D.u
  end
  ambient_ring(D)
  return D.u
end

@Markdown.doc """
    function denom(D::SpecPrincipalOpen)

Returns the element f of the ambient ring of the root 
of D in the localization tree that has last been inverted.
This is the inverse element of denom(D).
"""
function denom(D::SpecPrincipalOpen)
  if isdefined( D, :denom )
    return D.denom
  end
  return nothing
end

@Markdown.doc """
    function root_variables( D::SpecPrincipalOpen )

Returns a vector containing the images of the variables of the 
ambient ring R of the root of D in the localization tree as 
elements in the ambient ring R[t] of D.
"""
function root_variables( D::SpecPrincipalOpen )
  phi = pullback_from_root(D)
  return [ phi(x) for x in gens(ambient_ring(root(D))) ]
end
  
@Markdown.doc """
    divide_by_units( h::ElemTypeAmbientRing, 
        D::SpecPrincipalOpen{BaseRingType,AmbientRingType,ElemTypeAmbientRing}
      ) where {
        BaseRingType<:Ring, AmbientRingType<:Ring, ElemTypeAmbientRing
      }

Takes a multivariate polynomial h = ∑ₖ pₖ⋅uᵏ in the ambient ring R[f⁻¹] ≅ R[u]/⟨1-u⋅f⟩
and writes it as h = uᵈ⋅p with p ∈ R for some d ∈ ℕ₀. It then returns the triple
(f, d, p).
"""
function divide_by_units( h::ElemTypeAmbientRing, D::SpecPrincipalOpen{BaseRingType,AmbientRingType,ElemTypeAmbientRing}) where {BaseRingType<:Ring, AmbientRingType<:Ring, ElemTypeAmbientRing}
  parent(h) == ambient_ring(D) || error("polynomial does not belong to the ambient ring of the scheme")
  X = root(D)
  R = ambient_ring(X) 
  R_f = ambient_ring(D) 
  vars = gens(R_f)
  f = prod(div_denoms(D))
  p_k = coefficients(h, length(vars)) # assume the last variable is the one used for Rabinowitsch's trick.
  p = zero(R_f)
  for a in p_k
    p = p*pullback_from_root(D)(f) + a
  end
  return (f, length(p_k)-1, evaluate(p, push!(gens(R), zero(R))))
end


# Setter functions

@Markdown.doc """
    function set_name!( X::AffineScheme, name::String )

Each affine scheme X can be assigned an optional 'name'. 
This will then be used for printing, once the field is set.
"""
function set_name!( X::AffineScheme, name::String )
  X.name = name
end

@Markdown.doc """
    function subscheme( 
        X::SpecPrincipalOpen{S,T,U}, 
        J::MPolyIdeal{U} 
      ) where{ S<:Ring, T<:MPolyRing, U<:MPolyElem }

Returns the subscheme of X = Spec (R/I)[f⁻¹] defined by the ideal I + J.
"""
function subscheme( 
    X::SpecPrincipalOpen{S,T,U}, 
    J::MPolyIdeal{U} 
  ) where{ S<:Ring, T<:MPolyRing, U<:MPolyElem }
  if base_ring(J) != ambient_ring(root(X))
    error( "Ideal does not belong to the ambient ring of the variety" )
  end
  Y_root = subscheme( root(X), J )
  v = denoms(X)
  Y = Y_root
  for f in v
    Y = localize( Y, f )
  end
  return Y
end

@Markdown.doc """
    function subscheme( 
        X::SpecPrincipalOpen{S,T,U}, 
        f::U
      ) where{ S<:Ring, T<:MPolyRing, U<:MPolyElem }

Returns the subscheme of X = Spec (R/I)[f⁻¹] defined by the ideal I + ⟨f⟩.
"""
function subscheme( 
    X::SpecPrincipalOpen{S,T,U}, 
    f::U
  ) where{ S<:Ring, T<:MPolyRing, U<:MPolyElem }
  return subscheme( X, ideal( ambient_ring(root(X)), f ))
end

@Markdown.doc """
    function denoms(D::SpecPrincipalOpen{S,T,U}) where {S <: Ring, T <: MPolyRing, U<:MPolyElem}

Returns a collection of all the elements fᵢ that have been inverted along the 
path in the localization tree from D up to its root.
"""
function denoms(D::SpecPrincipalOpen{S,T,U}) where {S <: Ring, T <: MPolyRing, U<:MPolyElem}
  result = U[]
  P = D
  while typeof(P) <: SpecPrincipalOpen
    # Make sure that the denominators are collected in the same manner 
    # as the variables for the inverses appear.
    push!( result, P.denom )
    P = parent(P)
  end
  return result
end 

@Markdown.doc """
    denoms_up_to(D::SpecPrincipalOpen{S,T,U}, X::AffineScheme{S,T,U}) where {
	S <: Ring, T <: MPolyRing, U<:MPolyElem}

Returns the list of functions that have been localized in D compared to X, assuming 
X is a parent of D in the localization tree. 
"""
function denoms_up_to(D::SpecPrincipalOpen{S,T,U}, X::AffineScheme{S,T,U}) where {
	S <: Ring, T <: MPolyRing, U<:MPolyElem}
  result = U[]
  P = D
  while D != X
    if typeof(D)<:Spec 
      error( "the second argument is not a parent of the first" )
    end
    # Make sure that the denominators are collected in the same manner 
    # as the variables for the inverses appear.
    push!( result, P.denom )
    P = parent(P)
  end
  return result
end 

@Markdown.doc """
    function div_denoms(D::SpecPrincipalOpen{S,T,U}) where {S <: Ring, T <: MPolyRing, U<:MPolyElem}

When localizing, we internally do a cheap check for obvious redundance: When inverting f = u ⋅ f' 
where u is already a unit, we may replace f by f' in the localization process. 
This routine returns the actual values for the f's.
"""
function div_denoms(D::SpecPrincipalOpen{S,T,U}) where {S <: Ring, T <: MPolyRing, U<:MPolyElem}
  result = U[]
  P = D
  while typeof(P) <: SpecPrincipalOpen
    # Make sure that the denominators are collected in the same manner 
    # as the variables for the inverses appear.
    ambient_ring(P)
    push!( result, P.d )
    P = parent(P)
  end
  return result
end 

@Markdown.doc """
    function inverses(D::SpecPrincipalOpen{S,T,U}) where {S <: Ring, T <: MPolyRing, U<:MPolyElem}

Returns the inverses uᵢ= 1/fᵢof the elements that have been localized along the path 
in the localization tree from D up to its root.
"""
function inverses(D::SpecPrincipalOpen{S,T,U}) where {S <: Ring, T <: MPolyRing, U<:MPolyElem}
  result = [inverted_element(D)]
  R = ambient_ring(D)
  phi = pullback_from_parent(D)
  P = parent(D)
  while typeof(P) <: SpecPrincipalOpen
    push!(result, phi(inverted_element(P)))
    phi = compose( pullback_from_parent(P), phi )
    P = parent(P)
  end
  return result
end 
 
@Markdown.doc """
    function ambient_ring(D::SpecPrincipalOpen)

Returns the ambient ring R[t] of D where D = Spec (R/I)[(f₁)⁻¹,…,(fᵣ)⁻¹] with 
R a free polynomial ring and t projecting to the element 1/∏ᵢfᵢ.
"""
function ambient_ring(D::SpecPrincipalOpen)
  if isdefined( D, :R )
    return D.R
  end
  if typeof(parent(D))<:Spec 
    R0 = ambient_ring(parent(D))
    I0 = defining_ideal(parent(D))
    R, phi, t = add_variables( R0, [(isdefined( D, :loc_var_name) ? D.loc_var_name : "t" )] ) 
    D.R = R
    D.pullbackFromParent = phi
    D.pullbackFromRoot = phi
    I = ideal( R, [ phi(g) for g in gens(I0) ] ) + ideal( R, [ one(R) - t[1]*phi(D.denom) ] )
    D.u = t[1]
    D.d = D.denom
    D.I = I
  else
    R0 = ambient_ring(root(D))
    R = ambient_ring(parent(D)) # Since we have only one localization variable, the ring stays the same.
    # We have to modify the last entry t from this list. 
    # To this end, we compute the gcd of the new denominator 
    # with the previous ones in order to throw out the obvious 
    # superfluous factors. Note, however, that this might not 
    # catch all of them, since we're not calculating this modulo 
    # the defining ideal. Doing the latter would be algorithmically 
    # too expensive, so that we decided for this tradeoff at this 
    # point.
    #
    # Say we have already localized in elements d₁,…,dᵣand we now 
    # want to also localize in f. Then we would naively do this by 
    # adding the relation 1 - t⋅d₁⋅…⋅dᵣ⋅f to the original ideal. 
    # Then the pullback in the t-variable is given by sending 
    #   t ↦ t⋅f.
    # However, the polynomial f and the dᵢmight have common factors. 
    # We look for them using the gcd( f, dᵢ) =: gᵢ. If we found 
    # any such factor, then f = fᵢ'⋅gᵢ and dᵢ = dᵢ' ⋅gᵢ. We can 
    # then use the relation 1 - t⋅d₁⋅…⋅dᵣ⋅f' and the pullback 
    # becomes t ↦ t⋅f' where f' is the product of the fᵢ's.
    # In this case, the inverse of f is represented by the element
    #   1/f = t⋅d₁'⋅…⋅dᵣ'.
    # By induction we can assume that neither of the dᵢs have 
    # pairwise common factors.
    f = copy(D.denom)
    u = last(gens(R))
    phi = parent(D).pullbackFromRoot # We can recycle that since there's no t involved.
    for d in div_denoms(parent(D))
      g = gcd(f,d)
      u *= phi(div(d,g))
      # this code could be used to also cancel higher powers. However, 
      # u needs to be adjusted for that in an inconvenient way, so 
      # I leave it in the standard way for now.
     # while true
     #   h = div(f,g)
     #   if h != zero(parent(h))
     #     f = h
     #   else
     #     break
     #   end
     # end
      f = div(f,g)
    end
    D.I = ideal( R, [ phi(g) for g in gens(root(D).I) ] ) + ideal( R, [ one(R) - last(gens(R))*phi(prod(div_denoms(parent(D)))*f) ] )
    images = copy(gens(R))
    push!( images, pop!(images)*phi(f) )
    D.pullbackFromParent = AlgebraHomomorphism( R, R, images )
    D.pullbackFromRoot = phi
    D.u = u
    D.d = f
  end
  return R
end

@Markdown.doc """
    function pullback_from_root( D::SpecPrincipalOpen )

Returns the algebra homomorphism R → R[t] of the ambient 
rings of D and its root X lifting the pullback of functions 
for the open inclusion D ↪ X.
"""
function pullback_from_root( D::SpecPrincipalOpen )
  if isdefined( D, :pullbackFromRoot )
    return D.pullbackFromRoot
  end
  ambient_ring( D ) # This also stores the homomorphism
  return D.pullbackFromRoot
end

@Markdown.doc """
    pullback_from_root( X::Spec ) = AlgebraHomomorphism( X.R, X.R, gens(X.R) )

Returns the identity of the ambient ring of X in order to be able to complete 
inductive implementations on the localization tree conveniently.
"""
pullback_from_root( X::Spec ) = AlgebraHomomorphism( X.R, X.R, gens(X.R) )

@Markdown.doc """
    pullback_from_parent( X::Spec ) = AlgebraHomomorphism( X.R, X.R, gens(X.R) )

Returns the identity of the ambient ring of X in order to be able to complete 
inductive implementations on the localization tree conveniently.
"""
pullback_from_parent( X::Spec ) = AlgebraHomomorphism( X.R, X.R, gens(X.R) )

@Markdown.doc """
    function pullback_from_parent( D::SpecPrincipalOpen )
    
Returns the ring homomorphism of the ambient rings D' = parent(D) 
and D lifting the pullback of functions for the inclusion D ↪ D'.
"""
function pullback_from_parent( D::SpecPrincipalOpen )
  if isdefined( D, :pullbackFromParent )
    return D.pullbackFromParent
  end
  ambient_ring( D ) # This also stores the homomorphism
  return D.pullbackFromParent
end

@Markdown.doc """
    function defining_ideal(D::SpecPrincipalOpen)

Returns the ideal I + ⟨1-t⋅∏ᵢfᵢ⟩ of the ambient ring R[t] of D 
where D = Spec (R/I)[(f₁)⁻¹,…,(fᵣ)⁻¹] with 
R a free polynomial ring and t projecting to the element 1/∏ᵢfᵢ.
"""
function defining_ideal(D::SpecPrincipalOpen)
  if isdefined( D, :I )
    return D.I
  end
  ambient_ring(D) # This now also prepares the ideal.
  return D.I
end

@Markdown.doc """
    function inclusion_in_parent( D::SpecPrincipalOpen ) -> AffSchMorphism

Return the inclusion D ↪ D' of a principal open subscheme D into its parent D'.
"""
function inclusion_in_parent( D::SpecPrincipalOpen )
  # return AffSchMorphism( D, parent(D), pullback_from_parent(D) )
  F = FractionField(ambient_ring(root(D)))
  return AffSchMorphism( D, parent(D), [ F(x) for x in gens(ambient_ring(root(D))) ] )
end

@Markdown.doc """
    function inclusion_in_parent( D::SpecPrincipalOpen ) -> AffSchMorphism

Return the inclusion D ↪ X of a principal open subscheme D into its root X in 
the localization tree.
"""
function inclusion_in_root( D::SpecPrincipalOpen )
  # return AffSchMorphism( D, root(D), pullback_from_root(D) )
  F = FractionField(ambient_ring(root(D)))
  return AffSchMorphism( D, root(D), [ F(x) for x in gens(ambient_ring(root(D))) ] )
end

# outer constructors

@Markdown.doc """
    function Spec( k::S, R::T ) where{S <: Ring, T <:MPolyRing} -> Spec

Returns the scheme X = Spec R with R an instance of MPolyRing over k.
"""
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

@doc Markdown.doc"""
    Spec(R::MPolyRing) -> Spec

Return the affine scheme corresponding to the ring $R/I$.
"""
function Spec( R::T, I::MPolyIdeal{U} ) where{ T <: MPolyRing, U <: MPolyElem }
  k = coefficient_ring(R)
  return Spec(k, R, I )
end

# Construct affine n-space over the ring k.
@Markdown.doc """
    function affine_space( k::Ring, n::Int; var_name::String="x" )

Return the Affine scheme 𝔸ⁿ over the base ring k. An optional argument 
`var_name` can be 
passed on to determine how the variables will be printed.
"""
function affine_space( k::Ring, n::Int; var_name::String="x" )
  R, x = PolynomialRing( k, var_name => (1:n))
  return Spec( R )
end

@Markdown.doc """
    function Base.show( io::Base.IO, X::AffineScheme )

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
    function Base.show( io::Base.IO, X::SpecPrincipalOpen )

This prints the information of X to the stream io. Whenever 
the field X.name is set, the output will be abbreviated to 
just the given name of X.
"""
function Base.show( io::Base.IO, X::SpecPrincipalOpen)
  if isdefined( X, :name )
    Base.print( io, X.name )
    return
  end
  Base.print( io, "Principal open subscheme of \n" )
  Base.print( io, root(X) )
  Base.print( io, "\n" )
  Base.print( io, "localized at\n" )
  Base.print( io,  denoms(X))
end

############################################################################
# Morphisms of affine schemes.
#
# These need to take into account both the Specs and their localizations. 
# The information on the morphism is stored either as an explicit 
# homomorphism of polynomial rings for the coordinate rings of the 
# associated specs, or as a list of rational functions, the images of 
# the variables. 

@Markdown.doc """
   mutable struct AffSchMorphism{S,Tdom, Udom, Tcod, Ucod}

A morphism of affine schemes f : X → Y. In case both X = Spec S/J and 
Y = Spec R/I are given explicitly as quotients of polynomial rings, 
this is uniquely determined by an algebra homomorphism ϕ : R → S, lifting 
the pullback f* : R/I → S/J and hence by storing the images of the 
generators xᵢ of R = k[x₁,…,xₙ]. 

For principal open subsets Spec (S/J)[g⁻¹] = U ⊂ X = Spec S/J the situation 
is similar. In this case, a morphism of schemes f : U → V = Spec (R/I)[h⁻¹] 
is completely determined by prescribing the images of the generators of R:

    xᵢ ↦ pᵢ/qᵢ= pᵢ'/gᵏ⁽ⁱ⁾∈ S[g⁻¹],

again for some lift of the pullback f* .
Hence, one can determine such a morphism either by storing these fractions 
(in the field f.imgs_frac), or the algebra homomorphism f* on explicit realizations 
of the localized rings (in the field f.pullback). 
"""
mutable struct AffSchMorphism{S,Tdom, Udom, Tcod, Ucod}
  domain::AffineScheme{S, Tdom, Udom}
  codomain::AffineScheme{S, Tcod, Ucod}

  # Variables for caching 
  #
  # The morphism must be given either in terms of an algebra homomorphism 
  # or in terms of fractions. Here we exploit the fact that to give a morphism 
  # of localized schemes, it is only necessary to prescribe the images of 
  # the generators of the ambient polynomial ring of the root.
  pullback::AlgHom
  imgs_frac::Vector{Frac{Ucod}}
  inverse::AffSchMorphism{S,Tdom,Udom,Tcod,Ucod} # optional field in case this is indeed an isomorphism
  sections::Vector{AffSchMorphism{S,Tdom,Udom,Tcod,Ucod}} # optional field for storing sections

  # optional variables used for caching
  name::String

  function AffSchMorphism( domain::AffineScheme{S,Td,Ud},
                           codomain::AffineScheme{S,Tc,Uc}, pullback::AlgHom
                           ) where {S,Td,Ud,Tc,Uc}

   if base_ring(domain) != base_ring(codomain)
      error( "the base rings of the domain and the codomain do not coincide" )
    end
    k = base_ring(domain)
    #if domain(pullback) != ambient_ring(R) || codomain(pullback) != domain.R
      #error( "the domain and codomain of the given ring homomorphism is not compatible with the affine schemes." )
    #end
    x = new{S,Td,Ud,Tc,Uc}(domain,codomain) # These first two are the only mandatory arguments
    x.pullback = pullback # This is a cached variable, so it gets it's own line
    # TODO: Implement the checks for the homomorphism to be well defined (on demand).
    return x

  end

  function AffSchMorphism( domain::AffineScheme{S,Td,Ud},
      codomain::AffineScheme{S,Tc,Uc}, imgs_frac::Vector{Frac{Ud}}
                           ) where {S,Td,Ud,Tc,Uc}
    if base_ring(domain) != base_ring(codomain)
      error( "the base rings of the domain and the codomain do not coincide" )
    end
    k = base_ring(domain)
    #if domain(pullback) != codomain.R || codomain(pullback) != domain.R
      #error( "the domain and codomain of the given ring homomorphism is not compatible with the affine schemes." )
    #end
    x = new{S,Td,Ud,Tc,Uc}()
    x.domain = domain
    x.codomain = codomain
    x.imgs_frac = imgs_frac
    # TODO: Implement the checks for the homomorphism to be well defined (on demand).
    return x
  end

  # Shortcut for the above to avoid conversion
  function AffSchMorphism( domain::AffineScheme{S,Td,Ud},
      codomain::AffineScheme{S,Tc,Uc}, 
      imgs_frac::Union{Vector{Union{Uc,Frac{Uc}}},Vector{Uc}}
    ) where {S,Td,Ud,Tc,Uc}
    R = ambient_ring(root(domain))
    Q = FractionField(R)
    new_imgs = [ Q(x) for x in imgs_frac ]
    return AffSchMorphism( domain, codomain, [ Q(x) for x in imgs_frac ] )
  end
end

function Base.show( io::Base.IO, f::AffSchMorphism )
  if isdefined( f, :name ) && isdefined( domain(f), :name ) && isdefined( codomain(f), :name )
    print( io, f.name * " : " * domain(f).name * " → " * codomain(f).name )
    return
  end

  println( io, "Morphism of affine schemes from" )
  println( io, domain(f) )
  println( io, "to" )
  println( io, codomain(f) )
  if isdefined( f, :pullback )
    println( io, "determined by the pullback of functions" ) 
    println( io, pullback(f) )
  else 
    println( io, "determined by the pullback of functions via substitution by" ) 
    println( io, f.imgs_frac )
  end
end

@Markdown.doc """
    domain(f::AffSchMorphism) -> AffineScheme

Returns the domain of definition X of f : X → Y.
"""
domain(f::AffSchMorphism) = f.domain
@Markdown.doc """
    codomain(f::AffSchMorphism) -> AffineScheme

Returns the codomain Y of f : X → Y.
"""
codomain(f::AffSchMorphism) = f.codomain

# The copy constructor
@Markdown.doc """
    function AffSchMorphism( f::AffSchMorphism ) 

Returns a copy of f.
"""
function AffSchMorphism( f::AffSchMorphism ) 
  if isdefined( f, :pullback )
    return AffSchMorphism( f.domain, f.codomain, f.pullback )
  end
  if isdefined( f, :imgs_frac ) 
    return AffSchMorphism( f.domain, f.codomain, f.imgs_frac )
  end
  error( "neither the pullback nor its fractional representation are defined" )
end

Base.copy( f::AffSchMorphism ) = AffSchMorphism( f ) 

# Construct the pullback on the level of coordinate rings 
# from the fractional representation 
@Markdown.doc """
    function pullback(h::AffSchMorphism)

Returns the pullback of functions h* associated to h : U → V.
Whenever U and V are principal open subsets of honest affine schemes, 
this involves first realizing the rings of both U and V as 
finitely generated algebras using Rabinowitsch's trick.
"""
function pullback(h::AffSchMorphism)
  if isdefined(h, :pullback)
    return h.pullback
  end
  if !isdefined(h, :imgs_frac )
    error( "Neither the fractional representation, nor the pullback is defined for this morphism" )
  end
  # Assume we have a diagram as follows 
  #    h
  #  U → V
  #  ∩   ∩
  #  X   Y
  # with both U ⊂ X and V ⊂ Y inclusions of principal 
  # open subsets. Then associated to that there is a diagram 
  # of ring homomorphisms 
  #                  h*
  #    Q' := Q[g⁻¹] ← P[f⁻¹] =: P'
  #           ↑        ↑
  #           Q        P
  # where both P = R/I and Q = S/J are quotients of 
  # polynomial rings and the sets of elements g = g₁⋅…⋅gᵣ, 
  # gⱼ ∈ Q and f = f₁⋅…⋅fₛ, fᵢ∈ P have been inverted to 
  # arrive at the open sets U and V.
  # The ring P_loc is realized as a quotient of a polynomial 
  # ring as indicated in the diagram 
  #      R[t] 
  #       ↓  ↘ π 
  #      P[t] → P', t ↦ 1/f
  # and similar for Q'
  U = domain(h)
  V = codomain(h)
  X = root(U)
  Y = root(V)
  R_t = ambient_ring(V)
  S_t = ambient_ring(U)
  R = ambient_ring(Y)
  S = ambient_ring(X)
  # Let x₁,...,xₙ be the variables of R.
  # Now the fractional representation of the morphism, i.e. 
  # the elements in imgs_frac, describes the images of these 
  # variables in Q' 
  #            pᵢ
  #   xᵢ ↦  -------          
  #            qᵢ
  # for some polynomials pᵢand qᵢin Q.
  # We need to describe the homomorphism 
  #
  #  ϕ' : R → S[t]
  # 
  # lifting the map h*. To do this, we need to bring the 
  # fractions above to a form which has a power of g in 
  # the denominator and write 
  #   pᵢ/ qᵢ = pᵢ⋅aᵢ/ gᵏ⁽ⁱ⁾ ≡ pᵢ⋅aᵢ⋅tᵏ⁽ⁱ⁾, aᵢ∈ S
  # so that the image of each variable becomes a polynomial
  # in t. According to the assumption that h(U) ⊂ V, 
  # or equivalently h⁻¹(Y∖ V) ⊂ X ∖ U, we necessarily 
  # have a relation 
  # 
  #  π(aᵢ) ⋅ π∘ϕ'(qᵢ) = π(g)ᵏ⁽ⁱ⁾
  # 
  # for some polynomials aᵢ∈ S and integers k(i), as desired. 

  n = length( h.imgs_frac )
  if n != length( gens( R ))
    error( "Number of variables in the ambient ring of the root does not coincide with the number of images provided for the homomorphism" )
  end
  images = Vector{elem_type(S_t)}()
  den_U = div_denoms(U)
  g = length(den_U)==0 ? one(S) : prod(den_U)
  Q, proj = quo( S, defining_ideal(X) )
  t = last(gens(S_t))
  # Work out the images of the variables in R first, i.e. 
  # establish ϕ'.
  for i in (1:n)
    p = numerator( h.imgs_frac[i] )
    q = denominator( h.imgs_frac[i] )
    # We need to replace the given denominators by 
    # powers of the functions which have been inverted 
    # in S to arrive at S_loc.
    (k,a) = divides_power( proj(q), proj(g) )
    if k == -1
      error( "the homomorphism is not well defined" )
    end
    push!( images, pullback_from_root(U)(p)*pullback_from_root(U)(lift(a))*t^k )
  end
  phi_prime = AlgebraHomomorphism( R, S_t, images )

  # In case Y = V we are done
  if typeof(V)<:Spec 
    h.pullback = phi_prime
    return phi_prime
  end
  
  # Now we need to extend the morphism ϕ' from above from 
  # the ring R to a map
  #   ϕ : R[t] → S[t]
  # in a way that lifts h*.
  # To this end, we need to compute the image of t. 
  # Theoretically, we find that 
  #   ϕ(t) = ϕ(1/f) = 1/ϕ'(f),
  # so we are in the same position as for the qᵢ
  # before. We compute each one of the factors fᵢby 
  # itself. This will produce a huge overhead for the 
  # numerators aᵢwhich appear, but it will hopefully 
  # be faster than performing the divides routine on big 
  # polynomials. 
  #
  # TODO: check this in practice.
  den_V = div_denoms(V)
  #f = length(den_V)==0 ? one(R) : prod(den_V)
  Q_loc, proj_loc = quo( S_t, defining_ideal(U))
  phi_of_f = one(S_t)
  for f_i in den_V
    (k,a) = divides_power( proj_loc(phi_prime(f_i)), proj_loc(pullback_from_root(U)(g)) )
    phi_of_f *= lift(a)*t^k
  end
  push!( images, phi_of_f )
  phi = AlgebraHomomorphism( R_t, S_t, images )
  h.pullback = phi
  return phi
end

# Construct the fractional representation of the morphism 
# from the explicit algebra homomorphism. 
@Markdown.doc """
    function imgs_frac(f::AffSchMorphism)

Returns the images xᵢ of the generators of the ambient ring R of the root Y 
of the codomain V of f : U → V in their fractional representation. 
"""
function imgs_frac(f::AffSchMorphism)
  # first check if this variable is already cached
  if isdefined(f, :imgs_frac)
    return f.imgs_frac
  end
  # if both forms of the morphisms are not to be found, it is not defined at all.
  if !isdefined( f, :pullback )
    error( "Neither the fractional representation, nor the pullback is defined for this morphism" )
  end

  # start reconstructing the fraction representation from 
  # the explicit ring homomorphism
  phi = pullback(f)
  # Set up the codomain of a lift of the pullback phi to the ambient ring of 
  # the domain of f. This is a localization of a *free* polynomial ring, 
  # considered as a subalgebra of the fraction field. 
  P = ambient_ring(root(domain(f)))
  F = FractionField(P)
  den = div_denoms(domain(f))
  d = length(den)
  n = length(gens(P))
  # prepare a list of the images of the generators 
  fracs = [ F(x) for x in gens(P) ]
  push!( fracs, prod( [ 1//g for g in den ] ))
  R = ambient_ring(root(codomain(f)))
  imgs_frac = elem_type(F)[]
  h = compose( pullback_from_root(codomain(f)), phi )
  for g in gens(R)
    push!(imgs_frac,evaluate(h(g), fracs))
  end
  f.imgs_frac = imgs_frac
  return imgs_frac
end

@Markdown.doc """
    function identity_morphism( X::Spec ) -> AffSchMorphism

Returns the identity morphism of a Spec.
"""
function identity_morphism( X::Spec )
  R = ambient_ring(X)
  phi = AlgebraHomomorphism( R, R, gens(R) )
  f = AffSchMorphism(X,X,phi)
  return f 
end

inclusion_in_root( X::Spec ) = identity_morphism( X )
inclusion_in_parent( X::Spec ) = identity_morphism( X )

@Markdown.doc """
    function identity_morphism( D::SpecPrincipalOpen ) -> AffSchMorphism 

Returns the identity morphism of a SpecPrincipalOpen.
"""
function identity_morphism( D::SpecPrincipalOpen )
  R = ambient_ring(D)
  phi = AlgebraHomomorphism( R, R, gens(R) )
  f = AffSchMorphism(D,D,phi)
  return f
end

#####################################################################
# Composition of affine scheme morphisms. 
# Consists basically of composition of the associated pullback 
# homomorphisms.
#
@Markdown.doc """
    function compose( f::AffSchMorphism, g::AffSchMorphism ) 

Returns the composition f ∘ g of two morphisms of affine schemes. 
NOTE: The ordering of the arguments follows the mathematical convention, 
unlike other functions for composition in Oscar! 
"""
function compose( f::AffSchMorphism, g::AffSchMorphism ) 
  if codomain(g) != domain(f) 
    error( "morphisms can not be composed" )
  end
  # TODO: Extend functionality to also work for the fractional 
  # representations of morphisms.
  phi = pullback(f)
  gamma = pullback(g)
  return AffSchMorphism( domain(g), codomain(f), compose( phi, gamma ) )
end

@Markdown.doc """
    function inverse( f::AffSchMorphism )

Computes the inverse morphism of f, assuming that f is an isomorphism 
in the first place. 
"""
function inverse( f::AffSchMorphism )
  if isdefined(f, :inverse)
    return f.inverse
  end

  X = domain(f)
  Y = codomain(f)
  R = ambient_ring(X)
  I = defining_ideal(X)
  x_vars = gens(R)
  S = ambient_ring(Y)
  J = defining_ideal(Y)
  y_vars = gens(S)
  phi = pullback(f)
  RmodI, proj_R = quo(R, I)
  SmodJ, proj_S = quo(S, J)
  phi_bar = AlgebraHomomorphism(SmodJ, RmodI, [proj_R(x) for x in phi.image])
  @assert isinjective(phi_bar) "the associated algebra homomorphism is not injective"
  @assert issurjective(phi_bar) "the associated algebra homomorphism is not surjective"
  x_preim = [ preimage(phi_bar, x)[1] for x in gens(RmodI) ]
  phi_inv = AlgebraHomomorphism( R, S, [lift(f) for f in x_preim])
  f_inv = AffSchMorphism(Y, X, phi_inv)
  f.inverse = f_inv
  f_inv.inverse = f
  return f_inv
end

@Markdown.doc """
    mutable struct Glueing

Glueing of two affine schemes X and Y along principal open subsets 
U ⊂ X and V ⊂ Y via an isomorphism Φ : U → V. 
A glueing consists of a diagram X ↩ U → V ↪ Y with the outer maps 
the open embeddings of principal open subsets. 
"""
mutable struct Glueing
  firstPatch::AffineScheme
  secondPatch::AffineScheme
  inclusionFirstPatch::AffSchMorphism
  inclusionSecondPatch::AffSchMorphism
  glueingIsomorphism::AffSchMorphism

  # field used for caching
  inverse::Glueing

  function Glueing( 
      inclusionFirstPatch, 
      inclusionSecondPatch, 
      glueingIsomorphism 
    )
    if domain( inclusionFirstPatch ) != domain( glueingIsomorphism ) 
      error( "Morphisms can not be composed for glueing" )
    end
    if codomain( glueingIsomorphism ) != domain( inclusionSecondPatch ) 
      error( "Morphisms can not be composed for glueing" )
    end
    G = new()
    G.firstPatch = codomain( inclusionFirstPatch )
    G.secondPatch = codomain( inclusionSecondPatch )
    G.inclusionFirstPatch = inclusionFirstPatch
    G.inclusionSecondPatch = inclusionSecondPatch
    G.glueingIsomorphism = glueingIsomorphism
    return G
  end
end

# Glueing for two principal open subsets of a common root scheme X
# along their intersection in X
@Markdown.doc """
    function Glueing( U::SpecPrincipalOpen, V::SpecPrincipalOpen )

Provides the glueing between two principal open subsets U and V of 
some common parent X in their localization trees. If no such common 
parent is found, this necessarily returns an error. Otherwise, the 
intersection U ∩ V ⊂ X is computed and the glueing returned. 
"""
function Glueing( U::SpecPrincipalOpen, V::SpecPrincipalOpen )
  root(U) == root(V) || error( "schemes do not belong to the same root" )
  X = root(U)
  R = ambient_ring(X)
  den_U = denoms(U)
  den_V = denoms(V)
  f = prod(den_U)
  g = prod(den_V)
  U_loc = localize(U,g)
  V_loc = localize(V,f)
  iso = AffSchMorphism( U_loc, V_loc, gens(R) )
  return Glueing( inclusion_in_parent(U_loc), inclusion_in_parent(V_loc), iso )
end

@Markdown.doc """
    first_patch( G::Glueing ) = G.firstPatch

Returns the first patch for the glueing G.
"""
first_patch( G::Glueing ) = G.firstPatch
@Markdown.doc """
    second_patch( G::Glueing ) = G.secondPatch

Returns the second patch for the glueing G.
"""
second_patch( G::Glueing ) = G.secondPatch
@Markdown.doc """
    overlap(G::Glueing) = domain(G.inclusionFirstPatch)

Returns the overlap U ∩ V of the two patches U and V of the glueing G.
Note, that one naturally has two realizations of this overlap as 
the Spec of affine algebras, one as a localization of U and another one 
as a localization of V. By convention, it is always the first realization 
of the overlap that is being returned. 
"""
overlap(G::Glueing) = domain(G.inclusionFirstPatch)
@Markdown.doc """
    first_inclusion( G::Glueing ) = G.inclusionFirstPatch

Returns the inclusion of the overlap U ∩ V ↪ U into the first patch of the glueing G.
"""
first_inclusion( G::Glueing ) = G.inclusionFirstPatch
@Markdown.doc """
    second_inclusion( G::Glueing ) = G.inclusionSecondPatch

Returns the inclusion of the overlap U ∩ V ↪ V into the second patch of the glueing G.
"""
second_inclusion( G::Glueing ) = G.inclusionSecondPatch
@Markdown.doc """
    glueing_isomorphism( G::Glueing ) = G.glueingIsomorphism

Returns the glueing isomorphism for the two realizations of the overlap U ∩ V 
of the two patches. These two realizations come from considering U ∩ V as a 
localization of both either U and of V.
"""
glueing_isomorphism( G::Glueing ) = G.glueingIsomorphism

@Markdown.doc """
    function inverse(G::Glueing)

Provides the 'inverse glueing' meaning that the the inverse of the 
glueing map on the overlap of the two patches is used. 
The result is cached in the inverse field.
"""
function inverse(G::Glueing)
  if isdefined(G, :inverse) 
    return G.inverse
  end

  G_inv = Glueing( second_inclusion(G), first_inclusion(G), inverse(glueing_isomorphism(G)))
  G.inverse = G_inv
  G_inv.inverse = G
  return G_inv
end

@Markdown.doc """
    function restict(G::Glueing, U::AffineScheme, V::AffineScheme)

For a glueing G : X ↩ X ∩ Y ↪ Y and two principal open subsets U ⊂ X and V ⊂ Y 
this computes the induced glueing U ↩ U ∩ V ↪ V of U and V. 
"""
function restict(G::Glueing, U::AffineScheme, V::AffineScheme)
  X = first_patch(G)
  Y = first_patch(G)
  is_parent_of(X, U) || error( "The first argument is not a principal open subset of the first patch of the glueing" )
  is_parent_of(Y, V) || error( "The second argument is not a principal open subset of the second patch of the glueing" )
  i = first_inclusion(G)
  j = second_inclusion(G)
  phi = glueing_isomorphism(G)
  denoms_U = denoms_up_to(U, X)
  denoms_V = denoms_up_to(V, Y) 
  phi_inv = inverse(phi)
  denoms_V_in_U = Vector{elem_type(ambient_ring(root(U)))}()
  
  U_loc = localize( U, [pullback(phi)(d) for d in denoms_V])
  V_loc = localize( V, [pullback(phi_inv)(d) for d in denoms_U])
  return Glueing( restrict(inclusion_in_root(U_loc), U_loc, U), 
		 restrict(inclusion_in_root(V_loc), V_loc, V), 
		 restriction( phi, U_loc, V_loc ) )
end

@Markdown.doc """
    function compose(G::Glueing, H::Glueing)

For two glueings G : X ↩ X ∩ Y ↪ Y and H : Y ↩ Y ∩ Z ↪ Z this computes the 
induced glueing G ∘ H : X ∩ Y ↩ X ∩ Y ∩ Z ↪  Y ∩ Z.
"""
function compose(G::Glueing, H::Glueing)
  X = first_patch(G)
  Y = second_patch(G)
  Y == first_patch(H) || error( "glueings can not be composed." )
  Z = second_patch(H)
  XnY = domain(first_inclusion(G))
  YnZ = domain(first_inclusion(H))
  error( "not implemented, yet" )
end



@Markdown.doc """
    function Base.show( io::Base.IO, G::Glueing )

Print the glueing G to io. If the field G.name is set, then the 
print will be abbreviated to the output of this name.
"""
function Base.show( io::Base.IO, G::Glueing )
  if isdefined( G, :name )
    Base.print( io, G.name )
    return
  end
  Base.println( io, "Glueing of" )
  Base.println( io, G.firstPatch )
  Base.println( io, "and" )
  Base.println( io, G.secondPatch )
  Base.println( io, "along" )
  Base.println( io, domain(G.inclusionFirstPatch) )
  return
end

function Base.copy( G::Glueing )
  return Glueing( G.inclusionFirstPatch, G.inclusionSecondPatch, G.glueingIsomorphism )
end

################################################################################
# Covering of a scheme.
# It consists of a list of affine patches together with their glueings.
# The glueings are stored in an array (which is usually half empty, since 
# glueings need only be provided for ordered pairs of indices). 
# The index methods for double arrays have been overloaded so that 
# OrderedMultiindex can be used to address them.
#
@Markdown.doc """
    mutable struct Covering

Structure for storing a covering of a CoveredScheme X by affine 
patches Uᵢ, i = 1,…,N. The necessary data consists of these affine 
patches and the information on how to glue them pairwise via 
glueins Gᵢⱼ: Uᵢ↩ Uᵢ∩ Uⱼ ↪ Uⱼ. 

Once created, a covering must not be altered, since in general several 
other objects, like e.g. vector bundles and morphisms of schemes, will 
refer to such coverings and they need to rely on them to not change. 
If one needs another covering derived from the given one, then a 
new covering needs to be set up and stored. 
"""
mutable struct Covering
  patches::Vector{AffineScheme}	# A list of open patches for this scheme
  glueings::Matrix{Union{Glueing,Nothing}}     # A list of glueings between the patches 

  # Variables for caching
  name::String
  						# listed above
  function Covering( patches::Vector{AffineScheme}, glueings::Matrix{Union{Glueing,Nothing}} )
    return new( patches, glueings )
  end

  function Covering()
    C = new()
    C.patches = AffineScheme[]
    C.glueings = Matrix{Union{Nothing,Glueing}}(nothing,0,0)
    #C.glueings = Matrix{Union{Nothing,Glueing}}
    return C
  end

  function Covering( X::AffineScheme ) 
    C = new()
    C.patches = [X]
    C.glueings = Array{Nothing,2}(nothing,1,1)  # no glueing on the diagonal.
    return C
  end
end

function Base.copy( C::Covering )
  # I'm afraid we need a copy here. Otherwise an alteration 
  # of the original covering will also affect the new one.
  return Covering( copy(C.patches), copy(C.glueings) )
end

@Markdown.doc """
    function patches( C::Covering ) -> Vector{AffineScheme}

Returns the patches of the covering.
"""
function patches( C::Covering )
  if isdefined( C, :patches )
    return C.patches
  else 
    return nothing 
  end 
end

@Markdown.doc """
    function patch( C::Covering, i::Int )
    
Returns the i-th patch in the list of patches for the covering C.
"""
function patch( C::Covering, i::Int )
  if isdefined( C, :patches )
    return( C.patches[i] )
  else
    return nothing
  end
end

@Markdown.doc """
    function glueings( C::Covering ) -> Matrix{Union{Glueing,Nothing}}

Returns the matrix of glueings for the given covering. Note that 
some entries (such as the diagonal ones) do not need to be assigned. 
"""
function glueings( C::Covering )
  if isdefined( C, :glueings )
    return C.glueings
  else
    return nothing
  end
end

function Base.show( io::Base.IO, C::Covering )
  if isdefined( C, :name )
    Base.print( io, C.name )
    return
  end
  # TODO: Implement another way to print it.
end

@Markdown.doc """
    function _add_patch!( C::Covering, X::AffineScheme, glueings::Vector{Glueing} )
    
**WARNING: This is an internal function. The user is not supposed to use this since 
it changes the covering.**

Add a patch X to the given covering C. To do this, one needs to not only provide 
the new patch, but also how it glues to those patches Uᵢin C that are already there. 
Thus one needs to also provide a vector (G₁,…,Gₙ) of glueings Gᵢ: Uᵢ↩ Uᵢ∩ X ↪ X.
"""
function _add_patch!( C::Covering, X::AffineScheme, glueings::Vector{Glueing} )
  push!( C.patches, X )
  n = length(glueings)
  if n != length(patches(C))-1
    error( "the number of glueings does not coincide with the number of patches so far" )
  end
  for i in (1:n)
    if first_patch( glueings[i] ) != patches(C)[i] 
      error( "Domain of the glueing does not coincide with the patch" )
    end
    if second_patch( glueings[i] ) != X
      error( "Codomain of the glueing does not coincide with the new patch" )
    end
  end
  C.glueings = vcat( hcat( C.glueings, glueings ), Array{Nothing,2}(nothing,1,n+1) )
  return C
end

@Markdown.doc """
    function intersect( X::AffineScheme, Y::AffineScheme, C::Covering )

Intersection of two patches X and Y in a given covering.
Assumes that both X and Y are patches of the given covering.
Returns the intersection as an affine scheme 
together with the two open inclusions.
"""
function intersect( X::AffineScheme, Y::AffineScheme, C::Covering )
  if !(X in patches(C)) 
    error( "X is not a patch in the covering C" )
  end
  if !(Y in patches(C)) 
    error( "Y is not a patch in the covering C" )
  end
  if X == Y 
    return ( X, identity_morphism(X), identity_morphism(X) )
  end
  i = indexin( X, patches(C) )[1]
  j = indexin( Y, patches(C) )[1]
  # Keep in mind that the glueing is realized 
  # in a directed way. Thus, we need to make a distinction.
  if i < j
    G = glueings(C)[i,j]
    return (overlap(G),first_inclusion(G),second_inclusion(G))
  else
    G = glueings(C)[j,i]
    return (overlap(G),second_inclusion(G),first_inclusion(G))
  end
end

@Markdown.doc """
    function intersect( U::AffineScheme, V::AffineScheme )

Returns the intersection of two principal open sets U and V in some common affine 
ambient scheme X in their localization trees.
"""
function intersect( U::AffineScheme, V::AffineScheme )
  root(U) == root(V) || error( "schemes do not have the same root" )
  R = ambient_ring(root(U))
  den_U = denoms(U)
  den_V = denoms(V)
  # A missing implementation for rings forces us to do a weird workaround here:
  f = length(den_U) == 0 ? one(R) : prod(den_U)
  g = length(den_V) == 0 ? one(R) : prod(den_V)
  U_loc = localize(U,g)
  V_loc = localize(V,f)
  iso = AffSchMorphism( U_loc, V_loc, gens(R) )
  return (U_loc, inclusion_in_parent(U_loc), compose( inclusion_in_parent(V_loc), iso ))
end

@Markdown.doc """
    function union_of_schemes( U::AffineScheme, V::AffineScheme ) -> CoveredScheme

Returns the union of two schemes U and V as a CoveredScheme W = U ∪ V, 
provided that both U and V have a common parent X in their localization tree. 
"""
function union_of_schemes( U::AffineScheme, V::AffineScheme )
  root(U) == root(V) || error( "schemes do not have the same root" )
  C = Covering(U)
  _add_patch!( C, V, [Glueing(U,V)] )
  return CoveredScheme(C)
end

@Markdown.doc """
    function union_of_schemes( C::Covering, V::AffineScheme ) -> CoveredScheme

Returns the union of a subscheme U ↪ X of some affine scheme X = Spec R/I 
provided by means of an explicit covering C of U with some further subscheme 
V ⊂ X of the same parent in the localization tree.
"""
function union_of_schemes( C::Covering, V::AffineScheme )
  U = patches(C)
  l = length(U)
  G = Glueing[] # empty vector for additional glueings
  for i in (1:l) 
    # Compute the glueing of U[i] with V 
    root(U[i]) == root(V) || error( "the $i-th patch of the covering does not have the same root as the new patch" )
    push!( G, Glueing( U[i], V ))
  end
  D = copy(C)
  _add_patch!(D,V,G)
  return CoveredScheme( D )
end

@Markdown.doc """
    function restrict( f::AffSchMorphism, U::SpecPrincipalOpen, V::SpecPrincipalOpen )

Return the restriction of a morphism of affine schemes f : X → Y to principal 
open subsets U ⊂ X and V ⊂ Y. 

**Caution: No plausibility check as to whether or not f(U) ⊂ V is being performed!**
"""
function restrict( f::AffSchMorphism, U::SpecPrincipalOpen, V::AffineScheme )
  i = inclusion_in_root(U)
  j = inclusion_in_root(V)
  domain(f) == codomain(i) || error( "the given map is not compatible with the open subsets" )
  codomain(f) == codomain(j) || error( "the given map is not compatible with the open subsets" )
  g = compose( f, i )
  f_res = AffSchMorphism( U, V, imgs_frac(g) )
  if isdefined( f, :name ) 
    f_res.name = f.name
  end
  return f_res
end

restrict( f::AffSchMorphism, X::Spec, V::Spec ) = f

function restrict( f::AffSchMorphism, X::Spec, V::SpecPrincipalOpen ) 
  f_res = compose( f, inclusion_in_root(V) )
  if isdefined( f, :name ) 
    f_res.name = f.name
  end
  return f_res
end


@Markdown.doc """
    function union_of_schemes( X::AffineScheme, Y::AffineScheme, C::Covering )
    
Provide the intersection of two affine patches X and Y in a covering C. This 
assumes that both X and Y can be found in the list of patches of C.

**Warning: This has not yet been implemented!**
"""
function union_of_schemes( X::AffineScheme, Y::AffineScheme, C::Covering )
  error( "Not implemented yet." )
  # TODO: This is not functional yet. We first need a struct 
  # for morphisms of CoveredScheme.
end

@Markdown.doc """
    mutable struct Refinement
    
Structure for a refinement of a covering. 

Say C = { Uᵢ}ᵢ₌₁ᴹ and D = { Vⱼ}ⱼ₌₁ᴺ are two coverings of a CoveredScheme X. Then D is a refinement 
of C if there exists a map of the index ranges 

  n : {1,…,N} → {1,…,M}, j ↦  i = n(j)

such that there is an inclusion ρⱼ: Vⱼ↪ Uₙ₍ⱼ₎. 
This is the datum being stored. 

**WARNING: This has not been fully implemented, yet!**
"""
mutable struct Refinement
  original_cover::Covering
  refined_cover::Covering
  index_map::Vector{Int}
  inclusions::Vector{AffSchMorphism}

  # TODO: Implement the constructors.
end

@Markdown.doc """
    function original_cover( ρ::Refinement ) -> Covering

Returns the original cover of the refinement.
"""
function original_cover( ρ::Refinement )
  if isdefined( ρ, :original_cover )
    return ρ.original_cover
  end
  return nothing
end

@Markdown.doc """
    function refined_cover( ρ::Refinement ) -> Covering

Returns the refined cover of the refinement.
"""
function refined_cover( ρ::Refinement )
  if isdefined( ρ, :refined_cover )
    return ρ.refined_cover
  end
  return nothing
end

@Markdown.doc """
    function map_index( ρ::Refinement, i::Int ) -> Int

The map 

  n : {1,…,N} → {1,…,M}, j ↦  i = n(j)

of the index sets for the original covering C = { Uᵢ}ᵢ₌₁ᴹ and its 
refinement D = { Vⱼ}ⱼ₌₁ᴺ such that ρⱼ: Vⱼ↪ Uₙ₍ⱼ₎ are the open 
inclusions for the refinements.
"""
function map_index( ρ::Refinement, i::Int )
  return ρ.index_map[i]
end

@Markdown.doc """
    function inclusions( ρ::Refinement ) -> Vector{AffSchMorphism}

Returns the vector of inclusion morphisms ρⱼ: Vⱼ↪ Uₙ₍ⱼ₎ 
for the refined covering D = { Vⱼ}ⱼ₌₁ᴺ of C = { Uᵢ}ᵢ₌₁ᴹ.
"""
function inclusions( ρ::Refinement )
  if isdefined( ρ, :inclusions ) 
    return ρ.inclusions
  end
  return nothing
end

@Markdown.doc """
    mutable struct CoveredScheme

Models a covered scheme X = ⋃ᵢ₌₁ᴺ Uᵢ obtained by glueing of affine 
schemes Uᵢ.

This stores a list of coverings for X where at least one original covering needs 
to be given. Within the life span of X there will most likely be the need to 
also have refinements of this covering (e.g. for the description of vector bundles 
that can not be trivialized on the original patches), so more coverings can be added. 
To also maintain the information of their interplay, there is a second field to store 
a matrix of refinements for these coverings.
"""
mutable struct CoveredScheme
  coverings::Vector{Covering}
  refinements::Array{Refinement,2}

  # fields for caching
  name::String

  function CoveredScheme()
    return CoveredScheme( Covering() )
  end

  function CoveredScheme( C::Covering )
    X = new()
    X.coverings = [C]
    #=
    if length( patches(C) ) == 0 
      error( "empty covering. Cannot determine base ring." )
      return
    end
    X.base_ring = base_ring( patches(C)[1] )
    =#
    return X
  end

end

# The copy constructor
function CoveredScheme( X::CoveredScheme )
  # Only copy the default covering for now. Maybe 
  # we also want to copy the other cached data.
  return CoveredScheme( copy(coverings(X)[1]) )
end

function Base.copy( X::CoveredScheme )
  return CoveredScheme( X ) 
end

@Markdown.doc """
    function CoveredScheme( U::AffineScheme ) -> CoveredScheme

Turn the affine scheme U into a CoveredScheme with a one-patch-cover.
"""
function CoveredScheme( U::AffineScheme )
  return CoveredScheme(Covering(U))
end

@Markdown.doc """
    function coverings( X::CoveredScheme )

Return the available coverings of X.
"""
function coverings( X::CoveredScheme )
  return X.coverings
end

@Markdown.doc """
    function set_name!( X::CoveredScheme, name::String )

Give the CoveredScheme X a name that will then be used for printing.
"""
function set_name!( X::CoveredScheme, name::String )
  X.name = name
  return name
end

@Markdown.doc """
    function name( X::CoveredScheme ) 

Return the given name of X whenever the field X.name is set.
"""
function name( X::CoveredScheme ) 
  if isdefined( X, :name )
    return X.name
  end
  return nothing
end

@Markdown.doc """
    function Base.show( io::Base.IO, X::CoveredScheme )

Print X to the stream io. Only the given name of X will be 
used once the field X.name is being set. 
"""
function Base.show( io::Base.IO, X::CoveredScheme )
  if isdefined( X, :name )
    Base.print( io, X.name )
    return
  end
  Base.print( io, "Covered Scheme with $(length(X.coverings)) covering." )
end

@Markdown.doc """
    mutable struct CovSchMorphism

A morphism f : X → Y of CoveredSchemes.
This needs to be described in affine coverings X = ⋃ᵢ₌₁ᴹ Uᵢand Y = ⋃ⱼ₌₁ᴺ Vⱼ
which have to be chosen in such a way that for every j ∈ {1,…,N} 
the preimage f⁻¹(Vⱼ) is a *union* of patches Uᵢ for some subset 
N(j) ⊂ {1,…,M}. Then for every i ∈ N(j) the restrictions 
f : Uᵢ→ Vⱼ can be stored as an instance of AffSchMorphism, 
subject to the obvious compatibility requirements.
"""
mutable struct CovSchMorphism
  domain::CoveredScheme
  codomain::CoveredScheme
  dom_covering::Covering
  cod_covering::Covering
  N::Vector{Set} 
  restrictions::Vector{Vector{AffSchMorphism}}
end

@Markdown.doc """
    function projective_space( k::Ring, n::Int )

Construct projective space ℙₖⁿ over the ring k as a covered 
scheme with its standard open cover.
"""
function projective_space( k::Ring, n::Int; var_name="s" )
  # Start the standard covering with the first patch for x_0 ≠ 0.
  # println( "assembling $n-dimensional projective space over $k." )
  X0 = affine_space(k,n, var_name=var_name)
  set_name!( X0, "𝔸^$(n)_0" )
  C = Covering(X0)
  # Add the other patches successively
  for j in (1:n) 
    # println("Variable in the outer loop j = $j")
    # The patch for x_j ≠ 0
    # println("Glueing in the following new patch:")
    Y = affine_space( k, n, var_name=var_name )
    set_name!(Y,"𝔸^$(n)_$(j)")
    # come up with the glueings
    G = Vector{Glueing}()
    # println("New glueings for this step so far:") 
    for i in (0:j-1)
      # println("Variable of the inner loop i = $i")
      # The patch for x_i ≠ 0 with i < j
      X = patches(C)[i+1]
      # println( "first patch" )
      R = ambient_ring( Y ) 
      V = localize( Y, gens(R)[i+1] )
      S = ambient_ring( X ) 
      U = localize( X, gens(S)[j] )
      R_loc = ambient_ring( V )
      S_loc = ambient_ring( U )
      R_vars = root_variables(V)
      S_vars = root_variables(U)
      R_unit = inverted_element(V)
      S_unit = inverted_element(U)
      Phi = AlgebraHomomorphism( R_loc, S_loc, 
			      vcat( [ S_vars[k]*S_unit for k in (1:i) ], # the first i-1 variables have the same indices
				   [S_unit], # this is inverted to pass to the other patch
				   [ S_vars[k-1]*S_unit for k in (i+2:j) ], # fill up with intermediate variables
				   [ S_vars[k]*S_unit for k in (j+1:n) ], # skip the j-th variable
				   [ S_vars[j] ] ) # map the unit to the appropriate coordinate
			      )

      f = AffSchMorphism( U, V, Phi )
      glueing = Glueing( inclusion_in_parent(U), inclusion_in_parent(V), f )
      push!( G, glueing )
    end
    _add_patch!( C, Y, G )
  end
  C.name = "standard cover"
  IP = CoveredScheme( C )
  set_name!( IP, "ℙ^$n" )
  return IP
end

##########################################################
# 
# A sheaf of ideals on a CoveredScheme X
#
# Just as every coherent sheaf, a sheaf of ideals 
# needs to refer to an explicit covering C of X which 
# has to be listed in the stored coverings for X. 
# The second datum necessary is a list of ideals 
# I_k ⊂ R_k for all affine schemes Spec R_k in C. 
# Note that, for simplicity, we simply store the 
# preimages of the ideals in the chosen ambient 
# rings of (the roots of) the affine patches. 
#
# We refrain from implementing this explicitly for 
# affine schemes since there ideal sheaves are just 
# in 1:1 correspondence with actual ideals.
#
@Markdown.doc """
    mutable struct IdealSheaf <: AbstractCoherentSheaf 

A sheaf of ideals on a CoveredScheme X

Just as every coherent sheaf, a sheaf of ideals 
needs to refer to an explicit covering C of X which 
has to be listed in the stored coverings for X. 
The second datum necessary is a list of ideals 
Iₖ ⊂ Rₖ for all affine schemes Uₖ = Spec Rₖ in C. 
Note that, for simplicity, we simply store the 
preimages of the ideals in the chosen ambient 
rings of (the roots of) the affine patches. 

We refrain from implementing this explicitly for 
affine schemes since there ideal sheaves are just 
in 1:1 correspondence with actual ideals.
"""
mutable struct IdealSheaf <: AbstractCoherentSheaf 
  parent::CoveredScheme
  covering::Covering
  ideals::Vector{MPolyIdeal}

  # optional variables used for caching
  name::String

  function IdealSheaf( X::CoveredScheme, C::Covering, I::Vector{MPolyIdeal} )
    if !( C in coverings(X) ) 
      error( "covering does not belong to the variety" )
    end
    n = length(patches(C))
    if n != length(I) 
      error( "number of ideals provided does not coincide with the number of patches in the covering" )
    end
    for i in (1:n)
      if base_ring(I[i]) != ambient_ring(root(patch(C,i)))
	error( "ideal number $i does not live in the base ring of the $i-th patch" )
      end
    end
    return new( X, C, I )
  end
end

@Markdown.doc """
    function IdealSheaf( I::IdealSheaf )

The copy constructor. Returns a copy of I.
"""
function IdealSheaf( I::IdealSheaf )
  return IdealSheaf( I.parent, I.covering, copy(I.ideals) )
end

@Markdown.doc """
    parent( I::IdealSheaf ) -> CoveredScheme

Returns the CoveredScheme X on which the IdealSheaf I is defined.
"""
parent( I::IdealSheaf ) = I.parent
@Markdown.doc """
    covering( I::IdealSheaf ) -> Covering

Returns the covering C of the CoveredSchem X with respect to which the 
ideal sheaf I is defined. 
"""
covering( I::IdealSheaf ) = I.covering
@Markdown.doc """
    ideals( I::IdealSheaf ) -> Vector{MPolyIdeal}

Return the vector of ideals representing I locally on the patches. 
"""
ideals( I::IdealSheaf ) = I.ideals

@Markdown.doc """
    function to_ideal_sheaf( I::MPolyIdeal )

Take a homogeneous ideal I in a graded homogeneous 
ring S and turn it into an ideal sheaf on the 
corresponding projective space.

**WARNING: This has not been implemented, yet!**
"""
function to_ideal_sheaf( I::MPolyIdeal )
  error( "not implemented yet" )
  # TODO: Implement this based on graded modules in 
  # Oscar. The routine is supposed to recognize automatically 
  # from the parent how to construct the associated projective 
  # space and return the ideal sheaf on it. 
end

@Markdown.doc """
    mutable struct CoherentSheaf <: AbstractCoherentSheaf 

A sheaf of modules on a CoveredScheme X

Everything is similar to the IdealSheaf,
only that instead of ideals we store modules over 
the chosen ambient rings of the patches.
"""
mutable struct CoherentSheaf <: AbstractCoherentSheaf 
  parent::CoveredScheme
  covering::Covering
  modules::Vector{FPModule{MPolyElem}}

  # optional variables used for caching
  name::String
  
  # TODO: Implement the constructors.
end

@Markdown.doc """
    parent( I::CoherentSheaf ) -> CoveredScheme

Returns the CoveredScheme on which the CoherentSheaf I is defined.
"""
parent( I::CoherentSheaf ) = I.parent
@Markdown.doc """
    covering( I::CoherentSheaf ) -> Covering

Return the covering C of the CoveredScheme X on which the CoherentSheaf
I is defined.
"""
covering( I::CoherentSheaf ) = I.covering
@Markdown.doc """
    modules( I::CoherentSheaf ) -> Vector{FPModule{MPolyElem}}

Returns the vector of finitely generated modules on the 
patches of the CoveredScheme X on which I is defined. 
"""
modules( I::CoherentSheaf ) = I.modules


@Markdown.doc """
    mutable struct LineBundle <: AbstractCoherentSheaf

A line bundle on a CoveredScheme X, given by its 
local trivializations on the patches of a specific Covering C.
"""
mutable struct LineBundle <: AbstractCoherentSheaf
  parent::CoveredScheme
  covering::Covering
  transitions::Array{MPolyElem,2}

  # optional variables used for caching
  name::String

  function LineBundle( parent::CoveredScheme, covering::Covering, transitions::Array{MPolyElem,2} )
    # TODO: Implement cheap plausibility checks.
    return new( parent, covering, transitions )
  end
end

@Markdown.doc """
    function LineBundle( L::LineBundle )

Returns a copy of L.
"""
function LineBundle( L::LineBundle )
  return LineBundle( L.parent, L.covering, copy(L.transitions) )
end

function Base.copy( L::LineBundle )
  return LineBundle( L )
end

@Markdown.doc """
    parent( L::LineBundle ) -> CoveredScheme

Returns the CoveredScheme X on which L is defined.
"""
parent( L::LineBundle ) = L.parent
@Markdown.doc """
    covering( L::LineBundle ) -> Covering

Returns the Covering C of the CoveredScheme X on which L is defined.
"""
covering( L::LineBundle ) = L.covering
@Markdown.doc """
    transitions( L::LineBundle ) -> Matrix{MPolyElem}

Returns the matrix of transition functions for L for 
pairs of affine patches in the Covering C of the 
CoveredScheme X on which L is defined. 
"""
transitions( L::LineBundle ) = L.transitions

@Markdown.doc """
    function OO( k::Ring, n::Int, d::Int ) -> LineBundle

Returns Serre's twisting sheaf 𝒪(d) on ℙₖⁿ, trivialized in 
the canonical covering.
"""
function OO( k::Ring, n::Int, d::Int )
  IPn = projective_space( k, n )
  C = coverings( IPn )[1]
  transitions = Array{MPolyElem,2}(undef,n+1,n+1)
  if d >= 0
    for i in OrderedMultiindex(n+1,2)
      a = inverted_element(overlap(glueings(C)[getindex(i,1),getindex(i,2)]))
      transitions[getindex(i,1), getindex(i,2)] = a^d
      #transitions[getindex(i,1), getindex(i,2)] = (inverted_element(overlap(glueings[i])))^d
      transitions[i[2],i[1]] = (denom(overlap(glueings(C)[i[1],i[2]])))^d
    end
  else
    d = -d
    for i in OrderedMultiindex(n+1,2)
      a = inverted_element(overlap(glueings(C)[getindex(i,1),getindex(i,2)]))
      transitions[getindex(i,2), getindex(i,1)] = a^d
      #transitions[getindex(i,1), getindex(i,2)] = (inverted_element(overlap(glueings(C)[i])))^d
      transitions[i[1],i[2]] = (denom(overlap(glueings(C)[getindex(i,1),getindex(i,2)])))^d
    end
  end
  L = LineBundle( IPn, C, transitions )
  L.name = "OO($n,$d)"
  return L
end

@Markdown.doc """
    OO( n::Int, d::Int ) = OO( QQ, n, d )
"""
OO( n::Int, d::Int ) = OO( QQ, n, d )

@Markdown.doc """
    function Base.show( io::Base.IO, L::LineBundle )

Print the LineBundle L to the stream io. Whenever the 
field L.name is set, this output is abbreviated to printing 
that name. 
"""
function Base.show( io::Base.IO, L::LineBundle )
  if isdefined( L, :name )
    Base.print( io, L.name )
    return
  end
  Base.print( io, "Line bundle on " )
  Base.print( io, L.parent )
  #Base.print( io, " with respect to the covering " )
  #Base.print( io, L.covering )
  Base.print( io, " given by the transition functions\n" )
  n = length( patches(L.covering) )
  for i in (1:n)
    for j in (1:n)
      if i == j 
	Base.println( io, "($i,$j): 1" )
      else 
	Base.println( io, "($i,$j): $(L.transitions[i,j])" )
      end
    end
  end
end

@Markdown.doc """
    mutable struct VectorBundle <: AbstractCoherentSheaf

An algebraic vector bundle on a CoveredScheme X trivialized 
on a specific Covering C.
"""
mutable struct VectorBundle <: AbstractCoherentSheaf
  rank::Int
  parent::CoveredScheme
  covering::Covering
  transitions::Array{MatSpaceElem,2}
  
  # optional variables used for caching
  name::String

  function VectorBundle( 
      rank::Int,
      parent::CoveredScheme, 
      covering::Covering, 
      transitions::Matrix{MatSpaceElem} 
    )
    E = new()
    E.rank = rank
    E.parent = parent
    E.covering = covering
    E.transitions = transitions
    # TODO: Implement the cheap plausibility checks.
    return E
  end
end

@Markdown.doc """
    function VectorBundle( E::VectorBundle ) -> VectorBundle

Returns a copy of E.
"""
function VectorBundle( E::VectorBundle )
  return VectorBundle( E.rank, E.parent, E.covering, copy(E.transitions) )
end

function Base.copy( E::VectorBundle )
  return VectorBundle( E )
end

@Markdown.doc """
    function rank( E::VectorBundle ) -> Int

Returns the rank of the VectorBundle E.
"""
function rank( E::VectorBundle )
  return E.rank
end

@Markdown.doc """
    parent( E::VectorBundle ) -> Covering

Returns the CoveredScheme X on which E is defined.
"""
parent( E::VectorBundle ) = E.parent
@Markdown.doc """
    covering( E::VectorBundle ) -> Covering

Returns the Covering C of the CoveredScheme X on which E is defined.
"""
covering( E::VectorBundle ) = E.covering

@Markdown.doc """
transitions( E::VectorBundle ) -> Matrix{MatSpaceElem}

Returns the matrix of transition matrices for E.
"""
transitions( E::VectorBundle ) = E.transitions


@Markdown.doc """
    mutable struct HomVecBundle <: AbstractCoherentSheaf

Model for the vector bundle Hom(E,F) for two vector bundles E and F on the 
same scheme X. 
This deserves a struct of its own because of the special functionality. 
"""
mutable struct HomVecBundle <: AbstractCoherentSheaf
  parent::CoveredScheme
  covering::Covering
  domain::VectorBundle
  codomain::VectorBundle

  function HomVecBundle( E::VectorBundle, F::VectorBundle )
    parent(E) == parent( F ) || error( "The two input vector bundles are not defined on the same scheme" )
    covering(E) == covering( F ) || error( "Constructor is not yet implemented for vector bundles which are defined on different coverings" )
    return new(parent(E),covering(E),E,F)
  end
end

@Markdown.doc """
    function VectorBundle( H::HomVecBundle )

Turns the HomVecBundle H into an honest vector bundle, forgetting 
about its hom structure. 
"""
function VectorBundle( H::HomVecBundle )
  error( "not implemented, yet" )
end

@Markdown.doc """
    mutable struct HomBundleSection

Section in a bundle Hom(E,F) for vector bundles E and F over a common base scheme X. 
"""
mutable struct HomBundleSection
  parent::HomVecBundle
  loc_sections::Vector{MatSpaceElem}

  # optional fields for caching
  name::String

  function HomBundleSection( H::HomVecBundle, loc_sections::Vector{MatSpaceElem} )
    C = covering(H)
    n = length(patches(C))
    n == length(loc_sections) || error( "The number of local sections does not coincide with the number of patches of the covering" )
    for i in (1:n)
      parent(loc_sections[i]) == ambient_ring(patches(C)[i]) || error( "The $i-th local section is not an element of the ambient ring of the patch" )
    end
    return new( E, loc_sections )
  end
end

@Markdown.doc """
    function pullback( L::LineBundle, ρ::Refinement ) -> LineBundle

Pullback of a line bundle along a refinement. 
This returns an isomorphic line bundle but defined on a refinement 
of the original covering. 
"""
function pullback( L::LineBundle, ρ::Refinement )
  if covering(L) != original_cover(ρ) 
    error( "The line bundle is defined on a different cover than the one being refined" )
  end
  m = length( patches( covering( L )))
  n = length( patches( original_cover( ρ )))
  inclusions = inclusions( ρ )
  new_transitions = Array{MPolyElem,2}(undef,n,n)
  glueings_original = glueings(domain(ρ))
  glueings_refined = glueings(codomain(ρ))
  for i in (1:n-1)
    for j in (i+1:n)
      if i != j
	f_i = inclusions[i] # the inclusion of X_i into Y_{n(i)}
	phi_i = pullback(f_i) # the associated pullback on the level of rings
	S = codomain( phi_i )
	S_loc = ambient_ring(overlap(glueings_refined[i,j])) 
	f_j = inclusions[j] # the inclusion of X_j into Y_{n(j)}
	phi_j = pullback(f_j) # the associated pullback on the level of rings
	# We need to extract the pullback on the overlap 
	# X_i ∩ X_j ↪ Y_{n(i)} ∩ Y_{n(j)} in order to pull back the transition 
	# function. The order of i and j determines whether the intersection 
	# on the left hand side is realized via a localization of X_i or 
	# X_j. Similar on the other side, only that the ordering is not 
	# necessarily preserved by the function n, i.e. it can happen that 
	#   n(i) ≥ n(j) despite i < j. 
	n_i = map_index(ρ,i)
	n_j = map_index(ρ,j)

	if n_i == n_j 
	  # in case both patches are contained in a single one of the original cover, 
	  # just keep it simple.
	  new_transitions[i,j] = one(S_loc)
	elseif n_i < n_j 
	  # the next best case. Here, the overlap of Y_{n(i)} ∩ Y_{n(j)} 
	  # is realized via a localization of Y_{n(i)} and hence the 
	  # pullback can be inferred from X_i → Y_{n(i)}.
	  R = domain(phi) 
	  R_loc = ambient_ring(overlap(glueings_original[n_i,n_j]))
	  # phi_loc = 
	else
	  # In this last case, one needs to explicitly invert one 
	  # of the glueing isomorphisms.
	end 

      end
    end
  end
  M = LineBundle( parent(L), refined_cover(ρ), new_transitions )
  return M
end


@Markdown.doc """
    mutable struct VectorBundleSection 

Struct for a section in a VectorBundle. This will store a list of 
local sections in the distinct patches of the covering with respect 
to which the VectorBundle is defined. 
"""
mutable struct VectorBundleSection 
  parent::VectorBundle
  loc_sections::Vector{Vector{MPolyElem}}
  
  # optional variables used for caching
  name::String

  function VectorBundleSection( E::VectorBundle, loc_sections::Vector{Vector{MPolyElem}} )
    C = covering(E)
    n = length(patches(C))
    n == length(loc_sections) || error( "The number of local sections does not coincide with the number of patches of the covering" )
    for i in (1:n)
      parent(loc_sections[i]) == ambient_ring(patches(C)[i]) || error( "The $i-th local section is not an element of the ambient ring of the patch" )
    end
    return new( E, loc_sections )
  end
end

@Markdown.doc """
    function VectorBundleSection( s::VectorBundleSection )

Returns a copy of s.
"""
function VectorBundleSection( s::VectorBundleSection )
  return VectorBundleSection( s.parent, copy(s.loc_sections) )
end

function Base.copy( s::VectorBundleSection )
  return VectorBundleSection( s )
end

@Markdown.doc """
loc_sections(s::VectorBundleSection) -> Vector{Vector{MPolyElem}}

Returns a list of local sections of s on the distinct patches in 
the covering with respect to which the underlying vector bundle is 
defined.
"""
loc_sections(s::VectorBundleSection) = s.loc_sections
@Markdown.doc """
    parent(s::VectorBundleSection) -> VectorBundle

Returns the VectorBundle E of which s is a section.
"""
parent(s::VectorBundleSection) = s.parent

@Markdown.doc """
    function zero_locus( s::VectorBundleSection ) -> CoveredScheme

Returns the CoveredScheme Y obtained as the zero locus of the section s.

**WARNING: This is not implemented, yet!**
"""
function zero_locus( s::VectorBundleSection )
  error( "not implemented yet" )
  n = length(loc_sections(s))
  E = parent(s)
  C = covering(E)
  if rank(E) == 0 
    return parent(C)
  end
  U = patches(C)
  D = Covering()
  for i in (1:n)
    X = subscheme( U[i], ideal( ambient_ring(U[i]), s.loc_sections[i] ) )
    _add_patch!( D, X, glueings(C)[i,(1:i-1)] )
  end
  Y = CoveredScheme( D )
end

##########################################################
#
# Section in a line bundle
#
@Markdown.doc """
    mutable struct LineBundleSection 

Struct for a section in a LineBundle. This will store a list of 
local sections in the distinct patches of the covering with respect 
to which the LineBundle is defined. 
"""
mutable struct LineBundleSection 
  parent::LineBundle
  loc_sections::Vector{MPolyElem}
  
  # optional variables used for caching
  name::String

  function LineBundleSection( E::LineBundle, loc_sections::Vector{MPolyElem} )
    C = covering(E)
    n = length(patches(C))
    n == length(loc_sections) || error( "The number of local sections does not coincide with the number of patches of the covering" )
    for i in (1:n)
      parent(loc_sections[i]) == ambient_ring(patches(C)[i]) || error( "The $i-th local section is not an element of the ambient ring of the patch" )
    end
    return new( E, loc_sections )
  end
end

@Markdown.doc """
    function LineBundleSection( s::LineBundleSection )

Returns a copy of s.
"""
function LineBundleSection( s::LineBundleSection )
  return LineBundleSection( s.parent, copy(s.loc_sections) )
end

function Base.copy( s::LineBundleSection )
  return LineBundleSection( s )
end

@Markdown.doc """
loc_sections(s::VectorBundleSection) -> Vector{MPolyElem}

Returns a list of local sections of s on the distinct patches in 
the covering with respect to which the underlying vector bundle is 
defined.
"""
loc_sections(s::LineBundleSection) = s.loc_sections

@Markdown.doc """
    parent(s::LineBundleSection) -> LineBundle

Returns the LineBundle E of which s is a section.
"""
parent(s::LineBundleSection) = s.parent

function Base.show( io::Base.IO, s::LineBundleSection )
  if isdefined( s, :name )
    Base.print( s.name )
    return
  end
  Base.println( io, "Section in the line bundle $(parent(s)) on the scheme $(parent(parent(s)))" )
  Base.println( io, "given in the covering $(covering(parent(s))) with local sections" )
  Base.println( io, loc_sections(s))
end


@Markdown.doc """
    function zero_locus( s::LineBundleSection ) -> CoveredScheme

Returns the CoveredScheme Y obtained as the zero locus of the section s.

**WARNING: This is not implemented, yet!**
"""
function zero_locus( s::LineBundleSection )
  error( "not implemented yet" )
  n = length(loc_sections(s))
  E = parent(s)
  C = covering(E)
  if rank(E) == 0 
    return parent(C)
  end
  U = patches(C)
  D = Covering()
  for i in (1:n)
    X = subscheme( U[i], s.loc_sections[i] )
    _add_patch!( D, X, glueings(C)[i,(1:i-1)] )
  end
  Y = CoveredScheme( D )
end

@Markdown.doc """
    function VectorBundle( L::LineBundle )

Turns a LineBundle L into an object of type VectorBundle.
"""
function VectorBundle( L::LineBundle )
  C = covering(L)
  U = patches(C)
  n = length(U)
  transitions = Matrix{MatSpaceElem}(undef,n,n)
  for i in (1:n)
    for j in (1:n) 
      if i == j
        continue
      end
      a = L.transitions[i,j]
      M = MatrixSpace(parent(a),1,1)
      A = zero(M)
      A[1,1] = a
      transitions[i,j] = A
    end
  end
  return VectorBundle( 1, parent(L), covering(L), transitions )
end

@Markdown.doc """
    function tautological_sections( k::Ring, n::Int ) -> Vector{LineBundleSection}

Returns a vector with the canonical sections of 𝒪(1) on ℙⁿ corresponding 
to the homogeneous coordinates x₀,…xₙ in the ring 
S = k[x₀,…,xₙ] ≅ ⨁ ₖ H⁰(ℙⁿ, 𝒪(k)) of ℙⁿ.

**WARNING: The returned elements are not elements of any ring at this point!**
"""
function tautological_sections( k::Ring, n::Int )
  L = OO( k, n, 1 )
  IP = parent( L )
  C = covering(L)
  U = patches(C)
  result = LineBundleSection[]
  for i in (1:n+1)
    loc_sections = MPolyElem[]
    for j in (1:i-1)
      push!( loc_sections, gens(ambient_ring(U[j]))[i-1] )
    end
    push!( loc_sections, one(ambient_ring(U[i])) )
    for j in (i+1:n+1)
      push!( loc_sections, gens(ambient_ring(U[j]))[i] )
    end
    x = LineBundleSection( L, loc_sections )
    x.name = "x$(i-1)"
    push!( result, x )
  end
  return result
end

  
@Markdown.doc """
    function direct_sum( E::VectorBundle, F::VectorBundle )    

Return the direct sum E ⊕ F of the vector bundles E and F.

**WARNING: This is only implemented for vector bundles which 
are trivialized on the patches of the *same* covering.**
"""
function direct_sum( E::VectorBundle, F::VectorBundle )    
  covering(E) == covering(F) || error( "Direct sums of vector bundles on different coverings are not implemented, yet" )
  C = covering(E)
  r = rank(E)
  s = rank(F)
  U = patches(C)
  n = length(U)
  # Initialize an empty matrix for the transitions.
  G_transitions = Matrix{MatSpaceElem}(undef,n,n)
  for i in (1:n)
    for j in (1:n)
      if i == j 
	continue
      end
      G_transitions[i,j] = diagonal_matrix( [transitions(E)[i,j], transitions(F)[i,j]] )
    end
  end
  return VectorBundle( rank(E)+rank(F), parent(E), C, G_transitions )
end

direct_sum( L::LineBundle, E::VectorBundle ) = direct_sum( VectorBundle(L), E )
direct_sum( E::VectorBundle, L::LineBundle )  = direct_sum( E, VectorBundle(L) )
direct_sum( L::LineBundle, M::LineBundle) = direct_sum( VectorBundle(L), VectorBundle(M) )

@Markdown.doc """
    function cartesian_product( X::Spec, Y::Spec ) -> (::Spec, ::AffSchMorphims, ::AffSchMorphism)

Returns the product X ← X × Y → Y and its two projections.
"""
function cartesian_product( X::Spec, Y::Spec )
  #println( "Call to cartesian product with arguments $X and $Y" )
  base_ring(X) == base_ring(Y) || error( "The two given schemes are not defined over the same base ring" )
  k = base_ring(X)
  RX = ambient_ring(X)
  RY = ambient_ring(Y)
  m = length(gens(RX))
  n = length(gens(RY))
  R, prod_var = PolynomialRing( k, vcat( String.(symbols(RX)), String.(symbols(RY) )))
  pb_pi1 = AlgebraHomomorphism( RX, R, gens(R)[(1:m)] )
  pb_pi2 = AlgebraHomomorphism( RY, R, gens(R)[(m+1:n+m)] )
  IX = ideal( R, [ pb_pi1(g) for g in gens(defining_ideal(X)) ])
  IY = ideal( R, [ pb_pi2(g) for g in gens(defining_ideal(Y)) ])
  XxY = Spec( k, R, IX + IY )
  if isdefined( X, :name ) && isdefined( Y, :name )
    XxY.name = X.name * " × " * Y.name 
  end
  pi1 = AffSchMorphism( XxY, X, pb_pi1 )
  pi1.name = "π₁"
  pi2 = AffSchMorphism( XxY, Y, pb_pi2 )
  pi2.name = "π₂"
  return XxY, pi1, pi2
end

@Markdown.doc """
    function cartesian_product( U::SpecPrincipalOpen, Y::AffineScheme )

Returns the product U × Y and the two projections of the product of the root 
schemes of U and Y in their respective localization trees. 
The projections for U and Y can then be obtained by applying `restrict`.
"""
function cartesian_product( U::SpecPrincipalOpen, Y::AffineScheme )
  #println( "Call to cartesian product with arguments $U and $Y" )
  base_ring(U) == base_ring(Y) || error( "The two given schemes are not defined over the same base ring" )
  k = base_ring(U)
  X = parent(U)
  XxY, pi1, pi2 = cartesian_product( X, Y ) 
  UxY = localize( XxY, pullback(pi1)(denom(U)) )
  if isdefined( U, :name ) && isdefined( Y, :name )
    UxY.name = U.name * " × " * Y.name 
  end
  return UxY, pi1, pi2 
end

@Markdown.doc """
    function cartesian_product( X::Spec, V::SpecPrincipalOpen )

Returns the product X × V and the two projections of the product of X 
with the root scheme of Y in its localization tree.
The projection to V can then be obtained by applying `restrict`.
"""
function cartesian_product( X::Spec, V::SpecPrincipalOpen )
  #println( "Call to cartesian product with arguments $X and $V" )
  base_ring(X) == base_ring(V) || error( "The two given schemes are not defined over the same base ring" )
  k = base_ring(X)
  Y = parent(V)
  XxY, pi1, pi2 = cartesian_product( X, Y ) 
  XxV = localize( XxY, pullback(pi2)(denom(V)) )
  if isdefined( X, :name ) && isdefined( V, :name )
    XxV.name = X.name * " × " * V.name 
  end
  return XxV, pi1, pi2 
end

function cartesian_product( X::CoveredScheme, Y::CoveredScheme )
end

function cartesian_product( U::AffineScheme, Y::CoveredScheme )
  C = coverings(Y)[1]
  V = patches(C)
  n = length(V) 
  G = glueings(C)
  D = Covering() # the new covering of the product
  UxV = Vector{AffineScheme}() # the patches for the product 
  first_projections = Vector{AffSchMorphism}() # the first projections of the patches
  second_projections = Vector{AffSchMorphism}() # the second projections of the patches
  G_prod = Matrix{Union{Glueing,Nothing}}(undef,0,0) # the glueings for the product
  for j in (1:n) 
    W, p1, p2 = cartesian_product( U, V[j] )
    p1_res = restrict( p1, W, U )
    p2_res = restrict( p2, W, V[j] )
    G_new = Vector{Glueing}() # glueings for this particular patch to those already present
    push!( UxV, W )
    push!( first_projections, p1 )
    push!( second_projections, p2 )
    for i in (1:j-1)
      g_ij = G[i,j] # the glueing of these two patches 
      V_ij = overlap( g_ij )
      # the next two affine schemes are the ones to be glued together
      Wj_loc = localize( W, pullback(p2)(denom(domain(second_inclusion(g_ij)))) )
      Wi_loc = localize( UxV[i], pullback(second_projections[i])(denom(domain(first_inclusion(g_ij)))) )
      phi = glueing_isomorphism( g_ij )
      glueing_images = imgs_frac(phi)
      R = ambient_ring(root(W))
      Q = FractionField(R)
      imgs_U_vars = [Q(pullback(p1)(u)) for u in gens(ambient_ring(root(U)))]
      imgs_glueing = [ pullback(second_projections[i])(numerator(v))//pullback(second_projections[i])(denominator(v)) for v in glueing_images ]
      push!( G_new, Glueing( inclusion_in_parent(Wi_loc), 
			    inclusion_in_parent(Wj_loc), 
			    AffSchMorphism( Wi_loc, Wj_loc, vcat( imgs_U_vars, imgs_glueing ) )
			    ))
    end
    _add_patch!( D, W, G_new )
  end
  XxY = CoveredScheme( D )
  return XxY, 0, 0 # TODO: Implement that the two projections are also returned!
end


@Markdown.doc """
    function projectivize( E::VectorBundle ) -> ( Y::CoveredScheme, p::CovSchMorphism, L::LineBundle )

Constructs the projectivization Y = ℙ(E) of the VectorBundle E over X together 
with its projection p : ℙ(E) → X and the (relative) tautological bundle L = 𝒪(1) on 
ℙ(E).
"""
function projectivize( E::VectorBundle; fiber_var_name="u" )
  X = parent(E)
  k = base_ring(X)
  C = covering(E) 
  U = patches(C)
  G = glueings(C)
  n = length(U)
  r = rank(E)
  A = transitions(E)
  V = Vector{AffineScheme}()

# D = Covering()
# for i in (1:n)
#   UxIAr, p_base, p_fiber = cartesian_product( U[i], affine_space( k, r )) 
#   T = Vector{Glueing}()
#   for j in (1:i)
#     g = glueing_isomorphism(G[j,i])
#     T[j] = 

  #TODO: Finish the implementation.

  error( "Not implemented, yet" )
end

@Markdown.doc """
    projective_bundle_on_IP1(k::FieldType, a::Vector{Int}) where {FieldType<:Field}

Construct the bundle ℙ(𝒪(a₀)⊕𝒪(a₁)⊕ … ⊕ 𝒪(aᵣ)) → ℙₖ¹ as a CoveredScheme. 
"""
function projective_bundle_on_IP1(k::FieldType, a::Vector{Int}) where {FieldType<:Field}
  C = Covering()
  r = length(a)-1
  U = cartesian_product(affine_space(k,1,var_name="s"), projective_space(k,r,var_name="x"))
  V = cartesian_product(affine_space(k,1,var_name="t"), projective_space(k,r,var_name="x"))
  # Collect the 'diagonal' glueings
  diag_glueings = Vector{Glueing}()
  for i in (1:r+1)
    R_U = ambient_ring(U[i])
    s = gens(R_U)[1]
    U_loc = localize(U[i], s)
    R_V = ambient_ring(V[i])
    t = gens(R_V)[1]
    V_loc = localize(V[i], s)
    imgs_frac = [ 1//s ]
    imgs_frac = vcat(imgs_frac, [gens(R_U)[j+1]*(a[j]-a[i]>0 ? (1//s)^(a[j]-a[i]) : 1//(1//s)^(a[i]-a[j])) for j in (1:i-1)])
    imgs_frac = vcat(imgs_frac, [gens(R_U)[j+1]*(a[j]-a[i]>0 ? (1//s)^(a[j]-a[i]) : 1//(1//s)^(a[i]-a[j])) for j in (i+1:r)])
    phi = AffSchMorphism(U_loc, V_loc, imgs_frac)
    push!(diag_glueings, Glueing(inclusion_in_parent(U_loc), inclusion_in_parent(V_loc), phi))
  end
  # Now collect all the other glueings
  overall_glueings = Vector{Glueing}()
  G_U = glueings(covering(U)[1])
  G_V = glueings(covering(V)[1])
  for i in (1:r+1)
    R_Ui = ambient_ring(U[i])
    s = gens(R_Ui)[1]
    Ui_loc = localize(U[i], s)
    for j in (1:j-1) # in this range we have to use the inverse glueings on the other side.
      xj = gens(R_Ui)[j+1] # +1 for the base variable.
      Uij_loc = localize(Ui_loc, xj)
      g_U = inverse(G_U[j,i]) # the glueing of U[i] with U[j] for i>j, as is the case here.
      R_Vi = ambient_ring(V[i])
      t = gens(R_Vi)[1]
      Vi_loc = localize(V[i], t)
      yi = gens(R_Vi)[i+1]
      Vij_loc = localize(Vi_loc, yi)

      R_Vj = ambient_ring(V[j])
      t = gens(R_V)[1]
      Vj_loc = localize(V[j], s)
      y = gens(R_V)[i+1]
      U_lloc = localize(U_loc, x)
      V_lloc = localize(V_loc, y)
      #phi = restrict(glueing_isomorphism(diag_glueings[i]), U_lloc, V_lloc
      #push!(overall_glueings, Glueing( inclusion_in_root(U_lloc, V_lloc, compose(
    end
    push!(overall_glueings, diag_glueings[i])
    for j in (j+1:r+1)
    end
  end
end
