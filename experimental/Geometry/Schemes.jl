import AbstractAlgebra.Ring, Oscar.AlgHom, Oscar.compose, AbstractAlgebra.Generic.Frac
import Base: ∘
import Oscar: base_ring
import AbstractAlgebra: FPModule, FPModuleElem
import Base.copy

include( "./Multiindices.jl" )
using Oscar.Multiindices

include( "./Misc.jl" )
using Oscar.Misc

export AffineScheme, Spec, SpecPrincipalOpen
export base_ring, ambient_ring, defining_ideal, pullback, pullback_from_parent, pullback_from_root, inclusion_in_parent, inclusion_in_root, set_name!, inverted_element, identity_morphism, denoms, inverses
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

export AbstractCoherentSheaf

export IdealSheaf
export parent, covering, ideals

export CoherentSheaf
export parent, covering, modules

export LineBundle
export parent, covering, transitions, OO

export VectorBundleSection, LineBundleSection
export tautological_sections

abstract type Scheme{ S <: Ring }end
abstract type AffineScheme{S, T <: MPolyRing, U <: MPolyElem} <: Scheme{S} end
abstract type SchemeMorphism end

abstract type SchematicPoint end
abstract type Sheaf end
abstract type AbstractCoherentSheaf <: Sheaf end

####################################################################################
# The classical affine scheme, explicitly given as the quotient 
# ring of a polynomial algebra.
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
Spec( X::Spec ) = Spec( X.k, X.R, X.I )

Base.copy( X::Spec ) = Spec(X)

###################################################################################
# Getter functions
#
# No caching is needed in this case, since all these variables need to be assigned 
# at instantiation
#
function base_ring(A::Spec)
  return A.k
end

function ambient_ring(A::Spec)
  return A.R
end

function defining_ideal(A::Spec)
  return A.I
end

function denoms(A::Spec)
  return typeof(A.R)[]
end

function subscheme( X::Spec, J::MPolyIdeal ) 
  if base_ring(J) != X.R 
    error( "Ideal does not belong to the ambient ring of the variety" )
  end
  Y = Spec( X.k, X.R, X.I )
  Y.I = Y.I + J
  return Y
end

function subscheme( X::Spec, f::MPolyElem )
  return subscheme( X, ideal( ambient_ring(X), f ))
end

##################################################################################
# Affine scheme that arises as the localization of a Spec at a specific 
# element 'denom' of the coordinate ring, i.e. a principal open subset. 
#
# These are implemented in a recursive fashion meaning there is a parent 
# which is either a classical Spec or another principal open subset, from 
# which this instance is derived by inverting an element 'denom'. 
#
# Such parent structures form a tree which necessarily has a Spec at its 
# root. Currently, the implementation only allows localizations at elements 
# denom from the coordinate ring of this root. 
#
mutable struct SpecPrincipalOpen{S,T,U} <: AffineScheme{S,T,U}
  parent::AffineScheme{S,T,U}
  denom::U  # element of the "ambient" polynomial ring of the root
  	    # which has been localized in the parent to get this scheme.

  # Fields for caching. These provide the data as a Spec for this 
  # affine scheme
  k::S  # the ground field
  R::T	# the ambient polynomial ring
  I::MPolyIdeal{U} # the defining ideal
  pullbackFromParent::AlgHom
  pullbackFromRoot::AlgHom
  u::U # The inverse of denom in R
  name::String
  #inclusionInParent::AffSchMorphism
  #inclusionInRoot::AffSchMorphism

  
  function SpecPrincipalOpen(parent::Union{Spec{S,T,U},SpecPrincipalOpen{S,T,U}}, denom::U) where {S <: Ring, T<:MPolyRing, U<:MPolyElem}
    x = new{S,T,U}() 
    x.parent = parent
    # TODO: Implement a plausibility check: Does denom belong to the coordinate ring of the root?
    x.denom = denom
    return x 
  end 

  function SpecPrincipalOpen( k::S, R::T, I::MPolyIdeal{U} ) where {S <: Ring, T<:MPolyRing, U<:MPolyElem}
    if k != base_ring(R) 
      error( "basering not compatible with ambient ring" )
    end
    if base_ring(I) != R
      error( "ideal does not live in the given ambient ring" )
    end
    X = new{S,T,U}(k,R,I)
    return X
  end
end

# Copy constructor
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

###############################################################################
# Getter functions

function localize(parent, denom)
  return SpecPrincipalOpen(parent, denom)
end

function parent(S::AffineScheme)
  return S
end

# Check whether X is a parent of U and return the inclusion morphism 
# if true. Otherwise, return nothing.
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
    


function parent(D::SpecPrincipalOpen)
  return D.parent
end 

function root(D::SpecPrincipalOpen)
  return root(parent(D))
end

root(A::Spec) = A

function base_ring(D::SpecPrincipalOpen)
  if isdefined(D,:k)
    return D.k
  end
  k = base_ring(root(D))
  D.k = k
  return k
end

# Returns the unit that has last been adjoint to the 
# coordinate ring, i.e. the inverse of denom
function inverted_element(D::SpecPrincipalOpen)
  if isdefined( D, :u )
    return D.u
  end
  ambient_ring(D)
  return D.u
end

function denom(D::SpecPrincipalOpen)
  if isdefined( D, :denom )
    return D.denom
  end
  return nothing
end

# Return a vector containing the variables of the root's ambient 
# ring, but as elements of the new ambient ring.
function root_variables( D::SpecPrincipalOpen )
  phi = pullback_from_root(D)
  return [ phi(x) for x in gens(ambient_ring(root(D))) ]
end
  


# Setter functions

function set_name!( X::AffineScheme, name::String )
  X.name = name
end

# Obtain a closed subscheme
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

function subscheme( 
    X::SpecPrincipalOpen{S,T,U}, 
    f::U
  ) where{ S<:Ring, T<:MPolyRing, U<:MPolyElem }
  return subscheme( X, ideal( ambient_ring(root(X)), f ))
end

# Collect the denominators from this localization up to the root. 

function denoms(D::SpecPrincipalOpen{S,T,U}) where {S <: Ring, T <: MPolyRing, U<:MPolyElem}
  result = U[]
  P = D
  while typeof(P) <: SpecPrincipalOpen
    # Make sure that the denominators are collected in the same manner 
    # as the variables for the inverses appear.
    pushfirst!( result, P.denom )
    P = parent(P)
  end
  return result
end 

# Collect the inverses for this localization up to the root.

function inverses(D::SpecPrincipalOpen{S,T,U}) where {S <: Ring, T <: MPolyRing, U<:MPolyElem}
  result = [inverted_element(D)]
  R = ambient_ring(D)
  phi = pullback_from_parent(D)
  P = parent(D)
  while typeof(P) <: SpecPrincipalOpen
    pushfirst!(result, phi(inverted_element(P)))
    phi = compose( pullback_from_parent(P), phi )
    P = parent(P)
  end
  return result
end 


# Set up an ambient ring for this localization. 
#
# The current implementation adds a variable 'denom$i' for 
# every element whose inverse has been added to the 
# coordinate ring of the root. Again, the implementation is 
# recursive, but flattened in the sense that all variables
# are on the same level over the base ring.
 
function ambient_ring(D::SpecPrincipalOpen)
  if isdefined( D, :R )
    return D.R
  end
  if typeof(parent(D))<:Spec 
    R0 = ambient_ring(parent(D))
    I0 = defining_ideal(parent(D))
    R, phi, t = add_variables( R0, ["t"] ) 
    D.R = R
    D.pullbackFromParent = phi
    D.pullbackFromRoot = phi
    I = ideal( R, [ phi(g) for g in gens(I0) ] ) + ideal( R, [ one(R) - t[1]*phi(D.denom) ] )
    D.u = t[1]
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
    d = denoms(D)
    r = length(d)
    f = copy(D.denom)
    u = last(gens(R))
    phi = parent(D).pullbackFromRoot # We can recycle that since there's no t involved.
    for d in denoms(D)
      g = gcd(f,d)
      u *= phi(div(d,g))
      f = div(f,g) # TODO: Replace by the in-place call div! once it's there.
    end
    I = ideal( R, [ phi(g) for g in gens(root(D).I) ] ) + ideal( R, [ one(R) - last(gens(R))*phi(prod(d)*f) ] )
    images = copy(gens(R))
    push!( images, pop!(images)*phi(f) )
    D.pullbackFromParent = AlgebraHomomorphism( R, R, images )
    D.pullbackFromRoot = phi
    D.u = u
  end
  return R
end

function pullback_from_root( D::SpecPrincipalOpen )
  if isdefined( D, :pullbackFromRoot )
    return D.pullbackFromRoot
  end
  ambient_ring( D ) # This also stores the homomorphism
  return D.pullbackFromRoot
end

pullback_from_root( X::Spec ) = AlgebraHomomorphism( X.R, X.R, gens(X.R) )
pullback_from_parent( X::Spec ) = AlgebraHomomorphism( X.R, X.R, gens(X.R) )


function pullback_from_parent( D::SpecPrincipalOpen )
  if isdefined( D, :pullbackFromParent )
    return D.pullbackFromParent
  end
  ambient_ring( D ) # This also stores the homomorphism
  return D.pullbackFromParent
end

#function ring_hom(D::SpecPrincipalOpen)
#  R = ambient_ring(root(D))
#  S = ambient_ring(D)
#  d = length(denoms(D))
#  n = length(gens(R))
#  return AlgebraHomomorphism(R,S, gens(S)[d+1:n+d])
#end

function defining_ideal(D::SpecPrincipalOpen)
  if isdefined( D, :I )
    return D.I
  end
  ambient_ring(D) # This now also prepares the ideal.
  return D.I
end

function inclusion_in_parent( D::SpecPrincipalOpen )
  # return AffSchMorphism( D, parent(D), pullback_from_parent(D) )
  F = FractionField(ambient_ring(root(D)))
  return AffSchMorphism( D, parent(D), [ F(x) for x in gens(ambient_ring(root(D))) ] )
end

function inclusion_in_root( D::SpecPrincipalOpen )
  # return AffSchMorphism( D, root(D), pullback_from_root(D) )
  F = FractionField(ambient_ring(root(D)))
  return AffSchMorphism( D, root(D), [ F(x) for x in gens(ambient_ring(root(D))) ] )
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

domain(f::AffSchMorphism) = f.domain
codomain(f::AffSchMorphism) = f.codomain

# The copy constructor
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
  #                       h*
  #    Q_loc := Q[g^{-1}] ← P[f^{-1}] =: P_loc
  #                 ↑          ↑
  #                 Q          P
  # where both P = R/I and Q = S/J are both quotients of 
  # polynomial rings and the sets of elements g = g_1,...,g_ν, 
  # g_j ∈ Q and f = f_1,...,f_μ, f_i ∈ P have 
  # been inverted to arrive at the open sets U and V.
  # The ring P_loc is realized as a quotient of a polynomial 
  # ring as indicated in the diagram 
  #      R[s_1,...,s_μ] 
  #            ↓        ↘ π 
  #      P[s_1,...,s_μ] → P_loc, s_i ↦ 1/f_i
  # and similar for Q_loc.
  U = domain(h)
  V = codomain(h)
  X = root(U)
  Y = root(V)
  R_s = ambient_ring(V)
  S_t = ambient_ring(U)
  R = ambient_ring(Y)
  S = ambient_ring(X)
  @show U
  @show V
  @show X
  @show Y
  @show R_s
  @show S_t
  @show R
  @show S
  # Let x₁,...,xₙ be the variables of R.
  # Now the fractional representation of the morphism, i.e. 
  # the elements in imgs_frac, describe the images of these 
  # variables in Q_loc: 
  #            pᵢ
  #   xᵢ ↦  -------          
  #          b^{kᵢ}
  # for some polynomials p_i and non-negative integers k_i.
  # We need to describe the homomorphism 
  #
  #  ϕ : R[s_1,...,s_μ] → S[t_1,...,t_ν]
  # 
  # lifting the map h*. The information we have so far 
  # suffices to determine the map 
  #
  #  ϕ' : R → S[t_1,...,t_ν], x_i ↦ p_i ⋅ t_i^{k_i},
  #
  # but we also need to determine the images of the 
  # variables s_i. Theoretically we find 
  #
  #  π∘ϕ( s_i ) = 1/π∘ϕ'(f_i) 
  #
  # but the right hand side does not yet have an admissible 
  # denominator (i.e. a product of powers of the g_j).
  # But according to the assumption that h(U) ⊂ V, 
  # or equivalently h^{-1}(Y∖ V) ⊂ X ∖ U, we necessarily 
  # have a relation 
  # 
  #  π(a_i) ⋅ π∘ϕ'(f_i) = π(g)^{k_i}
  # 
  # for some polynomials a_i and integers k_i. Here, g is 
  # the product of all the g_j's and all relevant projections 
  # to the quotients are denoted by π. With the help of these 
  # relations we can extend the above equalities towards 
  #
  #  π∘ϕ( s_i ) = 1/π∘ϕ'(f_i) = π(a_i)/π(g^{k_i}) = π(a_i)⋅t^{k_i}
  #
  # with, again, t the product of variables t_j.

  println( "Call to pullback" )
  @show R
  @show S
  for f in h.imgs_frac
    @show f
  end
  n = length( h.imgs_frac )
  if n != length( gens( R ))
    error( "Number of variables in the ambient ring of the root does not coincide with the number of images provided for the homomorphism" )
  end
  images = Vector{elem_type(S_t)}()
  den_U = denoms(U)
  inv_U = inverses(U)
  Q, proj = quo( S, defining_ideal(X) )
  g = prod( den_U )
  t = prod( inv_U )
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
  #phi_prime = AlgebraHomomorphism( R, S_t, images )
  # TODO: go back to the old line, once the bug is fixed.
  phi_prime = AlgebraHomomorphism( R, S_t, copy(images) )

  # In case Y = V we are done
  if typeof(V)<:Spec 
    h.pullback = phi_prime
    return phi_prime
  end
  
  # Now assemble the images for the remaining variables s_i 
  # of R_loc along the same pattern.
  den_V = denoms(V)
  inv_V = inverses(V)
  mu = length(inv_V) # Problem: Das ist nicht notwendigerweise die Anzahl 
  # der notwendigen Variablen! TODO: Fix it!
  Q_loc, proj_loc = quo( S_t, defining_ideal(U))
  for i in (1:mu) 
    q = phi_prime(den_V[i])
    # TODO: This necessarily computes a Groebner basis for the localized 
    # ideal. We probably don't want that by default, but look for a computationally 
    # less expensive way to solve this problem.
    (k,a) = divides_power( proj_loc(q), proj_loc(pullback_from_root(U)(g)) )
    if k == -1
      error( "the homomorphism is not well defined" )
    end
    push!( images, lift(a)*t^k )
  end
  @show R_s
  @show S_t
  @show images
  phi = AlgebraHomomorphism( R_s, S_t, images )
  h.pullback = phi
  return phi
end

# Construct the fractional representation of the morphism 
# from the explicit algebra homomorphism. 
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
  den = denoms(domain(f))
  d = length(den)
  n = length(gens(P))
  # prepare a list of the images of the generators 
  fracs = [ F(x) for x in gens(P) ]
  fracs = vcat( fracs, [ 1//g for g in den ])
  R = ambient_ring(root(codomain(f)))
  imgs_frac = elem_type(F)[]
  for g in gens(R)
    h = compose( phi, pullback_from_root(codomain(f)) )
    push!(imgs_frac,evaluate(h(g), fracs))
  end
  f.imgs_frac = imgs_frac
  return imgs_frac
end

function identity_morphism( X::Spec )
  R = ambient_ring(X)
  phi = AlgebraHomomorphism( R, R, gens(R) )
  f = AffSchMorphism(X,X,phi)
  return f 
end

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

#####################################################################
# Glueing of two affine schemes X and Y along an open subset 
# U ⊂ X and V ⊂ Y via an isomorphism Phi : U → V. 
# A glueing consists of a diagram X ↩ U → V ↪ Y with the outer maps 
# the open embeddings of principal open subsets. 
#
mutable struct Glueing
  firstPatch::AffineScheme
  secondPatch::AffineScheme
  inclusionFirstPatch::AffSchMorphism
  inclusionSecondPatch::AffSchMorphism
  glueingIsomorphism::AffSchMorphism

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


first_patch( G::Glueing ) = G.firstPatch
second_patch( G::Glueing ) = G.secondPatch
overlap(G::Glueing) = domain(G.inclusionFirstPatch)
first_inclusion( G::Glueing ) = G.inclusionFirstPatch
second_inclusion( G::Glueing ) = G.inclusionSecondPatch
glueing_isomorphism( G::Glueing ) = G.glueingIsomorphism

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
    C.patches = []
    C.glueings = []
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

function patches( C::Covering )
  if isdefined( C, :patches )
    return C.patches
  else 
    return nothing 
  end 
end

function patch( C::Covering, i::Int )
  if isdefined( C, :patches )
    return( C.patches[i] )
  else
    return nothing
  end
end

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

########################################################################
# add one affine patch X to an existing covering C.
# To this end, the glueings along all the affine patches Y
# in the covering C to X need to be provided. They are 
# passed on as a vector using the same indices as the 
# list of patches in the covering C.

function _add_patch!( C::Covering, X::AffineScheme, glueings::Vector{Glueing} )
  push!( C.patches, X )
  n = length(glueings)
  if n != length(patches(C))-1
    error( "the number of glueings does not coincide with the number of patches so far" )
  end
  for i in n
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

#####################################################################
# intersection of two patches X and Y in a given covering.
# Returns the intersection as an affine scheme 
# together with the two open inclusions.
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


###############################################################
# Restriction of a morphism to open inclusions:
#    f
#  X → Y
# i↑   ↑j
#  U → V
#    f'
# When f, i, and j are given where X and Y are affine schemes 
# and i and j are inclusions of principal open subsets with 
# f(U) ⊂  V, we would like to explicitly infer the restriction 
# f' of f to U.
function restrict( f::AffSchMorphism, U::SpecPrincipalOpen, V::SpecPrincipalOpen )
  # TODO: implement this.
end

##############################################################
# Take the union of two affine patches in a given covering.
function union( X::AffineScheme, Y::AffineScheme, C::Covering )
  # TODO: This is not functional yet. We first need a struct 
  # for morphisms of CoveredScheme.
end

#######################################################################
#
# Refinement of a covering C by another covering D.
# Let I and J be the index sets for the patches of the coverings 
# C and D, respectively. Then a refinement consists of 
# a function 
#   n : J → I 
# such that for every patch Y_j in D there exists a patch X_i in C with 
# i = n(j) such that X_i contains Y_j. 
# The datum of the refinement stored here is then the respective collection 
# of morphisms of affine schemes given by these containments. 
#
mutable struct Refinement
  original_cover::Covering
  refined_cover::Covering
  index_map::Vector{Int}
  inclusions::Vector{AffSchMorphism}

  # TODO: Implement the constructors.
end

function original_cover( ρ::Refinement )
  if isdefined( ρ, :original_cover )
    return ρ.original_cover
  end
  return nothing
end

function refined_cover( ρ::Refinement )
  if isdefined( ρ, :refined_cover )
    return ρ.refined_cover
  end
  return nothing
end

function map_index( ρ::Refinement, i::Int )
  return ρ.index_map[i]
end

function inclusions( ρ::Refinement )
  if isdefined( ρ, :inclusions ) 
    return ρ.inclusions
  end
  return nothing
end


######################################################################
# A scheme explicitly described by providing a covering 
# with glueings. 
#
# In general, a scheme might have different coverings since, 
# for example, different vector bundles might be given by means of 
# different local trivializations. 
#
# In practice, however, we usually expect to have a graph of coverings 
# ordered by refinements and with a single root accounting for the covering 
# originally used to describe our scheme in the first place. 
#
mutable struct CoveredScheme
  coverings::Vector{Covering}
  refinements::Array{Refinement,2}

  # fields for caching
  name::String

  function CoveredScheme()
    X = new()
    coverings = Vector{Covering}[]
    refinements
    return X 
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

function coverings( X::CoveredScheme )
  return X.coverings
end

function set_name!( X::CoveredScheme, name::String )
  X.name = name
  return name
end

function name( X::CoveredScheme ) 
  if isdefined( X, :name )
    return X.name
  end
  return nothing
end


function Base.show( io::Base.IO, X::CoveredScheme )
  if isdefined( X, :name )
    Base.print( io, X.name )
    return
  end
  Base.print( io, "Covered Scheme with $(length(X.coverings)) covering." )
end

####################################################################
#
# Construct projective space IP^n_k over the ring k as a covered 
# scheme with its standard open cover
#
function projective_space( k::Ring, n::Int )
  # Start the standard covering with the first patch for x_0 ≠ 0.
  # println( "assembling $n-dimensional projective space over $k." )
  X0 = affine_space(k,n)
  set_name!( X0, "IA^$(n)_0" )
  C = Covering(X0)
  # Add the other patches successively
  for j in (1:n) 
    # println("Variable in the outer loop j = $j")
    # The patch for x_j ≠ 0
    # println("Glueing in the following new patch:")
    Y = affine_space( k, n )
    set_name!(Y,"IA^$(n)_$(j)")
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
  set_name!( IP, "IP^$n" )
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

# The copy constructor.
# Refer to the same Scheme and covering, but copy the 
# ideals.
function IdealSheaf( I::IdealSheaf )
  return IdealSheaf( I.parent, I.covering, copy(I.ideals) )
end

parent( I::IdealSheaf ) = I.parent
covering( I::IdealSheaf ) = I.covering
ideals( I::IdealSheaf ) = I.ideals

# Take a homogeneous ideal I in a graded homogeneous 
# ring S and turn it into an ideal sheaf on the 
# corresponding projective space
function to_ideal_sheaf( I::MPolyIdeal )
  error( "not implemented yet" )
  # TODO: Implement this based on graded modules in 
  # Oscar. The routine is supposed to recognize automatically 
  # from the parent how to construct the associated projective 
  # space and return the ideal sheaf on it. 
end

#########################################################
#
# A sheaf of modules on a CoveredScheme X
#
# Everything is similar to the ideal sheaf above, 
# only that instead of ideals we store modules over 
# the chosen ambient rings of the patches
#
mutable struct CoherentSheaf <: AbstractCoherentSheaf 
  parent::CoveredScheme
  covering::Covering
  modules::Vector{FPModule{MPolyElem}}

  # optional variables used for caching
  name::String
  
  # TODO: Implement the constructors.
end

parent( I::CoherentSheaf ) = I.parent
covering( I::CoherentSheaf ) = I.covering
modules( I::CoherentSheaf ) = I.modules


########################################################
# 
# A line bundle on a CoveredScheme X
#
# Again, this refers to a specific covering C of X 
# which has to be fine enough so that the line bundle is 
# trivialized on the patches. Then, the bundle itself 
# is described by an "antisymmetric" matrix of transition 
# functions. 
#
# Note that not every line bundle on an affine 
# scheme is necessarily trivial! But in general one will 
# need to even cover an affine scheme by finer patches 
# in order to represent it as a line bundle of this form.
#
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

function LineBundle( L::LineBundle )
  return LineBundle( L.parent, L.covering, copy(L.transitions) )
end

function Base.copy( L::LineBundle )
  return LineBundle( L )
end

parent( L::LineBundle ) = L.parent
covering( L::LineBundle ) = L.covering
transitions( L::LineBundle ) = L.transitions

########################################################
# The tautological bundle on projective space
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

OO( n::Int, d::Int ) = OO( QQ, n, d )

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

#######################################################
# 
# An algebraic vector bundle on a CoveredScheme X
#
# Basically the same as LineBundle above. Only the transition 
# functions are now taking values in the space of invertible 
# matrices.
#
mutable struct VectorBundle <: AbstractCoherentSheaf
  rank::Int
  parent::CoveredScheme
  covering::Covering
  transitions::Array{MatrixElem{MPolyElem},2}
  
  # optional variables used for caching
  name::String

  function VectorBundle( 
      rank::Int,
      parent::CoveredScheme, 
      covering::Covering, 
      transitions::Array{MatrixElem{MPolyElem},2} )
    E = new()
    E.rank = rank
    E.parent = parent
    E.covering = covering
    E.transitions = transitions
    # TODO: Implement the cheap plausibility checks.
    return E
  end
end

function VectorBundle( E::VectorBundle )
  return VectorBundle( E.rank, E.parent, E.covering, copy(E.transitions) )
end

function Base.copy( E::VectorBundle )
  return VectorBundle( E )
end

parent( L::VectorBundle ) = L.parent
covering( L::VectorBundle ) = L.covering
transitions( L::VectorBundle ) = L.transitions


#########################################################
# pullback of a line bundle along a refinement
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


##########################################################
#
# Section in a vector bundle
#
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

function VectorBundleSection( s::VectorBundleSection )
  return VectorBundleSection( s.parent, copy(s.loc_sections) )
end

function Base.copy( s::VectorBundleSection )
  return VectorBundleSection( s )
end

loc_sections(s::VectorBundleSection) = s.loc_sections
parent(s::VectorBundleSection) = s.parent

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

function LineBundleSection( s::LineBundleSection )
  return LineBundleSection( s.parent, copy(s.loc_sections) )
end

function Base.copy( s::LineBundleSection )
  return LineBundleSection( s )
end

loc_sections(s::LineBundleSection) = s.loc_sections
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

# The global sections of the dual of the tautological bundle 
# on IP^n which correspond to the homogeneous coordinates.
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

  
    

