import AbstractAlgebra.Ring, Oscar.AlgHom, Oscar.compose, AbstractAlgebra.Generic.Frac
import Base: ∘
import Oscar: base_ring
import AbstractAlgebra: FPModule, FPModuleElem

include( "./Multiindices.jl" )
using Oscar.Multiindices

include( "./Misc.jl" )
using Oscar.Misc

export AffineScheme, Spec, SpecPrincipalOpen, affine_space
export AffSchMorphism, base_ring,ambient_ring,defining_ideal, pullback, pullback_from_parent, pullback_from_root, inclusion_in_parent, inclusion_in_root, set_name!
export localize, get_unit

export AffSchMorphism
export domain, codomain, pullback

export Glueing
export first_patch, second_patch, overlap, first_inclusion, second_inclusion, glueing_isomorphism

export Covering
export get_patches, get_glueings, add_patch

export Refinement

export CoveredScheme
export get_coverings

export AbstractCoherentSheaf

export IdealSheaf
export parent, covering, ideals

export CoherentSheaf
export parent, covering, modules

export LineBundle
export parent, covering, transitions, OO

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
  u::U # The inverse of denom in R/I (one of the variables)
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
end

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
  ι = get_identity(U)
  if X == U
    return ι
  end
  V = U
  while typeof(V)<:SpecPrincipalOpen
    W = parent(V)
    ι = compose( inclusion_in_parent(V), ι )
    if W == X 
      return ι
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
# coordinate ring, i.e. the inverse of the inverted element
function get_unit(D::SpecPrincipalOpen)
  if isdefined( D, :u )
    return D.u
  end
  base_ring(D)
  return D.u
end

function get_denom(D::SpecPrincipalOpen)
  if isdefined( D, :denom )
    return D.denom
  end
  return nothing
end

# Return a vector containing the variables of the root's ambient 
# ring, but as elements of the new ambient ring.
function get_root_variables( D::SpecPrincipalOpen )
  ϕ = pullback_from_root(D)
  return [ ϕ(x) for x in gens(ambient_ring(root(D))) ]
end
  


# Setter functions

function set_name!( X::AffineScheme, name::String )
  X.name = name
end

##################################################################
# 
# Localize an affine scheme at an element f from its root
#
function hypersurface_complement( X::AffineScheme, f::MPolyElem )
  if parent( f ) != ambient_ring( root( X ))
    error( "The polynomial is not an element of the root ambient ring of the given affine scheme" )
  end
  #TODO: Complete this. Maybe it's redundant with the localize routine above, though.
end

# Collect the denominators from this localization up to the root. 

function denoms(D::SpecPrincipalOpen{S,T,U}) where {S <: Ring, T <: MPolyRing, U<:MPolyElem}
  result = U[]
  P = D
  while typeof(P) <: SpecPrincipalOpen
    push!(result, P.denom)
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
  #names = ["denom$i" for i in 1:length(denoms(D))]
  n = length( denoms( D ))
  R = ambient_ring(parent(D))
  S, ϕ, u = add_variables( R, ["u$n"] )
  D.R = S
  D.u = u[1]
  D.pullbackFromParent = ϕ
  if typeof(parent(D)) <: Spec 
    D.pullbackFromRoot = D.pullbackFromParent
  else 
    D.pullbackFromRoot = compose( pullback_from_root( parent( D )), ϕ )
  end
  return S
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

#function get_ring_hom(D::SpecPrincipalOpen)
#  R = ambient_ring(root(D))
#  S = ambient_ring(D)
#  d = length(denoms(D))
#  n = length(gens(R))
#  return AlgebraHomomorphism(R,S, gens(S)[d+1:n+d])
#end

function defining_ideal(D::SpecPrincipalOpen)
  ϕ = pullback_from_parent(D)
  ψ = pullback_from_root(D)
  I = defining_ideal(parent(D))
  R = ambient_ring(D)
  J = ideal( R, [ ϕ(f) for f in gens( I ) ])
  J = J + ideal( R, [ one(R)-D.u*ψ(D.denom) ] )
  return ( J )
end

function inclusion_in_parent( D::SpecPrincipalOpen )
  return AffSchMorphism( D, parent(D), pullback_from_parent(D) )
end

function inclusion_in_root( D::SpecPrincipalOpen )
  return AffSchMorphism( D, root(D), pullback_from_root(D) )
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
      codomain::AffineScheme{S,Tc,Uc}, imgs_frac::Vector{Frac{Uc}}
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
end

domain(f::AffSchMorphism) = f.domain
codomain(f::AffSchMorphism) = f.codomain

# Construct the pullback on the level of coordinate rings 
# from the fractional representation 
function pullback(f::AffSchMorphism)
  if isdefined(f, :pullback)
    return f.pullback
  end
  if !isdefined(f, :imgs_frac )
    error( "Neither the fractional representation, nor the pullback is defined for this morphism" )
  end
  S = ambient_ring(codomain(f))
  T = ambient_ring(domain(f))
  R = ambient_ring(root(codomain(f)))
  # TODO reconstruct the ring homomorphism S -> T from the fractional representation 
  # and describe it in terms of R -> T with R = ambient_ring(root(codomain(f)))

  n = length( f.imgs_frac )
  if n != length( gens( R ))
    error( "Number of variables in the ambient ring of the root does not coincide with the number of images provided for the homomorphism" )
  end
  for i in (1:n)
    p = numerator( f.imgs_frac[i] )
    q = denominator( f.imgs_frac[i] )
    u = prod( denoms( domain(f)))
    (k,a) = coeffs_in_radical( q, u )
    
  end
  #TODO: Finish this once we know how to deal with the radical membership in Oscar.

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
  ϕ = pullback(f)
  # Set up the codomain of a lift of the pullback ϕ to the ambient ring of 
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
    h = compose( ϕ, pullback_from_root(codomain(f)) )
    push!(imgs_frac,evaluate(h(g), fracs))
  end
  f.imgs_frac = imgs_frac
  return imgs_frac
end

function get_identity( X::Spec )
  R = ambient_ring(X)
  ϕ = AlgebraHomomorphism( R, R, gens(R) )
  f = AffSchMorphism(X,X,ϕ)
  return f 
end

function get_identity( D::SpecPrincipalOpen )
  R = ambient_ring(D)
  ϕ = AlgebraHomomorphism( R, R, gens(R) )
  f = AffSchMorphism(D,D,ϕ)
  return f
end

#####################################################################
# Composition of affine scheme morphisms. 
# Consists basically of composition of the associated pullback 
# homomorphisms.
#
function compose( f::AffSchMorphism, g::AffSchMorphism ) 
  if codomain(f) != domain(f) 
    error( "morphisms can not be composed" )
  end
  ϕ = pullback(f)
  γ = pullback(g)
  return AffSchMorphism( domain(g), codomain(f), compose( γ, ϕ ) )
end

#####################################################################
# Glueing of two affine schemes X and Y along an open subset 
# U ⊂ X and V ⊂ Y via an isomorphism Φ : U → V. 
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
    Γ = new()
    Γ.firstPatch = codomain( inclusionFirstPatch )
    Γ.secondPatch = codomain( inclusionSecondPatch )
    Γ.inclusionFirstPatch = inclusionFirstPatch
    Γ.inclusionSecondPatch = inclusionSecondPatch
    Γ.glueingIsomorphism = glueingIsomorphism
    return Γ
  end
end

first_patch( Γ::Glueing ) = Γ.firstPatch
second_patch( Γ::Glueing ) = Γ.secondPatch
overlap(Γ::Glueing) = domain(Γ.inclusionFirstPatch)
first_inclusion( Γ::Glueing ) = Γ.inclusionFirstPatch
second_inclusion( Γ::Glueing ) = Γ.inclusionSecondPatch
glueing_isomorphism( Γ::Glueing ) = Γ.glueingIsomorphism

function Base.show( io::Base.IO, Γ::Glueing )
  Base.println( io, "Glueing of" )
  Base.println( io, Γ.firstPatch )
  Base.println( io, "and" )
  Base.println( io, Γ.secondPatch )
  Base.println( io, "along" )
  Base.println( io, domain(Γ.inclusionFirstPatch) )
  return
end

################################################################################
# Covering of a scheme.
# It consists of a list of affine patches together with their glueings.
# The glueings are stored in an array (which is usually half empty, since 
# glueings need only be provided for ordered pairs of indices). 
# The get_index methods for double arrays have been overloaded so that 
# OrderedMultiindex can be used to address them.
#
mutable struct Covering
  patches::Vector{AffineScheme}	# A list of open patches for this scheme
  glueings::Matrix{Union{Glueing,Nothing}}     # A list of glueings between the patches 
  						# listed above
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

function get_patches( C::Covering )
  if isdefined( C, :patches )
    return C.patches
  else 
    return nothing 
  end 
end

function get_patch( C::Covering, i::Int )
  if isdefined( C, :patches )
    return( C.patches[i] )
  else
    return nothing
  end
end

function get_glueings( C::Covering )
  if isdefined( C, :glueings )
    return C.glueings
  else
    return nothing
  end
end

########################################################################
# add one affine patch X to an existing covering C.
# To this end, the glueings along all the affine patches Y
# in the covering C to X need to be provided. They are 
# passed on as a vector using the same indices as the 
# list of patches in the covering C.

function add_patch( C::Covering, X::AffineScheme, glueings::Vector{Glueing} )
  push!( C.patches, X )
  n = length(glueings)
  #@show n
  #@show length( get_patches(C))
  if n != length(get_patches(C))-1
    error( "the number of glueings does not coincide with the number of patches so far" )
  end
  for i in n
    if first_patch( glueings[i] ) != get_patches(C)[i] 
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
  if !(X in get_patches(C)) 
    error( "X is not a patch in the covering C" )
  end
  if !(Y in get_patches(C)) 
    error( "Y is not a patch in the covering C" )
  end
  if X == Y 
    return ( X, get_identity(X), get_identity(X) )
  end
  i = indexin( X, get_patches(C) )[1]
  j = indexin( Y, get_patches(C) )[1]
  # Keep in mind that the glueing is realized 
  # in a directed way. Thus, we need to make a distinction.
  if i < j
    G = get_glueings(C)[i,j]
    return (overlap(G),first_inclusion(G),second_inclusion(G))
  else
    G = get_glueings(C)[j,i]
    return (overlap(G),second_inclusion(G),first_inclusion(G))
  end
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
end

function get_original_cover( ρ::Refinement )
  if isdefined( ρ, :original_cover )
    return ρ.original_cover
  end
  return nothing
end

function get_refined_cover( ρ::Refinement )
  if isdefined( ρ, :refined_cover )
    return ρ.refined_cover
  end
  return nothing
end

function map_index( ρ::Refinement, i::Int )
  return ρ.index_map[i]
end

function get_inclusions( ρ::Refinement )
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

  function CoveredScheme( C::Covering ) where{ S <: Ring }
    X = new()
    X.coverings = [C]
    #=
    if length( get_patches(C) ) == 0 
      error( "empty covering. Cannot determine base ring." )
      return
    end
    X.base_ring = base_ring( get_patches(C)[1] )
    =#
    return X
  end

end

function get_coverings( X::CoveredScheme )
  return X.coverings
end

function set_name!( X::CoveredScheme, name::String )
  X.name = name
  return name
end

function get_name( X::CoveredScheme ) 
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
    #@show Y
    # come up with the glueings
    G = Vector{Glueing}()
    # println("New glueings for this step so far:") 
    #@show G
    #@show typeof( G )
    for i in (0:j-1)
      # println("Variable of the inner loop i = $i")
      # The patch for x_i ≠ 0 with i < j
      X = get_patches(C)[i+1]
      # println( "first patch" )
      #@show X
      R = ambient_ring( Y ) 
      V = localize( Y, gens(R)[i+1] )
      S = ambient_ring( X ) 
      U = localize( X, gens(S)[j] )
      R_loc = ambient_ring( V )
      S_loc = ambient_ring( U )
      R_vars = get_root_variables(V)
      S_vars = get_root_variables(U)
      R_unit = get_unit(V)
      S_unit = get_unit(U)
      Φ = AlgebraHomomorphism( R_loc, S_loc, 
			      vcat( [ S_vars[k]*S_unit for k in (1:i) ], # the first i-1 variables have the same indices
				   [S_unit], # this is inverted to pass to the other patch
				   [ S_vars[k-1]*S_unit for k in (i+2:j) ], # fill up with intermediate variables
				   [ S_vars[k]*S_unit for k in (j+1:n) ], # skip the j-th variable
				   [ S_vars[j] ] ) # map the unit to the appropriate coordinate
			      )

      f = AffSchMorphism( U, V, Φ )
      #@show f
      glueing = Glueing( inclusion_in_parent(U), inclusion_in_parent(V), f )
      #@show glueing
      push!( G, glueing )
      #@show G
      #@show typeof( G )
    end
    add_patch( C, Y, G )
  end
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
# rings of the affine patches. 
#
# We refrain from implementing this explicitly for 
# affine schemes since there ideal sheaves are just 
# in 1:1 correspondence with actual ideals.
#
mutable struct IdealSheaf <: AbstractCoherentSheaf 
  parent::CoveredScheme
  covering::Covering
  ideals::Vector{MPolyIdeal}
end

parent( I::IdealSheaf ) = I.parent
covering( I::IdealSheaf ) = I.covering
ideals( I::IdealSheaf ) = I.ideals

# Take a homogeneous ideal I in a graded homogeneous 
# ring S and turn it into an ideal sheaf on the 
# corresponding projective space
function to_ideal_sheaf( I::MPolyIdeal )
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

  function LineBundle( parent::CoveredScheme, covering::Covering, transitions::Array{MPolyElem,2} )
    L = new()
    L.parent = parent
    L.covering = covering
    L.transitions = transitions
    # TODO: Implement cheap plausibility checks.
    return L
  end
end

parent( L::LineBundle ) = L.parent
covering( L::LineBundle ) = L.covering
transitions( L::LineBundle ) = L.transitions

########################################################
# The tautological bundle on projective space
function OO( k::Ring, n::Int, d::Int )
  IPn = projective_space( k, n )
  C = get_coverings( IPn )[1]
  patches = get_patches(C)
  glueings = get_glueings(C)
  transitions = Array{MPolyElem,2}(undef,n+1,n+1)
  if d >= 0
    for i in OrderedMultiindex(n+1,2)
      a = get_unit(overlap(glueings[getindex(i,1),getindex(i,2)]))
      transitions[getindex(i,1), getindex(i,2)] = a^d
      #transitions[getindex(i,1), getindex(i,2)] = (get_unit(overlap(glueings[i])))^d
      transitions[i[2],i[1]] = (get_denom(overlap(glueings[i[1],i[2]])))^d
    end
  else
    d = -d
    for i in OrderedMultiindex(n+1,2)
      a = get_unit(overlap(glueings[getindex(i,1),getindex(i,2)]))
      transitions[getindex(i,2), getindex(i,1)] = a^d
      #transitions[getindex(i,1), getindex(i,2)] = (get_unit(overlap(glueings[i])))^d
      transitions[i[1],i[2]] = (get_denom(overlap(glueings[getindex(i,1),getindex(i,2)])))^d
    end
  end
  L = LineBundle( IPn, C, transitions )
  return L
end

OO( n::Int, d::Int ) = OO( QQ, n, d )

function Base.show( io::Base.IO, L::LineBundle )
  Base.print( io, "Line bundle on " )
  Base.print( io, L.parent )
  #Base.print( io, " with respect to the covering " )
  #Base.print( io, L.covering )
  Base.print( io, " given by the transition functions\n" )
  n = length( get_patches(L.covering) )
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

parent( L::VectorBundle ) = L.parent
covering( L::VectorBundle ) = L.covering
transitions( L::VectorBundle ) = L.transitions


#########################################################
# pullback of a line bundle along a refinement
function pullback( L::LineBundle, ρ::Refinement )
  if covering(L) != original_cover(ρ) 
    error( "The line bundle is defined on a different cover than the one being refined" )
  end
  m = length( get_patches( get_covering( L )))
  n = length( get_patches( get_original_cover( ρ )))
  inclusions = get_inclusions( ρ )
  new_transitions = Array{MPolyElem,2}(undef,n,n)
  glueings_original = get_glueings(domain(ρ))
  glueings_refined = get_glueings(codomain(ρ))
  for i in (1:n-1)
    for j in (i+1:n)
      if i != j
	f_i = inclusions[i] # the inclusion of X_i into Y_{n(i)}
	ϕ_i = pullback(f_i) # the associated pullback on the level of rings
	S = codomain( ϕ_i )
	S_loc = ambient_ring(overlap(glueings_refined[i,j])) 
	f_j = inclusions[j] # the inclusion of X_j into Y_{n(j)}
	ϕ_j = pullback(f_j) # the associated pullback on the level of rings
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
	  R = domain(ϕ) 
	  R_loc = ambient_ring(overlap(glueings_original[n_i,n_j]))
	  # ϕ_loc = 
	else
	  # In this last case, one needs to explicitly invert one 
	  # of the glueing isomorphisms.
	end 

      end
    end
  end
  M = LineBundle( parent(L), get_refined_cover(ρ), new_transitions )
  return M
end


# Todo adapt the code below to the new data structure

#=
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

=#

