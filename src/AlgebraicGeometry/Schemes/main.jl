import AbstractAlgebra.Ring
import Base: intersect

########################################################################
# Type declarations                                                    #
#                                                                      #
# These come first as the types will be needed for the internals of    #
# what follows                                                         #
########################################################################
include("Types.jl")
include("Methods.jl")

include("AffineSchemes/Objects/Types.jl")
include("AffineSchemes/Morphisms/Types.jl")
include("PrincipalOpenSubset/Objects/Types.jl")
include("SpecOpen/Objects/Types.jl")
include("SpecOpen/Rings/Types.jl")
include("SpecOpen/Morphisms/Types.jl")
include("ClosedEmbedding/Types.jl") # Needs SpecOpen for their complements
include("Glueing/Types.jl")
include("Covering/Objects/Types.jl")
include("Covering/Morphisms/Types.jl")
include("CoveredSchemes/Objects/Types.jl")
include("CoveredSchemes/Morphisms/Types.jl")
include("AbstractTypes.jl")
include("AffineVariety/Objects/Types.jl")
include("AffineAlgebraicSet/Objects/Types.jl")
include("ProjectiveSchemes/Objects/Types.jl")
include("ProjectiveSchemes/Morphisms/Types.jl")
include("ProjectiveAlgebraicSet/Objects/Types.jl")
include("ProjectiveVariety/Objects/Types.jl")

########################################################################
# Affine schemes                                                       #
########################################################################
include("AffineSchemes/Objects/Constructors.jl")
include("AffineSchemes/Objects/Properties.jl")
include("AffineSchemes/Objects/Attributes.jl")
include("AffineSchemes/Objects/Methods.jl")

include("AffineSchemes/Morphisms/Constructors.jl")
include("AffineSchemes/Morphisms/Properties.jl")
include("AffineSchemes/Morphisms/Attributes.jl")
include("AffineSchemes/Morphisms/Methods.jl")

########################################################################
# Projective Schemes                                                   #
########################################################################
include("ProjectiveSchemes/Objects/Constructors.jl")
include("ProjectiveSchemes/Objects/Properties.jl")
include("ProjectiveSchemes/Objects/Attributes.jl")
include("ProjectiveSchemes/Objects/Methods.jl")

include("ProjectiveSchemes/Morphisms/Constructors.jl")
include("ProjectiveSchemes/Morphisms/Attributes.jl")
include("ProjectiveSchemes/Morphisms/Methods.jl")

########################################################################
# Principal open subsets of affine schemes                             #
########################################################################
include("PrincipalOpenSubset/Objects/Constructors.jl")
include("PrincipalOpenSubset/Objects/Properties.jl")
include("PrincipalOpenSubset/Objects/Attributes.jl")
include("PrincipalOpenSubset/Objects/Methods.jl")

########################################################################
# Open subsets of affine schemes                                       #
########################################################################
include("SpecOpen/Objects/Constructors.jl")
include("SpecOpen/Objects/Properties.jl")
include("SpecOpen/Objects/Attributes.jl")
include("SpecOpen/Objects/Methods.jl")

include("SpecOpen/Rings/Constructors.jl")
include("SpecOpen/Rings/Properties.jl")
include("SpecOpen/Rings/Attributes.jl")
include("SpecOpen/Rings/Methods.jl")

include("SpecOpen/Morphisms/Constructors.jl")
include("SpecOpen/Morphisms/Attributes.jl")
include("SpecOpen/Morphisms/Methods.jl")

########################################################################
# ClosedEmbeddings of affine schemes                                   #
########################################################################

include("ClosedEmbedding/Attributes.jl")

########################################################################
# Glueings of affine schemes                                           #
########################################################################
include("Glueing/Constructors.jl")
include("Glueing/Attributes.jl")
include("Glueing/Methods.jl")

########################################################################
# Coverings of schemes                                                 #
########################################################################
include("Covering/Objects/Constructors.jl")
include("Covering/Objects/Attributes.jl")
include("Covering/Objects/Methods.jl")

include("Covering/Morphisms/Constructors.jl")
include("Covering/Morphisms/Attributes.jl")
include("Covering/Morphisms/Methods.jl")

########################################################################
# Covered schemes                                                      #
########################################################################
include("CoveredSchemes/Objects/Constructors.jl")
include("CoveredSchemes/Objects/Properties.jl")
include("CoveredSchemes/Objects/Attributes.jl")
include("CoveredSchemes/Objects/Methods.jl")

include("CoveredSchemes/Morphisms/Constructors.jl")
include("CoveredSchemes/Morphisms/Attributes.jl")
include("CoveredSchemes/Morphisms/Methods.jl")

########################################################################
# Affine Algebraic Sets                                                #
########################################################################
include("AffineAlgebraicSet/Objects/Constructors.jl")
include("AffineAlgebraicSet/Objects/Properties.jl")
include("AffineAlgebraicSet/Objects/Attributes.jl")
include("AffineAlgebraicSet/Objects/Methods.jl")

########################################################################
# Affine Varietes
########################################################################
include("AffineVariety/Objects/Constructors.jl")
include("AffineVariety/Objects/Properties.jl")
include("AffineVariety/Objects/Attributes.jl")
include("AffineVariety/Objects/Methods.jl")

########################################################################
# Projective Algebraic Sets                                                #
########################################################################
include("ProjectiveAlgebraicSet/Objects/Constructors.jl")
include("ProjectiveAlgebraicSet/Objects/Properties.jl")
include("ProjectiveAlgebraicSet/Objects/Attributes.jl")
include("ProjectiveAlgebraicSet/Objects/Methods.jl")

########################################################################
# Projective Varietes
########################################################################
include("ProjectiveVariety/Objects/Constructors.jl")
include("ProjectiveVariety/Objects/Properties.jl")
include("ProjectiveVariety/Objects/Attributes.jl")
include("ProjectiveVariety/Objects/Methods.jl")

# A global variable to throw an error whenever an enabled internal test was found.
#
# Synopsis: Many methods in Oscar have an additional optional argument `check::Bool=true`
# which is set to `true` by default and which forces several internal checks to be run 
# to check for sane input. 
#
# An example is the map from a quotient ring f : A -> R to another ring R: For f to be 
# well defined, one has to check that the lift F : P -> R of this map to the polynomial 
# ring P with A = P/I takes all elements in I to zero. In the default case, the constructor 
# for f will check that the given map is well defined, but one does not want this check 
# to be run when calling the constructor internally on verified sane input. 
#
# To this end, one should do the function calls with `check=false` internally. However, 
# it is very likely that one forgets to include this statement in all calls when programming. 
# To this end, we introduce a global flag here which does the following. 
#
# Every function using the optional `check::Bool=true` argument should have the line 
#
#   check && CHECK_ERROR && error("internal checks were enabled")
#
# Then the developer can set the flag below to `true` and run their test code. This will 
# throw the error whenever an unexpected internal check would have been triggered together 
# with a subsequent stacktrace. Then the internal checks can be switched off. 

CHECK_ERROR = false

