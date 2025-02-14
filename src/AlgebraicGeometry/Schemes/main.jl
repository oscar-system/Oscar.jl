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
include("AffineSchemes/SimplifiedAffineScheme/Types.jl")
include("PrincipalOpenSubset/Objects/Types.jl")
include("PrincipalOpenInclusion/Types.jl")
include("AffineSchemeOpenSubscheme/Objects/Types.jl")
include("AffineSchemeOpenSubscheme/Rings/Types.jl")
include("AffineSchemeOpenSubscheme/Morphisms/Types.jl")
include("ClosedEmbedding/Types.jl") # Needs AffineSchemeOpenSubscheme for their complements
include("Gluing/Types.jl")
include("Gluing/LazyGluing/Types.jl")
include("Covering/Objects/Types.jl")
include("Covering/Morphisms/Types.jl")
include("CoveredSchemes/Objects/Types.jl")
include("CoveredSchemes/Morphisms/Types.jl")
include("CoveredSchemes/Morphisms/CompositeCoveredSchemeMorphism/Types.jl")
include("AbstractTypes.jl")
include("AffineVariety/Objects/Types.jl")
include("AffineAlgebraicSet/Objects/Types.jl")
include("ProjectiveSchemes/Objects/Types.jl")
include("ProjectiveSchemes/Morphisms/Types.jl")
include("ProjectiveAlgebraicSet/Objects/Types.jl")
include("ProjectiveVariety/Objects/Types.jl")
include("Sheaves/Types.jl")
include("CoveredSchemes/Morphisms/CoveredClosedEmbedding/Types.jl") #needs IdealSheaf
include("Sheaves/PushforwardIdealSheaf.jl")
include("FunctionField/Types.jl")
include("Divisors/Types.jl")
include("CoveredSchemes/Morphisms/MorphismFromRationalFunctions/Types.jl") 


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
include("AffineSchemes/SimplifiedAffineScheme/Methods.jl")

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

include("PrincipalOpenInclusion/Constructors.jl")
include("PrincipalOpenInclusion/Methods.jl")
include("PrincipalOpenInclusion/Attributes.jl")

########################################################################
# Open subsets of affine schemes                                       #
########################################################################
include("AffineSchemeOpenSubscheme/Objects/Constructors.jl")
include("AffineSchemeOpenSubscheme/Objects/Properties.jl")
include("AffineSchemeOpenSubscheme/Objects/Attributes.jl")
include("AffineSchemeOpenSubscheme/Objects/Methods.jl")

include("AffineSchemeOpenSubscheme/Rings/Constructors.jl")
include("AffineSchemeOpenSubscheme/Rings/Properties.jl")
include("AffineSchemeOpenSubscheme/Rings/Attributes.jl")
include("AffineSchemeOpenSubscheme/Rings/Methods.jl")

include("AffineSchemeOpenSubscheme/Morphisms/Constructors.jl")
include("AffineSchemeOpenSubscheme/Morphisms/Attributes.jl")
include("AffineSchemeOpenSubscheme/Morphisms/Methods.jl")

########################################################################
# ClosedEmbeddings of affine schemes                                   #
########################################################################

include("ClosedEmbedding/Attributes.jl")

########################################################################
# Gluings of affine schemes                                           #
########################################################################
include("Gluing/Constructors.jl")
include("Gluing/Attributes.jl")
include("Gluing/Methods.jl")
include("Gluing/LazyGluing/Methods.jl")

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
include("CoveredSchemes/Morphisms/Properties.jl")
include("CoveredSchemes/Morphisms/CoveredClosedEmbedding/Methods.jl") 
include("CoveredSchemes/Morphisms/MorphismFromRationalFunctions/Methods.jl") 
include("CoveredSchemes/Morphisms/CompositeCoveredSchemeMorphism/Methods.jl")


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

########################################################################
# Sheaves
########################################################################
include("Sheaves/Sheaves.jl")
include("Sheaves/CoherentSheaves.jl")
include("Sheaves/StructureSheaf.jl")
include("Sheaves/IdealSheaves.jl")
include("Sheaves/Methods.jl")

########################################################################
# Rational functions
########################################################################
include("FunctionField/FunctionFields.jl")

########################################################################
# Divisors
########################################################################
include("Divisors/AlgebraicCycles.jl")
include("Divisors/WeilDivisor.jl")
include("Divisors/CartierDivisor.jl")
include("Divisors/base_change.jl")

########################################################################
# Blowups
########################################################################
include("CoveredProjectiveScheme/Types.jl")
include("CoveredProjectiveScheme/CoveredProjectiveScheme.jl")
include("BlowupMorphism/Types.jl")
include("BlowupMorphism/BlowupMorphism.jl")
