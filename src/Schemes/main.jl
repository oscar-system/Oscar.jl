import AbstractAlgebra.Ring
import Base: intersect

########################################################################
# Type declarations                                                    #
#                                                                      #
# These come first as the types will be needed for the internals of    #
# what follows                                                         #
########################################################################
include("Types.jl")

include("AffineSchemes/Objects/Types.jl")
include("AffineSchemes/Morphisms/Types.jl")
include("PrincipalOpenSubset/Objects/Types.jl")
include("SpecOpen/Objects/Types.jl")
include("SpecOpen/Rings/Types.jl")
include("SpecOpen/Morphisms/Types.jl")

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

