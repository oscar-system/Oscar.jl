include("Types.jl")
include("CoveredScheme.jl")
include("FunctionFields.jl")
include("ProjectiveModules.jl")
include("SpaceGerms.jl")
include("Sheaves.jl")
include("IdealSheaves.jl")
include("AlgebraicCycles.jl")
include("WeilDivisor.jl")
include("CoveredProjectiveSchemes.jl")
include("StructureSheaf.jl")

include("SimplifiedAffineScheme.jl")
include("CoherentSheaves.jl")
include("LazyGluing.jl")
include("CartierDivisor.jl")
include("Auxiliary.jl")
include("BlowupMorphism.jl")
include("duValSing.jl")
include("elliptic_surface.jl")
include("MorphismFromRationalFunctions.jl")

include("ToricIdealSheaves/auxiliary.jl")
include("ToricIdealSheaves/constructors.jl")
include("ToricIdealSheaves/attributes.jl")
include("ToricIdealSheaves/methods.jl")

include("ToricDivisors/constructors.jl")
include("ToricDivisors/attributes.jl")

include("NormalToricVarieties/attributes.jl")

include("ToricBlowups/types.jl")
include("ToricBlowups/constructors.jl")
include("ToricBlowups/attributes.jl")
include("ToricBlowups/methods.jl")

include("DerivedPushforward.jl")

# Deprecated after 0.15
Base.@deprecate_binding _compute_inherited_glueing _compute_inherited_gluing
Base.@deprecate_binding base_glueing base_gluing
Base.@deprecate_binding inherit_glueings! inherit_gluings!
Base.@deprecate_binding AbsProjectiveGlueing AbsProjectiveGluing
Base.@deprecate_binding CoveredProjectiveGlueingData CoveredProjectiveGluingData
Base.@deprecate_binding InheritGlueingData InheritGluingData
Base.@deprecate_binding LazyProjectiveGlueing LazyProjectiveGluing
Base.@deprecate_binding ProjectiveGlueing ProjectiveGluing
Base.@deprecate_binding ProjectiveGlueingData ProjectiveGluingData
