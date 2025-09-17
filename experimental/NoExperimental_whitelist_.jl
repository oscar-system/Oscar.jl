# List of all experimental packages that are currently needed for the CI job to succeed.
# Eventually, this list should be empty.
whitelist = String[
  "DoubleAndHyperComplexes",    # `MethodError: no method matching simplify(::SubquoModule{MPolyQuoRingElem{MPolyDecRingElem{FqMPolyRingElem, AbstractAlgebra.Generic.MPoly{FqMPolyRingElem}}}})`
  "ExteriorAlgebra",            # `test/Modules/PBWModules.jl` calls `exterior_algebra(::Field, ::Int)`
  "GModule",                    # `MethodError: no method matching (::FinGenAbGroup)(::FinGenAbGroupElem)`
  "ModStd",                     # `MethodError: no method matching monomial(::QQMPolyRing, ::Vector{Int64})` and many similar errors
  "Schemes",                    # TODO: untangle src/AlgebraicGeometry/Schemes/ and experimental/Schemes/
]
