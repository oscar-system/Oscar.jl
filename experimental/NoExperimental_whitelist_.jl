# List of all experimental packages that are currently needed for the CI job to succeed.
# Eventually, this list should be empty.
whitelist = String[
  "DoubleAndHyperComplexes",    # `MethodError: no method matching simplify(::SubquoModule{MPolyQuoRingElem{MPolyDecRingElem{FqMPolyRingElem, AbstractAlgebra.Generic.MPoly{FqMPolyRingElem}}}})`
  "ExteriorAlgebra",            # `undefined binding 'exterior_algebra' in `@docs` block in src/NoncommutativeAlgebra/PBWAlgebras/quotients.md:40-42` and `Error During Test at /home/runner/work/Oscar.jl/Oscar.jl/test/Rings/PBWAlgebraQuo.jl:41`
  "GaloisGrp",                  # `no docs found for 'fixed_field(C::Oscar.GaloisGrp.GaloisCtx, s::Vector{PermGroup})' in `@docs` block in src/NumberTheory/galois.md:275-278`
  "GModule",                    # many doctest failures and `MethodError: no method matching numerator(::QQPolyRingElem)`
  "InvariantTheory",            # `undefined binding 'linearly_reductive_group' in `@docs` block in src/InvariantTheory/reductive_groups.md:53-55` and more docs errors`
  "ModStd",                     # `MethodError: no method matching monomial(::QQMPolyRing, ::Vector{Int64})` and many similar errors
  "Rings",                      # used by Rings/binomial_ideals.jl, see https://github.com/oscar-system/Oscar.jl/blob/13282dfd07b6aee58e433a45353f48261cda787b/src/Oscar.jl#L268
  "Schemes",                    # TODO: untangle src/AlgebraicGeometry/Schemes/ and experimental/Schemes/
]
