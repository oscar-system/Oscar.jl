# We read experimental and filter out all packages that follow our desired
# scheme. Remember those packages to avoid doing this all over again for docs
# and test.
# We don't want to interfere with existing stuff in experimental though.
const expdir = joinpath(@__DIR__, "../experimental")
const oldexppkgs = [
  "ExteriorAlgebra",
  "GModule",
  "MatrixGroups",
  "ModStd",
  "Rings",
  "Schemes",
]
# DEVELOPER OPTION:
# If an experimental package A depends on another experimental package B, one
# can add `"B", "A"` to `orderedpkgs` below to ensure that B is loaded before A.
# DO NOT USE THIS UNLESS YOU KNOW THE CONSEQUENCES.
# For more background, see https://github.com/oscar-system/Oscar.jl/issues/2300.
const orderedpkgs = [
  "LieAlgebras",
  "BasisLieHighestWeight",   # nees code from LieAlgebras
]
exppkgs = filter(x->isdir(joinpath(expdir, x)) && !(x in oldexppkgs) && !(x in orderedpkgs), readdir(expdir))
append!(exppkgs, orderedpkgs)

# force trigger recompile when folder changes
include_dependency(".")

# for the NoExperimental CI job, remove all experimental packages from the inclusion lists
if haskey(ENV, "OSCAR_TEST_NOEXPERIMENTAL")
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
  filter!(in(whitelist), exppkgs)
  filter!(in(whitelist), oldexppkgs)
end

# Error if something is incomplete in experimental
for pkg in exppkgs
  if !isfile(joinpath(expdir, pkg, "src", "$pkg.jl"))
    error("experimental/$pkg is incomplete: $pkg/src/$pkg.jl missing. See the documentation at https://docs.oscar-system.org/dev/Experimental/intro/ for details.")
  end
  path = joinpath(expdir, pkg, "test")
  if !isdir(path) || length(filter(endswith(".jl"), readdir(path))) == 0
    error("experimental/$pkg is incomplete: $pkg/test/ missing or empty. See the documentation at https://docs.oscar-system.org/dev/Experimental/intro/ for details.")
  end
  # Load the package
  include("$pkg/src/$pkg.jl")
end

# Force some structure for `oldexppkgs`
for pkg in oldexppkgs
  if !isfile(joinpath(expdir, pkg, "$pkg.jl"))
    error("experimental/$pkg is incomplete: $pkg/$pkg.jl missing. Please fix this or remove $pkg from `oldexppkgs`.")
  end
  # Load the package
  include("$pkg/$pkg.jl")
end
