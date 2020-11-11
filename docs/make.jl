using Documenter, Oscar, DocumenterMarkdown
using AbstractAlgebra

AA_dir = pkgdir(AbstractAlgebra)
Oscar_dir = pkgdir(Oscar)

cp(joinpath(AA_dir, "docs", "src"), joinpath(Oscar_dir, "docs", "src", "AAsrc"), force=true)

# TODO: load this directly from AbstractAlgebra.jl
AA_pages = ["index.md",
           "constructors.md",
           "Rings" => [ "rings.md",
                        "ncrings.md",
                        "euclidean.md",
                        "integer.md",
                        "polynomial_rings.md",
                        "polynomial.md",
                        "ncpolynomial.md",
                        "mpolynomial_rings.md",
                        "mpolynomial.md",
                        "series_rings.md",
                        "series.md",
                        "puiseux.md",
                        "residue_rings.md",
                        "residue.md"],
           "Fields" => [ "fields.md",
                         "fraction_fields.md",
                         "fraction.md",
                         "rational.md",
                         "finfield.md",
                         "real.md",
                         "numberfield.md"],
           "Groups" => [ "perm.md",
                         "ytabs.md"],
           "Modules" => [ "module.md",
                          "free_module.md",
                          "submodule.md",
                          "quotient_module.md",
                          "direct_sum.md",
                          "module_homomorphism.md"],
           "Matrices" => [ "matrix_spaces.md",
                           "matrix.md",
                           "matrix_algebras.md"],
           "Maps" => [ "map.md",
                       "functional_map.md",
                       "map_cache.md",
                       "map_with_inverse.md"],
           "types.md" # Appendix A
           ]

# add "AAsrc" prefix to AA files
replace!(AA_pages) do s
    if s isa String
        "AAsrc/" * s
    else
        replace!(s.second) do t
            "AAsrc/" * t
        end
        s
    end
end

makedocs(
         format   = Documenter.HTML(),
#         format   = Markdown(),
         sitename = "Oscar.jl",
         push_preview = true,
         modules = [Oscar, AbstractAlgebra],
         clean = true,
         doctest = false,
         pages = ["index.md",
                  "Rings" => [ "Rings/integer.md",
			       "Rings/rational.md"],
                  "Groups" => [ "Groups/groups.md" ],
                  "AbstractAlgebra.jl" => AA_pages
                  ]
)

#deploydocs(
#   julia = "1.3",
#   repo   = "github.com/oscar-system/Oscar.jl.git",
#   target = "build",
#   deps = nothing,
#   make   = nothing,
#   osname = "linux"
#)

deploydocs(
   repo   = "github.com/oscar-system/Oscar.jl.git",
#  deps = Deps.pip("pymdown-extensions", "pygments", "mkdocs", "python-markdown-math", "mkdocs-material", "mkdocs-cinder"),
   deps = nothing,
   target = "build",
#  make = () -> run(`mkdocs build`),
   make = nothing
)
