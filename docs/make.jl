using Documenter, Oscar

include(normpath(joinpath(Oscar.oscardir, "docs", "make_work.jl")))

BuildDoc.doit(Oscar; warnonly=false, local_build=false, doctest=false)

deploydocs(
   repo   = "github.com/oscar-system/Oscar.jl.git",
#  deps = Deps.pip("pymdown-extensions", "pygments", "mkdocs", "python-markdown-math", "mkdocs-material", "mkdocs-cinder"),
   deps = nothing,
   target = "build",
   push_preview = true,
#  make = () -> run(`mkdocs build`),
   make = nothing,
   forcepush = true
)
