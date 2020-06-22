using Documenter, Oscar, DocumenterMarkdown

makedocs(
         format   = Documenter.HTML(),
#         format   = Markdown(),
         sitename = "Oscar.jl",
         modules = [Oscar],
         clean = true,
         doctest = false,
         pages    = [
             "index.md",
             "Rings" => [ "Rings/integer.md",
			  "Rings/rational.md"],
             "Groups" => [ "Groups/groups.md" ]
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
