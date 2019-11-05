using Documenter, Oscar

makedocs(
         format   = :html,
         sitename = "Oscar.jl",
         modules = [Oscar],
         clean = true,
         doctest = false,
         pages    = [
             "index.md",
             "Rings" => [ "integer.md"]
         ]
)

deploydocs(
   julia = "1.0",
   repo   = "github.com/oscar-system/Oscar.jl.git",
   target = "build",
   deps = nothing,
   make   = nothing,
   osname = "linux"
)
