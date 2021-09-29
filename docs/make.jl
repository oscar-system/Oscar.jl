push!(LOAD_PATH,"../src/")
using Documenter, JToric

makedocs(sitename="JToric -- Toric Geometry in Julia")

deploydocs(
   repo   = "github.com/oscar-system/JToric.jl.git",
   deps = nothing,
   target = "build",
   push_preview = true,
   make = nothing
)
