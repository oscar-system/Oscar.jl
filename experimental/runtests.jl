for pkg in Oscar.exppkgs
  include("$pkg/test/runtests.jl")
end
