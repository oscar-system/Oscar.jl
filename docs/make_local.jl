using Oscar
bla = normpath(joinpath(dirname(pathof(Oscar)), "..", "docs", "make_work.jl"))
local_build = true
include(bla)

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
   push_preview = true,
#  make = () -> run(`mkdocs build`),
   make = nothing
)


function open(filename; wait = false)
    @static if Sys.isapple()
        run(`open $(filename)`; wait = wait)
    elseif Sys.islinux() || Sys.isbsd()
        run(`xdg-open $(filename)`; wait = wait)
    elseif Sys.iswindows()
        cmd = get(ENV, "COMSPEC", "cmd.exe")
        run(`$(cmd) /c start $(filename)`; wait = wait)
    else
        @warn("Opening files the default application is not supported on this OS.",
              KERNEL = Sys.KERNEL)
    end
end

bla = normpath(joinpath(dirname(pathof(Oscar)), "..", "docs", "build", "index.html"))
open(bla)
