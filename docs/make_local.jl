using Oscar
bla = normpath(joinpath(dirname(pathof(Oscar)), "..", "docs", "make_work.jl"))
local_build = true
include(bla)

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
