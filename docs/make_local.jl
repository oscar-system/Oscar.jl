module BuildDoc
using Oscar

include(normpath(joinpath(Oscar.oscardir, "docs", "make_work.jl")))

function open_doc()
    filename = normpath(joinpath(dirname(pathof(Oscar)), "..", "docs", "build", "index.html"))
    @static if Sys.isapple()
        run(`open $(filename)`; wait = false)
    elseif Sys.islinux() || Sys.isbsd()
        run(`xdg-open $(filename)`; wait = false)
    elseif Sys.iswindows()
        cmd = get(ENV, "COMSPEC", "cmd.exe")
        run(`$(cmd) /c start $(filename)`; wait = false)
    else
        @warn("Opening files the default application is not supported on this OS.",
              KERNEL = Sys.KERNEL)
    end
end

end

using .BuildDoc
