#
# Looking at the manual once it has been built.
#

function open_doc(Oscar::Module)
  filename = normpath(Base.pkgdir(Oscar), "docs", "build", "index.html")
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

# Served from a separate process so that the REPL stays usable, and from its own
# temporary environment so that LiveServer stays out of the docs environment.
function start_doc_preview_server(Oscar::Module; open_browser::Bool = true, port::Int = 8000)
  build_dir = normpath(Base.pkgdir(Oscar), "docs", "build")
  cmd = """
        using Pkg;
        Pkg.activate(temp = true);
        Pkg.add("LiveServer");
        using LiveServer;
        LiveServer.serve(dir = "$build_dir", launch_browser = $open_browser, port = $port);
        """
  live_server_process = run(`$(Base.julia_cmd()) -e $cmd`, wait = false)
  atexit(_ -> kill(live_server_process))
  @info "Starting server with PID $(getpid(live_server_process)) listening on 127.0.0.1:$port"
  return nothing
end
