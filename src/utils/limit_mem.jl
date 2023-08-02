# try working around the CI issues with nightly + 1.10 (#2441) for now:
if haskey(ENV, "GITHUB_ACTIONS") && VERSION >= v"1.10.0-DEV"
  # 3.5GB might work on linux, at least sometimes
  # 10GB should be fine on macos
  maxmem = Sys.islinux() ? 3.5*(2^30) : 10*(2^30)
  println("OscarCI: Limiting memory to ", Base.format_bytes(maxmem));
  @ccall jl_gc_set_max_memory(maxmem::UInt64)::Cvoid
end

