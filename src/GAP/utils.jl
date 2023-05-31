function _write_gap_file(relpath::String, content::String)
  gaproot = GAP.Setup.gaproot()
  path = joinpath(gaproot, relpath)
  mkpath(dirname(path))
  # write to a temporary file and move later to avoid race conditions
  tmpfile = tempname(; cleanup=false)
  open(tmpfile, "w") do file
    write(file, content)
  end
  Base.Filesystem.mv(tmpfile, path; force=true)
end
