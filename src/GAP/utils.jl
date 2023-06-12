function _write_gap_file(relpath::String, content::String)
  gaproot = GAP.Setup.gaproot()
  path = joinpath(gaproot, relpath)
  mkpath(dirname(path))
  # write to a temporary file and move later to avoid race conditions
  tmpfile = tempname(dirname(path); cleanup=false)
  open(tmpfile, "w") do file
    write(file, content)
  end
  Base.Filesystem.rename(tmpfile, path)
end
