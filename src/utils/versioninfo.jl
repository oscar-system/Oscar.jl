# When a specific branch is loaded via `]add Package#branch` julia will only
# create a checkout and keep a bare git repo in a separate directory.
# In a bare repo HEAD will not point to the correct commit so we use the git
# tree-hash that Pkg.jl provides and manually map this to a corresponding
# commit.
function _lookup_commit_from_cache!(info::Dict, url::AbstractString, tree::AbstractString)
  if Sys.which("git") !== nothing
    try
      path = Pkg.Types.add_repo_cache_path(url)
      if isdir(path)
        commit = readchomp(`sh -c "git -C $path log --oneline --all --pretty='tree %T;%H' | grep \"^tree $tree\" | cut -d\; -f2 | head -n1"`)
        c = readchomp(`git -C $path show -s --format="%H#%ci" $commit`)
        (info[:commit], info[:date]) = split(c, "#")
      end
    catch
    end
  end
end

function _lookup_git_branch!(info::Dict, dir::AbstractString)
  # the .git entry might be a file instead of a dir when using git worktrees
  if Sys.which("git") !== nothing &&
    ispath(joinpath(dir, ".git"))
    try
      ref = readchomp(`git -C $dir rev-parse --abbrev-ref HEAD`)
      info[:branch] = ref
      c = readchomp(`git -C $dir show -s --format="%H#%ci" HEAD`)
      (info[:commit], info[:date]) = split(c, "#")
    catch
    end
  end
end

function _get_git_info(dep::Union{Pkg.API.PackageInfo,AbstractString})
  info = Dict{Symbol,String}()
  if dep isa Pkg.API.PackageInfo && dep.is_tracking_repo
    _lookup_commit_from_cache!(info, dep.git_source, dep.tree_hash)
    # this might be a branch, tag, or commit hash
    info[:branch] = dep.git_revision
  elseif dep isa Pkg.API.PackageInfo && dep.is_tracking_path
    _lookup_git_branch!(info, dep.source)
  elseif dep isa AbstractString
    _lookup_git_branch!(info, dep)
  end
  return info
end

function _get_oscar_git_info()
  # Oscar is either one of the dependencies or the active project.
  # For the active project we try to use the Oscar path as git directory.
  oscarinfo = get(Pkg.dependencies(), PROJECT_UUID, Oscar.oscardir)
  return _get_git_info(oscarinfo)
end

function _format_git_info(info::Dict; branch=true, commit=false)
  val = String[]
  if branch && haskey(info, :branch)
    push!(val, "#$(info[:branch])")
  end
  if commit && haskey(info, :commit)
    push!(val, "$(info[:commit][1:10]) -- $(info[:date])")
  end
  return length(val) > 0 ? " - $(join(val, ", "))" : ""
end

function _print_dependency_versions(io::IO, deps::AbstractArray{<:AbstractString}; padding="    ", suffix="", branch=false, commit=false)
  width = maximum(length.(deps))+length(suffix)+2
  deps = filter(d->d.name in deps, collect(values(Pkg.dependencies())))
  deps = sort!(deps; by=x->x.name)
  for dep in deps
    print(io, "$(padding)$(rpad(dep.name*suffix, width, ' ')) v$(dep.version)")
    if branch || commit
      print(io, _format_git_info(_get_git_info(dep); branch=branch, commit=commit))
    end
    println(io)
  end
end

@doc raw"""
    Oscar.versioninfo(io::IO=stdout; branch=false, jll=false, julia=false, commit=false, full=false)

Print the versions of all Oscar-related dependencies.

# Arguments
- `branch::Bool=false`: include git branch name in the output
- `commit::Bool=false`: include git commit hash and date where applicable
- `jll::Bool=false`   : include binary packages (jll) in the output
- `julia::Bool=false` : include julia `versioninfo` output
- `full::Bool=false`  : include all of the above
"""
function versioninfo(io::IO=stdout; branch=false, jll=false, julia=false, commit=false, full=false)
  if full
    branch = jll = julia = commit = true
  end
  print(io, "OSCAR version $(VERSION_NUMBER)")
  if branch || commit
    print(io, _format_git_info(_get_oscar_git_info(); branch=branch, commit=commit))
  end
  println(io)
  println(io, "  combining:")
  _print_dependency_versions(io, cornerstones; suffix=".jl", branch=branch, commit=commit)
  if jll
    println(io, "  building on:")
    _print_dependency_versions(io, jll_deps; branch=branch, commit=commit)
    println(io, "See `]st -m` for a full list of dependencies.")
  end
  if julia
    println(io, "")
    Main.InteractiveUtils.versioninfo(io)
    println(io, Base.TAGGED_RELEASE_BANNER)
  end
end

