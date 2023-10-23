# When a specific branch is loaded via `]add Package#branch` julia will only
# create a checkout and keep a bare git repo in a separate directory.
# In a bare repo HEAD will not point to the correct commit so we use the git
# tree-hash that Pkg.jl provides and manually map this to a corresponding
# commit.
function _lookup_commit_from_cache(url::AbstractString, tree::AbstractString)
   if Sys.which("git") !== nothing
      try
         path = Pkg.Types.add_repo_cache_path(url)
         if isdir(path)
            commit = readchomp(`sh -c "git -C $path log --oneline --all --pretty='tree %T;%H' | grep \"^tree $tree\" | cut -d\; -f2 | head -n1"`)
            return readchomp(`git -C $path show -s --format=", %h -- %ci" $commit`)
         end
      catch
      end
   end
   return ""
end

function _lookup_git_branch(dir::AbstractString; commit=false)
   info = ""
   if Sys.which("git") !== nothing &&
         isdir(joinpath(dir,".git"))
      try
         ref = readchomp(`git -C $dir rev-parse --abbrev-ref HEAD`)
         info = " - #$(ref)"
         if commit
            c = readchomp(`git -C $dir show -s --format="%h -- %ci" HEAD`)
            info = "$info, $c"
         end
      catch
      end
   end
   return info
end

function _deps_git_info(dep::Pkg.API.PackageInfo; commit=false)
   if dep.is_tracking_repo
      info = commit ? _lookup_commit_from_cache(dep.git_source, dep.tree_hash) : ""
      return " - #$(dep.git_revision)$info"
   elseif dep.is_tracking_path
      return _lookup_git_branch(dep.source; commit=commit)
   end
   return ""
end

function _print_dependency_versions(io::IO, deps::AbstractArray{<:AbstractString}; padding="    ", suffix="", branch=false, commit=false)
   width = maximum(length.(deps))+length(suffix)+2
   deps = filter(d->d.name in deps, collect(values(Pkg.dependencies())))
   deps = sort!(deps; by=x->x.name)
   for dep in deps
      print(io, "$(padding)$(rpad(dep.name*suffix, width, ' ')) v$(dep.version)")
      println(io, branch ? _deps_git_info(dep; commit=commit) : "")
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
   println(io, branch ? _lookup_git_branch(Oscar.oscardir; commit=commit) : "")
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

