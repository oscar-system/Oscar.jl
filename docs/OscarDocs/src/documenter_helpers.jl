# Doctest timings, keyed by source location. Filled by our method of
# `Documenter.eval_repl`, written out at exit by `setup_doctest_stats`.
const docstats = Dict{String,NamedTuple}()

const stats_registered = Ref(false)

# Statistics are only collected in CI or when explicitly requested.
function setup_doctest_stats(Oscar::Module)
  (haskey(ENV, "GITHUB_ACTION") || haskey(ENV, "OSCAR_TEST_STATS")) || return nothing
  stats_registered[] && return nothing
  stats_registered[] = true

  metadata = Oscar._test_stats_metadata()
  platform = Sys.islinux() ? "linux" : "macos"
  juliaVersion = join(split("$VERSION", ".")[1:2], ".")
  filename = "test-stats_$(metadata.timestamp)_$(platform)_$(juliaVersion)_doctests_$(metadata.commit).csv"

  Base.atexit() do
    open(filename, "a") do io
      @static if VERSION > v"1.11.0"
        println(io, "path,time,ctime,rctime,gctime,alloc")
        for (k, v) in docstats
          println(io, join((k, v.time, v.ctime, v.rctime, v.gctime, v.alloc), ","))
        end
      else
        println(io, "path,time,gctime,alloc")
        for (k, v) in docstats
          println(io, join((k, v.time, v.gctime, v.alloc), ","))
        end
      end
    end
  end
  return nothing
end

# this function is slightly more specific than the one from documenter, print the corresponding code location
# calls the original function via invoke and prints a timing for the doctest
function Documenter.eval_repl(block::Documenter.MarkdownAST.CodeBlock, sandbox::Module, meta::Dict, doc::Documenter.Document, page; kwargs...)
  src_lines = Documenter.find_block_in_file(block.code, meta[:CurrentFile])
  # skip stats if there was a failure
  if length(doc.internal.errors) > 0
    invoke(Documenter.eval_repl, Tuple{Documenter.MarkdownAST.CodeBlock,Any,Dict,Documenter.Document,Any}, block, sandbox, meta, doc, page; kwargs...)
  else
    loc = Documenter.locrepr(meta[:CurrentFile], src_lines)
    println("page: $loc")
    stats = @timed @time invoke(Documenter.eval_repl, Tuple{Documenter.MarkdownAST.CodeBlock,Any,Dict,Documenter.Document,Any}, block, sandbox, meta, doc, page; kwargs...)
    @static if VERSION > v"1.11.0"
      docstats[loc] = (time=stats.time, ctime=stats.compile_time-stats.recompile_time, rctime=stats.recompile_time, gctime=stats.gctime, alloc=stats.bytes/2^30)
    else
      docstats[loc] = (time=stats.time, gctime=stats.gctime, alloc=stats.bytes/2^30)
    end
  end
end
