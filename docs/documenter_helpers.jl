import DelimitedFiles

const docstats = Dict{String,NamedTuple}()
if haskey(ENV, "GITHUB_ACTION") || haskey(ENV, "OSCAR_TEST_STATS")
  Base.atexit() do
    open("test-stats-doctests.csv", "a") do io
      @static if VERSION > v"1.11.0"
        println(io, "path,time,ctime,rctime,gctime,alloc")
        DelimitedFiles.writedlm(io, ((k, v.time, v.ctime, v.rctime, v.gctime, v.alloc) for (k,v) in docstats), ",")
      else
        println(io, "path,time,gctime,alloc")
        DelimitedFiles.writedlm(io, ((k, v.time, v.gctime, v.alloc) for (k,v) in docstats), ",")
      end
    end
  end
end


# this function is slightly more specific than the one from documenter, print the corresponding code location
# calls the original function via invoke and prints a timing for the doctest
function Documenter.eval_repl(block::Documenter.MarkdownAST.CodeBlock, sandbox::Module, meta::Dict, doc::Documenter.Document, page)
  src_lines = Documenter.find_block_in_file(block.code, meta[:CurrentFile])
  # skip stats if there was a failure
  if length(doc.internal.errors) > 0
    invoke(Documenter.eval_repl, Tuple{Documenter.MarkdownAST.CodeBlock,Any,Dict,Documenter.Document,Any}, block, sandbox, meta, doc, page)
  else
    loc = Documenter.locrepr(meta[:CurrentFile], src_lines)
    println("page: $loc")
    stats = @timed @time invoke(Documenter.eval_repl, Tuple{Documenter.MarkdownAST.CodeBlock,Any,Dict,Documenter.Document,Any}, block, sandbox, meta, doc, page)
    @static if VERSION > v"1.11.0"
      docstats[loc] = (time=stats.time, ctime=stats.compile_time-stats.recompile_time, rctime=stats.recompile_time, gctime=stats.gctime, alloc=stats.bytes/2^30)
    else
      docstats[loc] = (time=stats.time, gctime=stats.gctime, alloc=stats.bytes/2^30)
    end
  end
end
