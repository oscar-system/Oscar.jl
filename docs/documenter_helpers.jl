# this function is slightly more specific than the one from documenter, print the corresponding code location
# calls the original function via invoke and prints a timing for the doctest
if isdefined(Documenter, :DocTests)
  function Documenter.DocTests.eval_repl(block, sandbox::Module, meta::Dict, doc::Documenter.Documents.Document, page)
    src_lines = Documenter.Utilities.find_block_in_file(block.code, meta[:CurrentFile])
    # skip stats if there was a failure
    if length(doc.internal.errors) > 0
      invoke(Documenter.DocTests.eval_repl, Tuple{Any,Any,Dict,Documenter.Documents.Document,Any}, block, sandbox, meta, doc, page)
    else
      println("page: $(Documenter.Utilities.locrepr(meta[:CurrentFile], src_lines))")
      @time invoke(Documenter.DocTests.eval_repl, Tuple{Any,Any,Dict,Documenter.Documents.Document,Any}, block, sandbox, meta, doc, page)
    end
  end
else
  function Documenter.eval_repl(block, sandbox::Module, meta::Dict, doc::Documenter.Document, page)
    src_lines = Documenter.find_block_in_file(block.code, meta[:CurrentFile])
    # skip stats if there was a failure
    if length(doc.internal.errors) > 0
      invoke(Documenter.eval_repl, Tuple{Any,Any,Dict,Documenter.Document,Any}, block, sandbox, meta, doc, page)
    else
      println("page: $(Documenter.locrepr(meta[:CurrentFile], src_lines))")
      @time invoke(Documenter.eval_repl, Tuple{Any,Any,Dict,Documenter.Document,Any}, block, sandbox, meta, doc, page)
    end
  end
end
