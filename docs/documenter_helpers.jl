# this function is slightly more specific than the one from documenter, print the corresponding code location
# calls the original function via invoke and prints a timing for the doctest
#
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
