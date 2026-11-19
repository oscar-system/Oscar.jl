# Credits for software and papers that Oscar uses. See `Credit`, `@credit` and
# `@collect_credits` in AbstractAlgebra. Define a credit by its key in
# docs/oscar_references.bib, as a `const` next to the code that uses it:
#
#   const CREDIT_Kem99 = _credit("Kem99")
#   ...
#   @credit CREDIT_Kem99

const _bibtex_entries = let file = joinpath(@__DIR__, "..", "..", "docs", "oscar_references.bib")
  include_dependency(file)
  text = read(file, String)
  entries = Dict{String, String}()
  for m in eachmatch(r"^@(\w+)\s*\{\s*([^,\s]+)\s*,"m, text)
    lowercase(m.captures[1]) in ("comment", "string", "preamble") && continue
    # the entry ends at the brace that closes the first `{`
    i, depth = m.offset, 0
    while true
      c = codeunit(text, i)
      depth += (c == UInt8('{')) - (c == UInt8('}'))
      c == UInt8('}') && depth == 0 && break
      i += 1
    end
    entries[m.captures[2]] = text[m.offset:i]
  end
  entries
end

function _credit(key::String)
  haskey(_bibtex_entries, key) || error("key $key is not in docs/oscar_references.bib")
  return Credit(key, _bibtex_entries[key])
end
