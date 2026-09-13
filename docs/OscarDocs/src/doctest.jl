#
# Running doctests without building the manual.
#

# Documenter only offers to doctest a whole package. Testing a single function
# or file means feeding its docstrings to `Documenter._doctest` one by one, and
# that needs a `Documenter.Document` to record results in, which this creates.
function get_document(Oscar::Module; fix::Bool, set_meta::Bool)
  doc = Documenter.Document(root = joinpath(Base.pkgdir(Oscar), "docs"); doctest = fix ? :fix : true,
                            doctestfilters=Oscar.doctestfilters(), meta=Oscar.docmeta())

  if Documenter.DocMeta.getdocmeta(Oscar, :DocTestSetup) === nothing || set_meta
    Documenter.DocMeta.setdocmeta!(Oscar, :DocTestSetup, Oscar.doctestsetup(); recursive=true)
  end

  return doc
end

# The expected output in the manual and in the doctests was produced on an
# 80x24 terminal with OSCAR's unicode printing switched off. Run `f` under the
# same conditions, so that results do not depend on the caller's terminal.
function with_docs_terminal(f, Oscar::Module)
  withenv("COLUMNS" => 80, "LINES" => 24) do
    Oscar.with_unicode(false) do
      f()
    end
  end
end

@doc """
    doctest(Oscar::Module; fix::Bool = false)

Run all doctests of OSCAR, both those in docstrings and those in the manual
sources under `docs/src`.
"""
function doctest(Oscar::Module; fix::Bool = false)
  Documenter.DocMeta.setdocmeta!(Oscar, :DocTestSetup, Oscar.doctestsetup(); recursive=true)
  setup_doctest_stats(Oscar)
  with_docs_terminal(Oscar) do
    Documenter.doctest(Oscar; fix=fix, doctestfilters=Oscar.doctestfilters(),
                       meta=Oscar.docmeta())
  end
end

@doc """
    doctest(Oscar::Module, f::Function; fix::Bool = false, set_meta::Bool = false)

Run all doctests for the given function `f`.
"""
function doctest(Oscar::Module, f::Function; fix::Bool = false, set_meta::Bool = false)
  doc = get_document(Oscar; fix, set_meta)
  setup_doctest_stats(Oscar)

  with_docs_terminal(Oscar) do
    #essentially inspired by Documenter/src/DocTests.jl
    pm = parentmodule(f)
    bm = Base.Docs.meta(pm)
    md = bm[Base.Docs.Binding(pm, Symbol(f))]
    for s in md.order
      Documenter._doctest(md.docs[s], Oscar, doc)
    end
  end
end

@doc """
    doctest(Oscar::Module, path::String; fix::Bool = false, set_meta::Bool = false)

Run all doctests for all files in OSCAR where `path` occurs in the full
pathname.
"""
function doctest(Oscar::Module, path::String; fix::Bool = false, set_meta::Bool = false)
  doc = get_document(Oscar; fix, set_meta)
  setup_doctest_stats(Oscar)

  with_docs_terminal(Oscar) do
    walkmodules(Oscar) do m
      #essentially inspired by Documenter/src/DocTests.jl
      bm = Base.Docs.meta(m)
      for (_, md) in bm
        for s in md.order
          if occursin(path, md.docs[s].data[:path])
            Documenter._doctest(md.docs[s], Oscar, doc)
          end
        end
      end
    end
  end
end

# copied from JuliaTesting/Aqua.jl
function walkmodules(f, x::Module)
  f(x)
  for n in names(x; all=true)
    # `isdefined` and `getproperty` can trigger deprecation warnings
    if Base.isbindingresolved(x, n) && !Base.isdeprecated(x, n)
      isdefined(x, n) || continue
      y = getproperty(x, n)
      if y isa Module && y !== x && parentmodule(y) === x
        walkmodules(f, y)
      end
    end
  end
end
