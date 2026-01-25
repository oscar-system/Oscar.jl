# This script checks for missing jldoctest end markers and unknown admonitions.
#
# It is run as part of the OSCAR CI tests, but it can also be run manually
#
# It can be run from the command line with:
# > julia --project=. -e 'using Oscar; include("etc/check_docstrings.jl")'
#
# it is also possible to test more packages at once (as long as Main.X exists)
# > julia --project=. -e 'using Oscar; include("etc/check_docstrings.jl")' Nemo Singular Hecke Polymake GAP Oscar
#
# note: the @main function requires julia 1.11 or newer


using Markdown

admonition_types = string.((:note, :info, :tip, :danger, :warning, :compat, :todo, :details))

# this might not catch all broken jldoctests but seems to work for now
function has_broken_doctest(md::Markdown.MD)
  todo = copy(md.content)
  while !isempty(todo)
    elem = popfirst!(todo)
    # nested
    if elem isa Markdown.MD
      append!(todo, elem.content)
    else
      # proper docstrings are in a Markdown.Code block
      if elem isa Markdown.Paragraph
        for block in elem.content
          # unterminated jldoctests seem to end up inside some string in a Paragraph block
          if contains(string(block),"```jldoctest")
            return "Unterminated jldoctest: $(string(block))"
          end
        end
      elseif elem isa Markdown.Admonition
        if !(elem.category in admonition_types)
          return "Unknown admonition category: $(string(elem.category))"
        end
      end
    end
  end
  return nothing
end

function find_docstr_src(docs)
  locs = []
  res = docs.meta[:results]
  for i in 1:length(docs.content)
    msg = has_broken_doctest(docs.content[i])
    if msg !== nothing
      file = res[i].data[:path]
      line = res[i].data[:linenumber]
      mod = res[i].data[:module]
      push!(locs, (mod, file, line, msg))
    end
  end
  return locs
end

function get_broken_docstrings(m::Module)
  allpairs = collect(zip(Iterators.repeated(m), names(m; all=true, imported=false)));
  broken = []
  locs = []
  done = [m]

  while !isempty(allpairs)
    (src, name) = popfirst!(allpairs)
    memberobj = try
      getfield(src,name)
    catch
      nothing
    end
    if memberobj isa Module
      if !(memberobj in done)
        push!(done, memberobj)
        # we only want to consider one pkg
        # but note that we still seem to pick up some reexported names
        # even though imported=false and the check below
        if (pkgdir(memberobj) == pkgdir(m))
          newnames = names(memberobj; all=true, imported=false)
          filter!(x->!((memberobj, x) in allpairs || x in done), newnames)
          append!(allpairs, collect(zip(Iterators.repeated(memberobj), newnames)))
        end
      end
    elseif !isnothing(memberobj)
      docs = Base.Docs.doc(memberobj)
      if docs isa Markdown.MD
        if has_broken_doctest(docs) !== nothing && !(memberobj in broken)
          loc = find_docstr_src(docs)
          push!(broken, memberobj)
          # we could filter out docs from other packages via startswith(x, pkgdir(m))
          append!(locs, loc)
        end
      end
    end
  end
  return locs
end

function (@main)(args)
  modnames = length(args) > 0 ? args : ["Oscar"]
  mod = getfield.(Ref(Main), Symbol.(modnames))
  locs = reduce(vcat, get_broken_docstrings.(mod))
  isempty(locs) && exit(0)
  for (mod, file, line, msg) in locs
    relfile = relpath(file, pkgdir(mod))
    println("::error file=$relfile,line=$line,title=$(string(mod))::$msg")
  end
  exit(-1)
end
