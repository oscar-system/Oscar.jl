### Essential getters
maps(f::CompositeCoveredSchemeMorphism) = f.maps
map(f::CompositeCoveredSchemeMorphism, i::Int) = f.maps[i]
domain(f::CompositeCoveredSchemeMorphism) = domain(first(f.maps))
codomain(f::CompositeCoveredSchemeMorphism) = codomain(f.maps[end])

### Forwarding essential functionality (to be avoided!)
function underlying_morphism(f::CompositeCoveredSchemeMorphism)
  if !isdefined(f, :composed_map)
    result = underlying_morphism(first(maps(f)))::CoveredSchemeMorphism
    for i in 2:length(maps(f))
      result = compose(result, underlying_morphism(maps(f)[i]))::CoveredSchemeMorphism
    end
    f.composed_map = result
  end
  return f.composed_map::CoveredSchemeMorphism
end

### Specialized functionality

# Casting into the minimal concrete type for AbsCoveredSchemeMorphism
function CoveredSchemeMorphism(f::CompositeCoveredSchemeMorphism)
  return underlying_morphism(f)
end

########################################################################
# Printing   CompositeCoveredSchemeMorphism
########################################################################
function Base.show(io::IO, f::CompositeCoveredSchemeMorphism)
  io = pretty(io)
  if is_terse(io)
    print(io, "Composite morphism")
  else
    print(io, "Composition of ", "$(domain(f)) -> ")
    for i in 2:length(maps(f))
      print(io, "$(domain(maps(f)[i])) -> ")
    end
    print(io, "$(codomain(maps(f)[end]))")
  end
end

function Base.show(io::IO, ::MIME"text/plain", f::CompositeCoveredSchemeMorphism)
  io = pretty(io)
  println(io, "Composite morphism of", Indent())
  for g in maps(f)
    println(io, g)
  end
  println(io, Dedent())
end
