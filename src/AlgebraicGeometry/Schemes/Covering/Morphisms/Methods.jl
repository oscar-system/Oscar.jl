########################################################################
# Composition of CoveringMorphisms                                     #
########################################################################
function compose(f::CoveringMorphism, g::CoveringMorphism)
  domain(g) === codomain(f) || error("morphisms can not be composed")
  morphism_dict = IdDict{AbsSpec, AbsSpecMor}()
  for U in patches(domain(f))
    morphism_dict[U] = compose(f[U], g[codomain(f[U])])
  end
  return CoveringMorphism(domain(f), codomain(g), morphism_dict, check=false)
end

########################################################################
# Simplification of Coverings                                          #
########################################################################
@doc raw"""
    simplify(C::Covering)

Given a covering ``C`` apply `simplify` to all basic affine patches 
in ``C`` and return a triple ``(C', f, g)`` consisting of the 
resulting covering ``C'`` and the identifying isomorphism 
``f : C' ↔ C``.
"""
function simplify(C::Covering)
  n = npatches(C)
  new_patches = AbsSpec[simplify(X) for X in patches(C)]
  GD = glueings(C)
  new_glueings = IdDict{Tuple{AbsSpec, AbsSpec}, AbsGlueing}()
  for (X, Y) in keys(GD)
    Xsimp = new_patches[C[X]]
    iX, jX = identification_maps(Xsimp)
    Ysimp = new_patches[C[Y]]
    iY, jY = identification_maps(Ysimp)
    G = GD[(X, Y)]
    #new_glueings[(Xsimp, Ysimp)] = restrict(G, jX, jY, check=false)
    new_glueings[(Xsimp, Ysimp)] = LazyGlueing(Xsimp, Ysimp, _compute_restriction, _compute_domains,
                                               RestrictionDataIsomorphism(G, jX, jY)
                                              )
  end
  iDict = IdDict{AbsSpec, AbsSpecMor}()
  jDict = IdDict{AbsSpec, AbsSpecMor}()
  for i in 1:length(new_patches)
    iDict[new_patches[i]] = identification_maps(new_patches[i])[1]
    jDict[C[i]] = identification_maps(new_patches[i])[2]
  end
  Cnew = Covering(new_patches, new_glueings, check=false)
  i_cov_mor = CoveringMorphism(Cnew, C, iDict, check=false)
  j_cov_mor = CoveringMorphism(C, Cnew, jDict, check=false)
  
  # Carry over the decomposition information.
  if has_decomposition_info(C)
    for U in new_patches
      V = codomain(i_cov_mor[U])
      pb = pullback(i_cov_mor[U])
      set_decomposition_info!(Cnew, U, elem_type(OO(U))[pb(a) for a in decomposition_info(C)[V]])
    end
  end

  return Cnew, i_cov_mor, j_cov_mor
end

########################################################################
# Base change
########################################################################
function base_change(phi::Any, f::CoveringMorphism;
    domain_map::CoveringMorphism=base_change(phi, domain(f))[2],
    codomain_map::CoveringMorphism=base_change(phi, codomain(f))[2]
  )
  D = domain(f)
  C = codomain(f)
  DD = domain(domain_map)
  CC = domain(codomain_map)
  mor_dict = IdDict{AbsSpec, AbsSpecMor}()
  for UU in patches(DD)
    U = codomain(domain_map[UU])
    V = codomain(f[U])
    g_V = first(maps_with_given_codomain(codomain_map, V)) # The result must be unique as it arises 
                                                           # from a base change.
    _, ff, _ = base_change(phi, f[U], domain_map=domain_map[UU], codomain_map=g_V)
    mor_dict[UU] = ff
  end

  return domain_map, CoveringMorphism(DD, CC, mor_dict, check=true), codomain_map # TODO: Set to false after testing.
end

###############################################################################
#
#  Printing
#
###############################################################################

function Base.show(io::IO, f::CoveringMorphism)
  io = pretty(io)
  if get(io, :supercompact, false)
    print(io, "Morphism")
  else
    print(io, "Morphism: ", Lowercase(), domain(f), " -> ", Lowercase(), codomain(f))
  end
end

# Like the detailed printing but we avoid to mention the domain and codomain -
# this printing is used in nested printed where we assumed that the domain and
# codomain of `f` are already in the nest.
function _show_semi_compact(io::IO, f::CoveringMorphism)
  io = pretty(io)
  mor = morphisms(f)
  lX = ndigits(length(domain(f)))
  lY = ndigits(length(codomain(f)))
  for i in 1:length(domain(f))
    li = ndigits(i)
    U = domain(f)[i]
    g = mor[U]
    j = findfirst(V -> codomain(g) === V, collect(codomain(f)))
    lj = ndigits(j)
    print(io, " "^(lX-li)*"$(i)a -> "*" "^(lY-lj)*"$(j)b")
    print(io, Indent())
    x = coordinates(codomain(g))
    pg = pullback(mor[U])
    co_str = String["$(y)" for y in x]
    pushfirst!(co_str, "")
    k = max(length.(co_str)...)
    for j in 1:length(x)
      kj = length(co_str[j+1])
      println(io)
      print(io, " "^(k-kj)*"$(x[j]) -> $(pg(x[j]))")
    end
    if i != length(patches(domain(f)))
      println(io)
      println(io, "----------------------------------------")
    end
    print(io, Dedent())
  end
end

# As for covered scheme morphisms, we mention the domain and the codomain with a
# fix covering. The charts in each covering are labeled by a number and a letter
# - "a" for the domain and "b" for the codomain.
#
# We need to take care of some offsets for the alignments:
# - when mentioning the coordinates of all charts of a given covering, we take
# care of the length of the longest printed system of coordinated - this gives
# us a sample to use to create our offsets. Doing so, we align our coordinates
# systems on the left, and the description of the charts are also aligned on the
# left, 3 spaces after the systems of coordinates;
# - for the labels, we have to take care that if we have more that 10 charts,
# the small labels will need a left offset in order to align all labels on the
# right
# - for the description on the map by pullbacks between given charts, we want to
# have all the arrows at the same level - to do so, we take care of some left
# offsets so that the elements in the "domain of the pullback morphism" are
# aligned on the right, and their respective images are aligned on the left.
function Base.show(io::IO, ::MIME"text/plain", f::CoveringMorphism)
  io = pretty(io)
  println(io, "Morphism")
  print(io, Indent(), "from ", Lowercase(), domain(f))
  print(io, Indent())
  co_str = String[""]
  lX = ndigits(length(domain(f)))
  lY = ndigits(length(codomain(f)))
  for i in 1:length(domain(f))
    U = domain(f)[i]
    co = coordinates(U)
    str = "["*join(co, ", ")*"]"
    push!(co_str, str)
  end
  k = max(length.(co_str)...)
  for i in 1:length(domain(f))
    li = ndigits(i)
    U = domain(f)[i]
    kc = length(co_str[i+1])
    println(io)
    print(io, " "^(lX-li)*"$(i)a: "*co_str[i+1]*" "^(k-kc+3), Lowercase(), U)
  end
  println(io, Dedent())
  print(io, "to   ", Lowercase(), codomain(f))
  print(io, Indent())
  co_str = String[""]
  for i in 1:length(codomain(f))
    U = codomain(f)[i]
    co = coordinates(U)
    str = "["*join(co, ", ")*"]"
    push!(co_str, str)
  end
  k = max(length.(co_str)...)
  for i in 1:length(codomain(f))
    li = ndigits(i)
    U = codomain(f)[i]
    kc = length(co_str[i+1]) 
    println(io)
    print(io, " "^(lY-li)*"$(i)b: "*co_str[i+1]*" "^(k-kc+3), Lowercase(), U)
  end
  print(io, Dedent(), Dedent())
  mor = morphisms(f)
  if length(mor) > 0
    println(io)
    print(io, "given by the pullback function")
    length(mor) > 1 && print(io, "s")
    println(io, Indent())
    for i in 1:length(domain(f))
      li = ndigits(i)
      U = domain(f)[i]
      g = mor[U]
      j = findfirst(V -> codomain(g) === V, collect(codomain(f)))
      lj = ndigits(j)
      print(io, " "^(lX-li)*"$(i)a -> "*" "^(lY-lj)*"$(j)b")
      print(io, Indent())
      x = coordinates(codomain(g))
      co_str = String["$(y)" for y in x]
      pushfirst!(co_str, "")
      k = max(length.(co_str)...)
      pg = pullback(mor[U])
      for j in 1:length(x)
        kj = length(co_str[j+1])
        println(io)
        print(io, " "^(k-kj)*"$(x[j]) -> $(pg(x[j]))")
      end
      print(io, Dedent())
      if i != length(patches(domain(f)))
        println(io)
        println(io, "----------------------------------------")
      end
    end
    print(io, Dedent(), Dedent())
  end
end


