################################################################################
#
#  Helper
#
################################################################################

op(A) = *

op(A::FinGenAbGroup) = +

################################################################################

function effective_presentation(A::FinGenAbGroup)
  forward = isomorphism(FPGroup, A)
  G = codomain(forward)
  set_order(G, order(A))
  S, GtoS = simplified_fp_group(G)
  return EffectivePresentation(A, S, forward * GtoS, inv(GtoS) * inv(forward))
end

function effective_presentation(A::GAPGroup)
  forward = isomorphism(FPGroup, A)
  G = codomain(forward)
  if Oscar.has_order(A)
    set_order(G, order(A))
  end
  return EffectivePresentation(A, G, forward, inv(forward))
end

# Assume that A is abstract presentation of a group A
# and f : A -> H, g : H -> A are isomorphisms
function EffectivePresentation(A::EffectivePresentation, H, HtoA, AtoH)
  return EffectivePresentation(H, A.G, x -> 
                                begin
                                  @assert parent(x) === H
                                 A.forward(HtoA(x)) 
                                end
                                , y -> AtoH(A.backward(y)))
end

# treat as map

function (A::EffectivePresentation)(x)
  return A.forward(x)
end

function preimage(A::EffectivePresentation, y)
  return A.backward(y)
end

################################################################################
#
#  Extensions of effective presentations
#
################################################################################

function extension(A::EffectivePresentation, C::EffectivePresentation, B, f, fpreim, g, gpreim)
                          
  #                       f    g
  # Assume we have 1 -> A -> B -> C -> 1
  #                     |         |
  #                     v         v
  #                1 -> N -> G -> H -> 1
  #                       i    pi
  #
  # A and C are finitely presented groups,
  # N, G, H are black box groups, i, pi morphisms and A -> N, C -> H isomorphisms
  # I know all the maps. I compute a finitely presented group B and
  # the canonical map B -> G making everything commute
  
  # This is Handbook, 2.4.3
  N = A.G
  H = C.G
  @vprintln :FiniteRings "Computing extension of "
  @vprintln :FiniteRings "  N: $(_show_fp(N))"
  @vprintln :FiniteRings "  H: $(_show_fp(H))"
  rH = ngens(H)
  rN = ngens(N)
  F = free_group(rH + rN, "x")
  gNinG = gens(F)[rH + 1:end]
  gHinG = gens(F)[1:rH]
  new_relators = elem_type(F)[]
  gensN_in_B = [f(A.backward(c)) for c in gens(N)]
  gensH_in_B = [gpreim(C.backward(c)) for c in gens(H)]
  gensH_in_B_inv = inv.(gensH_in_B)
  all_gens = vcat(gensH_in_B, gensN_in_B)
  all_gens_inv = inv.(all_gens)
  @vprintln :FiniteRings "Creating new relators ..."
  for r in relators(H)
    rpreim = map_word(r, gHinG) # this is the lifted word
    reval = map_word(r, gensH_in_B; genimgs_inv = gensH_in_B_inv) # this is the lifted word, evaluated in B
    # the lifted (evaluated) word is in ker(g) = im(f)
    kernel_word = A.forward(fpreim(reval)) # this is the word in N
    # now moveA to G
    kernel_word_in_B = map_word(kernel_word, gNinG; init = one(F))
    # take the difference which is now a relator for the generators in G
    push!(new_relators, rpreim * kernel_word_in_B^-1)
  end
  # the relators of N are mapped through
  for r in relators(N)
    push!(new_relators, map_word(r, gNinG))
  end
  for i in 1:length(gNinG)
    y = gNinG[i]
    yinA = f(A.backward(gen(N, i)))
    for j in 1:length(gHinG)
      x = gHinG[j]
      xinA = gensH_in_B[j]
      wxy = map_word(A.forward(fpreim(inv(xinA) * yinA * xinA)), gNinG)
      push!(new_relators, inv(x) * y * x * inv(wxy))
    end
  end
  @vprint :FiniteRings "Constructing new fp group "
  if ngens(F) == 0
    G, FtoG = F, id_hom(F)
  else
    G, FtoG =  quo(F, new_relators)
  end
  if Oscar.has_order(N) && Oscar.has_order(H)
    Oscar.set_order(G, order(N) * order(H))
  end
  @vprintln :FiniteRings "G: $(_show_fp(G))"
  GtoH = hom(G, H, vcat(gens(H), [one(H) for i in 1:rN], ); check = false)
  @vprintln :FiniteRings "Constructing homomorphisms of fp groups"
  NtoG = hom(N, G, FtoG.(gNinG); check = false)
  forward = function(b)
    @assert parent(b) === B
    binH = C.forward(g(b))
    binHlifted = map_word(binH, gHinG; init = one(F)) # this is a word X
     y = inv(b) * map_word(binH, gensH_in_B; init = one(B), genimgs_inv = gensH_in_B_inv)
     a = A.forward(fpreim(y))
     ia = map_word(a, gNinG; init = one(F))
     res = image(FtoG, binHlifted * inv(ia))
     @assert parent(res) === G
     return res
  end
  backward = function(g)
    uuuu = map_word(g, all_gens, genimgs_inv = all_gens_inv; init = one(B))
    return uuuu
  end
  if true ngens((G)) >= 40
    @vprintln :FiniteRings "Not simplifying"
    GG, GtoGG = G, id_hom(G)
  else
    @vprint :FiniteRings "Simplifying ... "
    GG, GtoGG = simplified_fp_group(G)
    @vprintln :FiniteRings "to: $(_show_fp(GG))"
  end
  return EffectivePresentation(B, GG, x -> 
                              begin
                                xx = forward(x)
                                xxx = GtoGG(xx)
                                return xxx
                              end, z -> backward(GtoGG\(z)))

end

function extension(A::EffectivePresentation, C::EffectivePresentation, B, f, g)
  fpreim = x -> begin
    fl, y = has_preimage_with_preimage(f, x; check = false)
    @assert fl
    return y
  end

  gpreim = x -> begin
    fl, y = has_preimage_with_preimage(g, x; check = false)
    @assert fl
    return y
  end

  return extension(A, C, B, f, fpreim, g, gpreim)
end

################################################################################
#
#  Printing
#
################################################################################

function Base.show(io::IO, A::EffectivePresentation)
  io = pretty(io)
  println(io, "Effective presentation")
  println(io, Indent(), "of ", Lowercase(), A.A)
  print(io, "as ", Lowercase(), A.G)
  print(io, Dedent())
end

function _show_fp(G)
  return "order: $(Oscar.has_order(G) ? order(G) : "??"), $(ngens(G)) generators, $(length(relators(G))) relators"
end

