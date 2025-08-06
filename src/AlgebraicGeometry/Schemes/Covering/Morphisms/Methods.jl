########################################################################
# Composition of CoveringMorphisms                                     #
########################################################################
function compose(f::CoveringMorphism, g::CoveringMorphism)
  domain(g) === codomain(f) || error("morphisms can not be composed")
  morphism_dict = IdDict{AbsAffineScheme, AbsAffineSchemeMor}()
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
  n = n_patches(C)
  new_patches = AbsAffineScheme[simplify(X) for X in patches(C)]
  GD = gluings(C)
  new_gluings = IdDict{Tuple{AbsAffineScheme, AbsAffineScheme}, AbsGluing}()
  for (X, Y) in keys(GD)
    Xsimp = new_patches[C[X]]
    iX, jX = identification_maps(Xsimp)
    Ysimp = new_patches[C[Y]]
    iY, jY = identification_maps(Ysimp)
    G = GD[(X, Y)]
    #new_gluings[(Xsimp, Ysimp)] = restrict(G, jX, jY, check=false)
    new_gluings[(Xsimp, Ysimp)] = LazyGluing(Xsimp, Ysimp,
                                               RestrictionDataIsomorphism(G, jX, jY)
                                              )
  end
  iDict = IdDict{AbsAffineScheme, AbsAffineSchemeMor}()
  jDict = IdDict{AbsAffineScheme, AbsAffineSchemeMor}()
  for i in 1:length(new_patches)
    iDict[new_patches[i]] = identification_maps(new_patches[i])[1]
    jDict[C[i]] = identification_maps(new_patches[i])[2]
  end
  Cnew = Covering(new_patches, new_gluings, check=false)
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
function base_change(phi::Any, f::CoveringMorphism{<:Any, <:Any, MorphismType, BaseMorType};
    domain_map::CoveringMorphism=base_change(phi, domain(f))[2],
    codomain_map::CoveringMorphism=base_change(phi, codomain(f))[2]
  ) where {MorphismType, BaseMorType}
  D = domain(f)
  C = codomain(f)
  DD = domain(domain_map)
  CC = domain(codomain_map)
  mor_dict = IdDict{AbsAffineScheme, MorphismType}()
  for UU in patches(DD)
    U = codomain(domain_map[UU])
    V = codomain(f[U])
    g_V = first(maps_with_given_codomain(codomain_map, V)) # The result must be unique as it arises 
                                                           # from a base change.
    _, ff, _ = base_change(phi, f[U], domain_map=domain_map[UU], codomain_map=g_V)
    mor_dict[UU] = ff
  end

  return domain_map, CoveringMorphism(DD, CC, mor_dict, check=false), codomain_map
end

###############################################################################
#
#  Printing
#
###############################################################################

function Base.show(io::IO, f::CoveringMorphism)
  io = pretty(io)
  if get(io, :show_semi_compact, false)
    _show_semi_compact(io, f)
  elseif is_terse(io)
    print(io, "Covering morphism")
  else
    print(io, "Hom: ", Lowercase(), domain(f), " -> ", Lowercase(), codomain(f))
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
  println(io, "Covering morphism")
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
  print(io, "to ", Lowercase(), codomain(f))
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

########################################################################
# fiber products of Coverings
#
# For any two `CoveringMorphism`s `f : A → C` and `g : B → C` we 
# define the `fiber_product` to be the covering `A×B` which makes 
# the obvious square commute:
#        f'
#    A×B → B
#   g'↓    ↓ g
#     A →  C
#       f
#
# i.e. the `patches` of `A×B` are Uᵢ× Vⱼ for the `patches` Uᵢ
# of `A` and `Vⱼ` of `B`.
#
# In special cases, e.g. when `f` is a refinement, the creation of 
# the patches Uᵢ× Vⱼ and its associated maps can be significantly 
# simplified. This is taken care of by the dispatch for fiber products
# of morphisms of affine schemes and it is assumed that the morphisms 
# of a refinement is of type `PrincipalOpenEmbedding`.
########################################################################

function fiber_product(f::CoveringMorphism, g::CoveringMorphism)
  A = domain(f)
  B = domain(g)
  C = codomain(f)
  @assert C === codomain(g) "codomains do not agree"

  new_patches = AbsAffineScheme[]
  cache = IdDict{AbsAffineScheme, Tuple{<:AbsAffineSchemeMor, <:AbsAffineSchemeMor}}()
  maps_to_A = IdDict{AbsAffineScheme, AbsAffineSchemeMor}()
  maps_to_B = IdDict{AbsAffineScheme, AbsAffineSchemeMor}()
  for U in patches(A)
    f_U = f[U]
    W = codomain(f_U)
    for V in patches(B) 
      g_V = g[V]
      codomain(g_V) === W || continue
      UxV, to_U, to_V = fiber_product(f_U, g_V)
      push!(new_patches, UxV)
      cache[UxV] = (to_U, to_V)
      maps_to_A[UxV] = to_U
      maps_to_B[UxV] = to_V
    end
  end
  result = Covering(new_patches)

  # construct all the gluings
  # TODO: Make these lazy!
  for UxV in new_patches
    to_U, to_V = cache[UxV]
    U = codomain(to_U)
    V = codomain(to_V)
    for UUxVV in new_patches
      to_UU, to_VV = cache[UUxVV]
      UU = codomain(to_UU)
      VV = codomain(to_VV)
      !haskey(gluings(A), (U, UU)) && continue
      !haskey(gluings(B), (V, VV)) && continue
      UUU = A[U, UU]::AbsGluing
      VVV = B[V, VV]::AbsGluing
      U_UU, UU_U = gluing_domains(UUU)
      V_VV, VV_V = gluing_domains(VVV)
      inc_U = inclusion_morphism(U_UU)
      inc_UU = inclusion_morphism(UU_U)
      inc_V = inclusion_morphism(V_VV)
      inc_VV = inclusion_morphism(VV_V)

      
      # assemble the new gluing domains:
      h_U = complement_equation(U_UU)
      h_V = complement_equation(V_VV)
      h_UV = [pullback(to_U)(h_U), pullback(to_V)(h_V)]
      U_UUxV_VV = PrincipalOpenSubset(UxV, h_UV)
      f_res_U = restrict(to_U, U_UUxV_VV, U_UU, check=false)
      g_res_V = restrict(to_V, U_UUxV_VV, V_VV, check=false)
      # TODO: Can we produce U_UUxV_VV as a PrincipalOpenSubset of UxV 
      # with a generic `fiber_product` routine? If so, what should the 
      # signature and functionality be? We need it as a PrincipalOpenSubset
      # because of the convention for the `SimpleGluing`s.
      #   U_UUxV_VV, f_res_U, g_res_V = fiber_product(restrict(f, U_UU, codomain(f)), restrict(g, V_VV, codomain(g)))

      h_UU = complement_equation(UU_U)
      h_VV = complement_equation(VV_V)
      h_UUVV = [pullback(to_UU)(h_UU), pullback(to_VV)(h_VV)]
      UU_UxVV_V = PrincipalOpenSubset(UUxVV, h_UUVV)
      f_res_UU = restrict(to_UU, UU_UxVV_V, UU_U, check=false)
      g_res_VV = restrict(to_VV, UU_UxVV_V, VV_V, check=false)

      simple_to_double_U, double_to_simple_U = gluing_morphisms(UUU)
      simple_to_double_V, double_to_simple_V = gluing_morphisms(VVV)

      # construct the gluing morphisms
      # Since the fiber products have not been constructed in the appropriate way, 
      # we can not use the generic method.
      #=
      simple_to_double = induced_map_to_fiber_product(compose(f_res_U, simple_to_double_U), compose(g_res_V, simple_to_double_V), 
                                                      restrict(f[UU], UU_U, codomain(f[UU]), check=true), 
                                                      restrict(g[VV], VV_V, codomain(g[VV]), check=true),
                                                      fiber_product=(UU_UxVV_V, f_res_UU, g_res_VV)
                                                     )
      double_to_simple = induced_map_to_fiber_product(compose(f_res_UU, double_to_simple_U), compose(g_res_VV, double_to_simple_V), 
                                                      restrict(f[U], U_UU, codomain(f[U]), check=true), 
                                                      restrict(g[V], V_VV, codomain(g[V]), check=true),
                                                      fiber_product=(U_UUxV_VV, f_res_U, g_res_V)
                                                     )
      =#
      pre_glue_double_to_simple = induced_map_to_fiber_product(
           compose(f_res_UU, compose(double_to_simple_U, inc_U)),
           compose(g_res_VV, compose(double_to_simple_V, inc_V)),
           f[U], g[V];
           fiber_product=(UxV, to_U, to_V),
           check=false
          )
      double_to_simple = restrict(pre_glue_double_to_simple, domain(pre_glue_double_to_simple), 
                                  U_UUxV_VV; check=false
                                 )
      pre_glue_simple_to_double = induced_map_to_fiber_product(
           compose(f_res_U, compose(simple_to_double_U, inc_UU)),
           compose(g_res_V, compose(simple_to_double_V, inc_VV)),
           f[UU], g[VV];
           fiber_product=(UUxVV, to_UU, to_VV),
           check=false
          )
      simple_to_double = restrict(pre_glue_simple_to_double, domain(pre_glue_simple_to_double), 
                                  UU_UxVV_V; check=false
                                 )

      new_gluing = Gluing(UxV, UUxVV, simple_to_double, double_to_simple; check=false)
      add_gluing!(result, new_gluing)
    end
  end
  to_A = CoveringMorphism(result, A, maps_to_A; check=false)
  to_B = CoveringMorphism(result, B, maps_to_B; check=false)
  return result, to_A, to_B
end

