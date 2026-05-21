@testset "finite extensions" begin
  # Create a `FiniteExtension`
  S, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z])
  P, (u, v) = graded_polynomial_ring(QQ, [:u, :v])
  f = x^2 + y^2 + z^2
  I = ideal(S, f)
  A, _ = quo(S, I)
  psi = Oscar.FiniteExtension(hom(P, A, A.([x, y])));

  # As a test case we recreate the codomain as a free module over itself
  # and push it forward to the domain.
  F = graded_free_module(A, [zero(grading_group(A))])
  pf_F, interp = pushforward(psi, F)
  @test base_ring(pf_F) === P

  # Now we create the bidual F** as a `P`-module.
  dd_pf_F, dd = double_dual(pf_F)

  # The bidual has a natural structure as an `A`-module which 
  # we will now recover. To this end, we first take the tensor 
  # product - ⊗_P A again.
  M, to_M = change_base_ring(psi, dd_pf_F)

  # We need additional relations on the raw tensor product, 
  # namely for each `x` in the `generating_system` of `A` over `P`
  # we should have 
  #
  #   g ⊗ x == x⋅g ⊗ 1
  #
  # for any element `g` in F**. Note that the multiplication 
  # on the right hand side is to be carried out as: 
  #
  #   F** ∋ λ : (φ : F → P) ↦ λ(φ) 
  #
  # goes to 
  # 
  #   x⋅λ : (φ : F → P) ↦ λ(φ∘(v↦ x⋅v))
  #
  # Here `F` is the original module seen as both, an `A` and a `P`-module.

  # We first construct the endomorphisms of `pf_F` given by multiplication 
  # with the elements of the `generating_system` of `A` over `P`.
  mult_maps = []
  for x in generating_system(psi)
    img_gens = gens(pf_F)
    img_gens = [preimage(interp, v) for v in img_gens]
    img_gens = [x*v for v in img_gens]
    img_gens = [interp(v) for v in img_gens]
    f = hom(pf_F, pf_F, img_gens)
    push!(mult_maps, f)
  end

  d_pf_F, cod = get_attribute(dd_pf_F, :hom)

  # For every generator `Phi` of F** and every `x` in the 
  # `generating_system` of `A` over `P` we build the element 
  # for `x ⋅Phi`.
  ms = []
  for Phi in gens(dd_pf_F)
    mms = []
    for f in mult_maps
      img_gens = elem_type(cod)[]
      for e in gens(d_pf_F)
        g = element_to_homomorphism(e)
        h = compose(f, g)
        v = homomorphism_to_element(d_pf_F, h)
        push!(img_gens, element_to_homomorphism(Phi)(v))
      end
      m = hom(d_pf_F, cod, img_gens)
      push!(mms, homomorphism_to_element(dd_pf_F, m))
    end
    push!(ms, mms)
  end

  # From this we can assemble the relations which need to be 
  # added to the raw tensor product F** ⊗_P A in order to 
  # recover the `A`-module structure on F**
  new_rels = []
  for (g, mms) in zip(gens(dd_pf_F), ms)
    gg = to_M(g)
    l = []
    for (x, f) in zip(generating_system(psi), mms)
      xx = A(x)
      ggg = xx*gg
      push!(l, ggg - to_M(f))
    end
    push!(new_rels, l)
  end

  nr = reduce(vcat, new_rels)

  N, _ = quo(M, elem_type(M)[v for v in nr])

  # We test whether this module is actually free
  NN, _ = simplify(N)
  @test ngens(NN) == 1
  @test length(relations(NN)) == 0
end

