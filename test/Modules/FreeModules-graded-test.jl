function _rand_polys_nonzero(R, n)
  polys = elem_type(R)[]
  while length(polys) < n
    f = MPolyBuildCtx(R)
    for z = 1:3
      e = rand(0:2, ngens(R))
      r = base_ring(R)(rand(0:22))
      push_term!(f, r, e)
    end
    g = finish(f)
    if !iszero(g)
      push!(polys, g)
    end
  end
  return polys
end

function _eq(A::Oscar.SubQuo_dec, B::Oscar.SubQuo_dec)
  if A.F != B.F
    return false
  end
  Oscar.singular_assure(A.sum)
  Oscar.singular_assure(B.sum)
  if !(iszero(Singular.lift(A.sum.S, B.sum.S)[2])) || !(iszero(Singular.lift(B.sum.S, A.sum.S)[2]))
    return false
  end
  if isdefined(A, :quo) || isdefined(B, :quo)
    if !isdefined(A, :quo) || !isdefined(B, :quo)
      return false
    end
    Oscar.singular_assure(A.quo)
    Oscar.singular_assure(B.quo)
    if !(iszero(Singular.lift(A.quo.S, B.quo.S)[2])) || !(iszero(Singular.lift(B.quo.S, A.quo.S)[2]))
      return false
    end
  end
  return true
end

function _sparse_to_array(t::Oscar.FreeModuleElem_dec, F::Oscar.FreeModule_dec)
  res = [zero(F) for i = 1:3]
  for (i,g) = t.r
    res[div(i-1, 3) + 1] += g * gen(F, (i-1)%3 + 1)
  end
  return res
end

@testset "Modules" begin
  Qx, (x,y,z) = PolynomialRing(QQ, ["x", "y", "z"])
  t = gen(Hecke.Globals.Qx)
  k1, l = number_field(t + 3)
  NFx = PolynomialRing(k1, ["x", "y", "z"])[1]
  k2 = Nemo.GF(23)
  GFx = PolynomialRing(k2, ["x", "y", "z"])[1]
  RNmodx = PolynomialRing(Nemo.ResidueRing(ZZ,17), :x => 1:3)[1]
  Rings = [Qx, NFx, GFx, RNmodx]

  A = abelian_group([0 3 0; 2 1 2])
  GrpElems = elem_type(A)[A(fmpz[0, 1, 1]), A(fmpz[0, 1, 0]), A(fmpz[1, 2, 0])]

  Rings_dec = []
  v = [1, 2, 4]

  for R in Rings
    
    decorated_rings = [decorate(R)[1],
                     grade(R, [v[i] for i=1:ngens(R)])[1],
                     filtrate(R, [v[i] for i=1:ngens(R)])[1],
                     grade(R, [GrpElems[i] for i = 1:ngens(R)])[1]]


    for (j, RR) in enumerate(decorated_rings)
      G = RR.D
      if j == 4
        Elems = [G([convert(fmpz,x) for x = [0,4,1]]), G([convert(fmpz,x) for x = [2,1,0]]), G([convert(fmpz,x) for x = [1,0,1]])]
      else
        Elems = [G([convert(fmpz,1)]), G([convert(fmpz,2)]), G([convert(fmpz,3)])]
      end
      
      F1, F2 = free_module(RR,3), free_module(RR, Elems)
      Mods = [F1, F2]

      polys = _rand_polys_nonzero(R, 18)
      @assert all(!iszero, polys)

      for F in Mods
        
        @test Oscar.is_graded(F) == Oscar.is_graded(RR)
        @test Oscar.is_filtered(F) == Oscar.is_filtered(RR)

        G = grading_group(F)
        if j == 4
          a = G([convert(fmpz,x) for x = [1,0,1]])
        else
          a = G([convert(fmpz,5)])
        end
        b = (F)(a)
        #@test parent(b) === F
        @test (F)() == zero(F)
        #@test ngens(R^3) == 3
        FreeModElems = [polys[c*3+1]*gen(F,1) + polys[c*3+2]*gen(F,2) + polys[c*3+3]*gen(F,3) for c = 0:5]
        @test parent_type(FreeModElems[1]) == typeof(F)
        #@test (5::Integer)*((4::Int)* (-(FreeModElems[1]))) == QQ(-20) * FreeModElems[1]
        Oscar.BiModArray(FreeModElems, F)
        Hom_FreeModElemst = []
        len = Dict{GrpAbFinGenElem, Int64}()
        for c = 1:6
          for k in keys(homogeneous_components(FreeModElems[c]))
            if haskey(len, k)
              len[k] += 1
            else
              len[k] = 1
            end
          end
        end
        for p = 1:6
          max = 0
          temp = undef
          for k in keys(len)
            if len[k] >= max
              max = len[k]
              temp = k
            end
          end
          len[temp] = 0
          comp = zero(F)
          for l = 1:6
            comp += homogeneous_component(FreeModElems[l], temp)
          end
          push!(Hom_FreeModElemst, comp)
        end
        order_old = [1,2,3,4,5,6]
        order_new = []
        for p = 1:6
          t = rand(1:(7-p))
          push!(order_new, order_old[t])
          deleteat!(order_old, t)
        end
        Hom_FreeModElems = []
        for p = 1:6
          push!(Hom_FreeModElems, Hom_FreeModElemst[order_new[p]])
        end
        hom_keys = [collect(keys(homogeneous_components(F.R(polys[s]))))[1] for s = 1:6]
        hom_pols = [homogeneous_component(F.R(polys[s]), hom_keys[s]) for s = 1:6]
        SubQuos = [sub(F, [Hom_FreeModElems[e] for e = 1:3]), quo(F, [Hom_FreeModElems[2*e-1] for e = 1:3])]
        Hom_SubQuoElems = [[SubQuos[1](SubQuos[1](hom_pols[e] * Hom_FreeModElems[e])) for e = 1:3], [SubQuos[2](Hom_FreeModElems[2*e]) for e = 1:3]]
        @test parent_type(Hom_SubQuoElems[1][1]) == typeof(SubQuos[1])
        #@test (5::Integer)*((4::Int)* (-(Hom_SubQuoElems[1][1]))) == QQ(-20) * Hom_SubQuoElems[1][1]
        if !iszero(Hom_SubQuoElems[1][1])
          @test degree(Hom_SubQuoElems[1][1].a, SubQuos[1]) == degree(Hom_SubQuoElems[1][1])
        end
        non_zero = true
        temp = F
        for b = 1:3
          if iszero(Hom_SubQuoElems[2][b])
            non_zero = false
          end
        end
        if !non_zero
          SubQuos[2] = quo(F, [gen(F.R, 2) * gen(F, e) for e = 1:3])
          Hom_SubQuoElems[2] = [SubQuos[2](gen(F.R, 1) * gen(F, e)) for e = 1:3]
          for e = 1:3
            Hom_FreeModElems[2*e - 1] = gen(F.R, 2) * gen(F, e)
          end
          SubQuos[1] = sub(F, [Hom_FreeModElems[e] for e = 1:3])
        end
        push!(SubQuos, quo(SubQuos[1], [Hom_FreeModElems[e] for e = 4:6]), quo(SubQuos[1], Hom_SubQuoElems[2]))
        push!(Hom_SubQuoElems, [SubQuos[3](hom_pols[e+3] * Hom_FreeModElems[e]) for e = 1:3])
        non_zero = true
        for b = 1:3
          if iszero(Hom_SubQuoElems[3][b]) || iszero(gen(SubQuos[3], b))
            non_zero = false
          end
        end
        if !non_zero
          SubQuos[1] = sub(F, [gen(F.R, 1) * gen(F, e) for e = 1:3])
          SubQuos[4] = quo(SubQuos[1], Hom_SubQuoElems[2])
          SubQuos[3] = quo(SubQuos[1], [gen(F.R, 2) * gen(SubQuos[1], e) for e = 1:3])
          Hom_SubQuoElems[3] = gens(SubQuos[3])
          Hom_SubQuoElems[1] = gens(SubQuos[1])
        end
					        
        @test _eq(sub(F, Hom_SubQuoElems[1]), sub(SubQuos[1], Hom_SubQuoElems[1]))
        @test _eq(sub(F, SubQuos[1]), SubQuos[1])
        @test _eq(quo(F, Hom_SubQuoElems[2]), quo(SubQuos[2], Hom_SubQuoElems[2]))
        @test _eq(quo(F, [Hom_FreeModElems[e] for e = 1:6]), quo(SubQuos[2], [Hom_FreeModElems[2*e] for e = 1:3]))
        @test _eq(quo(SubQuos[1], SubQuos[2]), quo(F,gens(F)))
        @test _eq(quo(sub(F, gens(F, F)), sub(F, [Hom_FreeModElems[2*e - 1] for e = 1:3])), quo(F, [Hom_FreeModElems[2*e - 1] for e = 1:3]))
        @test _eq(quo(F, gens(SubQuos[1])), quo(F, SubQuos[1]))

        FHoms = [Oscar.FreeModuleHom_dec(F, F, [Hom_FreeModElems[t] for t = 1:3]), Oscar.FreeModuleHom_dec(F, F, [Hom_FreeModElems[t] for t = 4:6])]

        for t in [:none, :sum, :prod, :both]
          direct_product(F, F, task = t)
        end

        h = hom(F,F)
        for a = 1:2
          k = h[2].header.preimage(FHoms[a])
          @test _sparse_to_array(k, F) == [Hom_FreeModElems[t] for t = (3*a - 2) : (3*a)]
          for e in FreeModElems
            @test h[2].header.image(k)(e) == FHoms[a](e)
          end
        end

        Image = []
        w = [1, 3, 1]
        if j <= 3
          for s = 1:3
            deg = []
            for t = 1:3
              push!(deg, degree(Hom_SubQuoElems[s][t]).coeff[1] - degree(gen(SubQuos[w[s]], t)).coeff[1])
            end
            m = max(deg[1], deg[2], deg[3])
            push!(Image, [gen(F.R, 1)^(m - deg[t]) * Hom_SubQuoElems[s][t] for t = 1:3])
          end
        else
          for s = 1:3
            deg = []
            for t = 1:3
              #push!(deg, degree(Hom_SubQuoElems[s][t]).coeff - degree(gen(SubQuos[w[s]], t)).coeff)
              push!(deg, (iszero(Hom_SubQuoElems[s][t]) ? id(G).coeff : degree(Hom_SubQuoElems[s][t]).coeff) -
                    (iszero(gen(SubQuos[w[s]], t)) ? id(G).coeff : degree(gen(SubQuos[w[s]], t)).coeff))
            end
            m = [max(deg[1][l], deg[2][l], deg[3][l]) for l = 1:3]
            push!(Image, [((gen(F.R, 2) * gen(F.R, 3))^(m[1] - deg[t][1]) * gen(F.R, 2)^(m[2] - deg[t][2]) * (gen(F.R, 1) * gen(F.R, 2)^2)^(m[3] - deg[t][3])) * Hom_SubQuoElems[s][t]  for t = 1:3])
          end
        end
        non_zero = true
        for y = Image
          for x = y
            if iszero(x)
              non_zero = false
            end
          end
        end

        non_zero || continue #now the types and parents in Image are wrong

        w = [1, 2, 3]
        for t = 1:3
          Image[t] = gens(SubQuos[t])
        end

        SQHoms = [Oscar.SubQuoHom_dec(SubQuos[w[t]], SubQuos[t], Image[t]) for t = 1:3]
        for t in [:none, :prod, :sum]
          direct_product(SubQuos[4], SubQuos[4], task = t)
        end
			        
        #if i == 1
        #  @test ngens(direct_product(tensor_product(SubQuos[1], SubQuos[2]), tensor_product(SubQuos[1], SubQuos[3]))) == ngens(tensor_product(SubQuos[1], direct_product(SubQuos[2], SubQuos[3])))
        #  @test ngens(direct_product(tensor_product(F, F), tensor_product(F, F))[1]) == ngens(tensor_product(F, direct_product(F, F)[1]))
        #  tensor_product(SubQuos[1], SubQuos[3], task = :map)
	      #end
			        
        f = SQHoms[3]
        k = kernel(f)
        @test iszero(k[2](gen(k[1], 1)))
        im = image(f)
        @test _eq(im[1], sub(codomain(f), [im[2](x.a) for x = gens(im[1])]))
        D = homogeneous_components(f)
        first = true
        res = 0
        t = gen(domain(f), 1)
        for deg in keys(D)
          gm = D[deg]
          @test is_homogeneous(gm)
          @test degree(gm) == deg
          if first
            res = gm(t)
          else
            res += gm(t)
          end
          first = false
        end
        @test iszero(res - f(t))

        #=
        for Q in SubQuos
          free_resolution(Q)
        end
        =#

        I = ideal([hom_pols[t] for t = 1:3])
        R_quo = Oscar.MPolyQuo(RR, I)

        free_resolution(I)
        free_resolution(R_quo)
        free_resolution(F)

        for x in gens(F)
          @test ((FHoms[1] - FHoms[2]) * FHoms[1] * identity_map(F))(x) == ((FHoms[1] * FHoms[1]) - (FHoms[2] * FHoms[1]))(x)
          @test ((FHoms[1] + FHoms[2]) * FHoms[1])(x) == ((FHoms[1] * FHoms[1]) + (FHoms[2] * FHoms[1]) * identity_map(F))(x)
        end

        for f in FHoms
          @test _eq(image(f + f)[1], image(f)[1])
          D = homogeneous_components(f)
          for deg in keys(D)
            gm = D[deg]
            @test is_homogeneous(gm)
            @test degree(gm) == deg
            @test _eq(kernel(gm + gm)[1], kernel(gm)[1])
          end
        end
			        
        if R isa AbstractAlgebra.Field
          g = hom_keys[rand(1:6)]
          Ob = [F, SubQuos[3], I]
          El = [gen(F, rand(1:3)), Hom_SubQuoElems[3][rand(1:3)], hom_pols[rand(1:3)]]
          for t = 1:2
            comp = homogeneous_component(Ob[t], g)[2]

            #=
            if j < 4
              temp = homogeneous_component(El[t], g)
              comp.header.image(comp.header.preimage(temp)) == temp
            end
            =#

          end
        end

      end
    end
  end
end
