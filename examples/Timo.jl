function timo(f::fmpq_poly, p_all::Vector{Int})

  G, S = galois_group(f)
  K = Oscar.SolveRadical._fixed_field(S, [sub(G, [one(G)])[1]])

  zk = [pmaximal_overorder(equation_order(K.fld), p) for p in p_all]
  i = [prime_decomposition(zk[i], p_all[i])[1][1] for i in 1:length(p_all)]

  all_g = [x for x = G]
  @show :aut
  @time au = Oscar.SolveRadical.recognise(S, K, [K.pe^g for g = all_g])
#  au = Hecke.closure([hom(k.fld, k.fld, x) for x = au])
#recognise seems to be faster than closure - and keeps the link
#between the Galois group and the automorphisms
  @time au = [hom(K.fld, K.fld, x, check = false) for x = au]
  
  data = []
  for P = i
    @show "doin' ", minimum(P)
    F, mF = ResidueField(order(P), P)
    #the uniformizer is a generator for the local maximal order
    # over the unramifed subfield
    #the generator of the residue class field generates the
    # unramified maximal order (locally - but that is all than counts)
    gns = [K.fld(preimage(mF, gen(F))), K.fld(uniformizer(P))]
    # g in G is in the i-th decomposition group iff
    # forall a in gens, valuation(a-g(a), P) >= i+1
    # so, I can get all valuation groups at the price of one
    # the expensive bit is to compute valuation(a-g(a), P)
    #it appears that the REALLY hard bit is to apply the
    #automorphisms - so
    # - look into it
    # - combine the gens for the different prime ideals (CRT)
    #   to minimize this...
    g = au
    gg = all_g
    i = 0
    d = Dict{Int, Vector{Int}}()
    for i=1:length(gg)
      #= aut(x) -> eval of polynomials/ Q which takes forever
         replace by computation in Qp() of apropriate precision - much
         faster
         Also: the discriminant should give a upper bound on the valuation
         this can get the precision
      =#
      @time a = [x-g[i](x) for x = gns]
      v = [iszero(t) ? 10^6 : valuation(t, P) for t = a]
      m = minimum(v)
      if haskey(d, m)
        push!(d[m], i)
      else
        d[m] = [i]
      end
    end

    k = sort(collect(keys(d)))
    @assert k[end] == 10^6
    Gi = [G]
    Ki = Any[QQ]
    for i=0:k[end-1]
      push!(Gi, sub(G, gg[vcat([d[j] for j = k if j>=i+1]...)])[1])
      push!(Ki, fixed_field(S, Gi[end]))
    end
    push!(data, (Gi, Ki))
  end
  return data
end

#= the missing bit in Timo's code os the local factorisation
  but that should just be a coset length or so...
=#
