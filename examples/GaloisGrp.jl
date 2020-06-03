module GaloisGrp

using Oscar
import Base: ^
import Hecke

function ^(f::MPolyElem, s::Oscar.GAPGroupElem{PermGroup})
  G = parent(s)
  @assert ngens(parent(f)) == degree(G)

  g = Generic.MPolyBuildCtx(parent(f))
  for (c, e) = Base.Iterators.zip(Generic.MPolyCoeffs(f), Generic.MPolyExponentVectors(f))
    s_e = zeros(Int, degree(G))
    for i=1:degree(G)
      s_e[s(i)] = e[i]
    end
    push_term!(g, c, s_e)
  end
  return finish(g)
end

function orbit(G::Oscar.PermGroup, f::MPolyElem)
  s = Set([f])
  while true
    n = Set(x^g for x = s for g = gens(G))
    sn = length(s)
    union!(s, n)
    if length(s) == sn
      break
    end
  end
  return s
end

function maximal_subgroup_reps(G::PermGroup)
  return Oscar._as_subgroups(GAP.Globals.MaximalSubgroupClassReps(G.X), G)
end

function galois_group(K::AnticNumberField, p::Int = 31)
  C = qAdicConj(K, p, splitting_field = true)
  G = symmetric_group(degree(K))
  Zxx, x = PolynomialRing(ZZ, :x=>1:degree(K))
  m = prod(gen(Zxx, i)^i for i = 1:degree(K)-1)

  Zx, x = PolynomialRing(ZZ)

  c = conjugates(gen(K), C, 10)

  nG = length(G)

  while true
    @show S = maximal_subgroup_reps(G)
    for s = S
      @show s
      @show I = sum(orbit(s, m))
      @assert all(I^x == I for x = s)
      @assert !any(I^x == I for x = right_transversal(G, s) if !isone(x))
      local fd
      cnt = 0
      while true
        cs = Set{typeof(c[1])}()
        fd = []
        for t = right_transversal(G, s)
          e = evaluate(I^t, c)
          if e in cs
            @show "darn, need tschirni", e, cs
            local ts
            c = conjugates(gen(K), C, 10)
            while true
              @show ts = rand(Zx, 2:degree(K), -2:2)
              if degree(ts) > 1
                break
              end
            end
            cnt += 1
            if cnt > 2
              error("bad")
            end
            c = map(ts, c)
            break
          end
          push!(cs, e)
          if e.length<2
            l = lift(coeff(e, 0))
            l = Hecke.mod_sym(l, fmpz(p)^10)
            if abs(l) < fmpz(p)^5
              println(t, " => ", l)
              push!(fd, t)
            end
          end
        end
        @show length(cs), index(G, s)
        if length(cs) == index(G, s)
          break
        end
      end
      @show "descent via", fd, s
      if length(fd)>0
        G = intersection([s^x for x = fd]...)[1]
        break
      end
    end
    if length(G) == nG
      return G
    else
      nG = length(G)
    end
  end
end

end
