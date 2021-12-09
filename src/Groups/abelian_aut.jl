export _isomorphic_gap_group

function _isomorphic_gap_group(A::GrpAbFinGen; T=PcGroup)

  eldiv = [Int(a) for a in elementary_divisors(A)]
  Agap = abelian_group(T, eldiv)
  Asnf, Asnf_to_A = snf(A)
  A_to_Asnf = inv(Asnf_to_A)

  function A_to_Agap(a)
      return _oscar_to_gap(A_to_Asnf(a), Agap)
  end
  G = GAP.Globals
  gensindep = G.IndependentGeneratorsOfAbelianGroup(Agap.X)
  Aindep = abelian_group([G.Order(g) for g in gensindep])
  exp = [GAP.gap_to_julia(Vector{fmpz},G.IndependentGeneratorExponents(Agap.X, a.X)) for a in gens(Agap)]
  Asnf_to_Aindep = hom(Asnf, Aindep, [Aindep(e) for e in exp])
  Aindep_to_Asnf = inv(Asnf_to_Aindep)
  Aindep_to_A = compose(Aindep_to_Asnf, Asnf_to_A)

  function Agap_to_A(a)
      return Aindep_to_A(_gap_to_oscar(a, Aindep))
  end

  to_gap = Hecke.map_from_func(A_to_Agap, A, Agap)
  to_oscar = Hecke.map_from_func(Agap_to_A, Agap, A)

  return Agap, to_gap, to_oscar
end

function _oscar_to_gap(a::GrpAbFinGenElem, B::GAPGroup)
  gensB = gens(B)
  n = length(a.coeff)
  @assert length(gensB) == n
  img = one(B)
  for i in 1:n
    img *= gensB[i]^a[i]
  end
  return img
end

function _gap_to_oscar(a::Oscar.BasicGAPGroupElem, B::GrpAbFinGen)
  A = parent(a)
  G = GAP.Globals
  exp = GAP.gap_to_julia(Vector{fmpz},G.IndependentGeneratorExponents(A.X, a.X))
  return B(exp)
end

(f::GAPGroupElem{AutomorphismGroup{GrpAbFinGen}})(x::GrpAbFinGenElem)  = apply_automorphism(f, x, true)
Base.:^(x::GrpAbFinGenElem,f::GAPGroupElem{AutomorphismGroup{GrpAbFinGen}}) = apply_automorphism(f, x, true)

function apply_automorphism(f::GAPGroupElem{AutomorphismGroup{GrpAbFinGen}}, x::GrpAbFinGenElem, check=true)
  aut = parent(f)
  if check
    @assert parent(x) == aut.G "Not in the domain of f!"
  end
  xgap = aut.to_gap(x)
  A = parent(f)
  domGap = parent(xgap)
  imgap = typeof(xgap)(domGap, GAP.Globals.Image(f.X,xgap.X))
  return aut.to_oscar(imgap)
end

function automorphism_group(G::GrpAbFinGen)
  Ggap, to_gap, to_oscar = _isomorphic_gap_group(G)
  AutGAP = GAP.Globals.AutomorphismGroup(Ggap.X)
  return AutomorphismGroup{typeof(G)}(AutGAP, G, to_gap, to_oscar)
end




