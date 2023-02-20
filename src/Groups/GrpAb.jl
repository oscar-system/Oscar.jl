# Something additional functions for abelian groups

# Restrict a morphism f : G -> H to the surjective morphism g : G -> im(G)
function restrict_codomain(f::GrpAbFinGenMap)
  G = domain(f)
  H, Htocd = image(f, false)
  imgs = elem_type(H)[]
  for g in gens(G)
    fl, h = haspreimage(Htocd, f(g))
    @assert fl
    push!(imgs, h)
  end
  return hom(G, H, imgs)
end
