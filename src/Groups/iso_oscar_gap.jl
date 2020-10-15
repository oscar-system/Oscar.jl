
# isomorphism from the Oscar field to the GAP field
mutable struct GenRingIso <: Map{FqNmodFiniteField,GapObj,SetMap,GenRingIso}
   domain::FqNmodFiniteField
   codomain::GapObj
   f
   f_inv

   function GenRingIso(F::FqNmodFiniteField, Fg::GapObj, g, ginv)
      z = new()
      z.domain = F
      z.codomain = Fg
      z.f = g
      z.f_inv = ginv
      return z
   end
end


function elem_f(x::fq_nmod; B=0,d=0)
   COEF = [Int64(coeff(x,i)) for i in 0:d-1]
   x_gap = sum([COEF[i]*B[i] for i in 1:d])
   return x_gap
end

function elem_g(x::FFE; B=0, z=0, d=0)
   L = GAP.gap_to_julia(GAP.Globals.List(GAP.Globals.Coefficients(B,x), y-> GAP.Globals.IntFFE(y)))
   return sum([L[i]*z^(i-1) for i in 1:d])
end

function gen_ring_iso(F::FqNmodFiniteField)
   p = Int64(characteristic(F))
   d = Int64(degree(F))
   z = gen(F)

   if d==1
      f(x::fq_nmod) = x -> Int(coeff(x,0))*GAP.Globals.One(GAP.Globals.GF(p))
      finv(x::FFE) = x -> F(GAP.Globals.IntFFE(x))
      return GenRingIso(F, GAP.Globals.GF(p), f, finv)
   end
   L = [Int64(lift(coeff(defining_polynomial(F),i))) for i in 0:d-1]
   L = vcat(L,[1])
   i = GAP.Globals.One(GAP.Globals.GF(p))
   L_gap = GAP.julia_to_gap([i*y for y in L])
   f_gap = GAP.Globals.UnivariatePolynomial(GAP.Globals.GF(p),L_gap)
   F_gap = GAP.Globals.GF(GAP.Globals.GF(p),f_gap)
   Basis_F = GAP.Globals.Basis(F_gap)
   homom(x::fq_nmod) = elem_f(x;B=GAP.gap_to_julia(GAP.Globals.BasisVectors(Basis_F)),d=d)
   homominv(x::FFE) = elem_g(x; B=Basis_F, z=z, d=d)
   return GenRingIso(F, F_gap, homom, homominv)
end

(g::GenRingIso)(x::fq_nmod) = g.f(x)
(g::GenRingIso)(x::FFE) = g.f_inv(x)

Base.show(io::IO, f::GenRingIso) = print(io, "Ring isomorphism between ", F, " and the corresponding GAP")

