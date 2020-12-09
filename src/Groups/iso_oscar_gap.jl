import AbstractAlgebra: MatElem
import Hecke: FqNmodFiniteField, FqNmodMatSpace, fq_nmod, fq_nmod_mat
import GAP: FFE

# isomorphism from the Oscar field to the GAP field
struct GenRingIso <: Map{FqNmodFiniteField,GapObj,SetMap,GenRingIso}
   domain::FqNmodFiniteField
   codomain::GapObj
   f
   f_inv
end

# turns the Oscar matrix into the GAP matrix and viceversa
struct GenMatIso <: Map{FqNmodMatSpace,GapObj,SetMap,GenMatIso}
   domain::FqNmodMatSpace
   codomain::GapObj
   fr::GenRingIso
   f
   f_inv
end


################################################################################
#
#  Ring isomorphism
#
################################################################################

function elem_f(x::fq_nmod; B=0,d=0)
   COEF = [Int64(coeff(x,i)) for i in 0:d-1]
   x_gap = sum([COEF[i]*B[i] for i in 1:d])
   return x_gap
end

function elem_g(x::FFE; B=0, z=0, d=0)
   L = Vector{Int}(GAP.Globals.List(GAP.Globals.Coefficients(B,x),GAP.Globals.IntFFE ))
   return sum([L[i]*z^(i-1) for i in 1:d])
end

function gen_ring_iso(F::FqNmodFiniteField)
   p = Int64(characteristic(F))
   d = Int64(degree(F))
   z = gen(F)

   if d==1
      f(x::fq_nmod) = Int(coeff(x,0))*GAP.Globals.One(GAP.Globals.GF(p))
      finv(x::FFE) = F(GAP.Globals.IntFFE(x))
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

Base.show(io::IO, f::GenRingIso) = print(io, "Ring isomorphism between ", f.domain, " and the corresponding GAP")



################################################################################
#
#  Matrix space isomorphism
#
################################################################################

function mat_oscar_gap(x::fq_nmod_mat; n=0, r=0)
   S = Vector{GapObj}(undef, n)
   for i in 1:n
      S[i] = GAP.julia_to_gap([r(x[i,j]) for j in 1:n])
   end

   return GAP.julia_to_gap(S)
end

function mat_gap_oscar(x::GapObj; n=0, r=0)
   Arr = [GAP.gap_to_julia(x[i]) for i in 1:n]
   L = [r(Arr[i][j]) for i in 1:n for j in 1:n]

   return matrix(r.domain, n, n, L)
end


function gen_mat_iso(deg::Int, F::FqNmodFiniteField)
   riso = gen_ring_iso(F)                                      # "riso" = Ring ISOmorphism
   homom(x::fq_nmod_mat) = mat_oscar_gap(x; n=deg, r=riso)
   homominv(x::GapObj) = mat_gap_oscar(x; n=deg, r=riso)
   return GenMatIso(MatrixSpace(F,deg,deg),GAP.Globals.MatrixAlgebra(riso.codomain, deg), riso, homom, homominv)
end

(g::GenMatIso)(x::fq_nmod_mat) = g.f(x)
(g::GenMatIso)(x::GapObj) = g.f_inv(x)

Base.show(io::IO, f::GenMatIso) = print(io, "Matrix algebra homomorphism from Oscar algebra to GAP algebra")
