import AbstractAlgebra: MatElem
import Hecke: FqNmodFiniteField, FqNmodMatSpace, fq_nmod, fq_nmod_mat
import GAP: FFE

# isomorphism from the Oscar field to the GAP field
struct GenRingIso <: Map{FqNmodFiniteField,GapObj,SetMap,GenRingIso}
   domain::FqNmodFiniteField      # Oscar field
   codomain::GapObj          # GAP field
   f                      # function from the Oscar field to the GAP field
   f_inv               # its inverse
end

# isomorphism from the Oscar matrix space to the GAP matrix space
struct GenMatIso <: Map{FqNmodMatSpace,GapObj,SetMap,GenMatIso}
   domain::FqNmodMatSpace        # Oscar matrix space
   codomain::GapObj             # GAP matrix space
   fr::GenRingIso             # see type above
   f                         # function from the Oscar matrix space to the GAP matrix space
   f_inv                  # its inverse
end


################################################################################
#
#  Ring isomorphism
#
################################################################################

# return the GAP element corresponding to the Oscar element x
function elem_f(x::fq_nmod, B,d)
   COEF = [Int64(coeff(x,i)) for i in 0:d-1]
   x_gap = sum([COEF[i]*B[i] for i in 1:d])
   return x_gap
end

# return the Oscar element corresponding to the GAP element x
function elem_g(x::FFE, B, z, d)
   L = Vector{Int}(GAP.Globals.List(GAP.Globals.Coefficients(B,x),GAP.Globals.IntFFE ))
   return sum([L[i]*z^(i-1) for i in 1:d])
end


# computes the isomorphism between the Oscar field F and the corresponding GAP field F_gap
# the output has type GenRingIso (see above)

function gen_ring_iso(F::FqNmodFiniteField)
   p = Int64(characteristic(F))
   d = Int64(degree(F))
   z = gen(F)

   if d==1
      myOne = GAP.Globals.One(GAP.Globals.GF(p))
      f(x::fq_nmod) = Int(coeff(x,0))*myOne
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
   homom(x::fq_nmod) = elem_f(x,GAP.gap_to_julia(GAP.Globals.BasisVectors(Basis_F)),d)
   homominv(x::FFE) = elem_g(x, Basis_F, z, d)
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

# return the GAP matrix corresponding to the Oscar matrix x
# assumes x is a square matrix
function mat_oscar_gap(x::fq_nmod_mat, riso)
   n = nrows(x)
   S = Vector{GapObj}(undef, n)
   for i in 1:n
      S[i] = GAP.julia_to_gap([riso(x[i,j]) for j in 1:n])
   end

   return GAP.julia_to_gap(S)
end

# return the Oscar matrix corresponding to the GAP matrix x
# assumes x is a square matrix
function mat_gap_oscar(x::GapObj, riso)
   n = GAP.Globals.NrRows(x)
   Arr = [GAP.gap_to_julia(x[i]) for i in 1:n]
   L = [riso(Arr[i][j]) for i in 1:n for j in 1:n]

   return matrix(riso.domain, n, n, L)
end


# computes the isomorphism between the Oscar matrix space of dimension deg over F and the corresponding GAP matrix space
# the output has type GenMatIso (see above)

function gen_mat_iso(deg::Int, F::FqNmodFiniteField)
   riso = gen_ring_iso(F)                                      # "riso" = Ring ISOmorphism
   homom(x::fq_nmod_mat) = mat_oscar_gap(x, riso)
   homominv(x::GapObj) = mat_gap_oscar(x, riso)
   return GenMatIso(MatrixSpace(F,deg,deg),GAP.Globals.MatrixAlgebra(riso.codomain, deg), riso, homom, homominv)
end

(g::GenMatIso)(x::fq_nmod_mat) = g.f(x)
(g::GenMatIso)(x::GapObj) = g.f_inv(x)

Base.show(io::IO, f::GenMatIso) = print(io, "Matrix algebra homomorphism from Oscar algebra to GAP algebra")
