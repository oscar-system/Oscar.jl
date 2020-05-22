import Hecke: order, base_ring, elements
import GAP: FFE


# TODO: at this moment, only finite fields are implemented; in future we want to use infinite fields
const TempMatType = fq_nmod_mat

function FieldGapToHecke(F::GapObj)
   p = GAP.Globals.Characteristic(F)
   q = GAP.Globals.Size(F)
   d = Base.Integer(log(p,q))

   return GF(p,d)
end

function FieldHeckeToGap(F::FqNmodFiniteField)
   p = Int64(characteristic(F))
   d = Int64(degree(F))

   return GAP.Globals.GF(p,d)
end

function FieldElemGapToHecke(x::FFE, F::GapObj)
   q = GAP.Globals.Size(F)
   K,z = FieldGapToHecke(F)
   if GAP.Globals.IsZero(x)
      return 0*z
   else
      d = Integer(GAP.Globals.LogFFE(x, GAP.Globals.Z(q)))
      return z^d
   end
end

function FieldElemHeckeToGap(x::fq_nmod, F::FqNmodFiniteField)
   q = Int64(order(F))
   d = degree(F)
   z = gen(F)
   v = [Int64(coeff(x,i)) for i in 0:d]
   y = sum([GAP.Globals.Z(q)^i*v[i+1] for i in 0:d])
   
   return y
end

function MatGapToHecke(x::GapObj, F::GapObj)
   n = GAP.Globals.Size(x)
   Arr = [GAP.gap_to_julia(x[i]) for i in 1:n]
   L = [FieldElemGapToHecke(Arr[i][j],F) for i in 1:n for j in 1:n]

   return matrix(FieldGapToHecke(F)[1],n,n,L)
end

function MatHeckeToGap(x::TempMatType, F::FqNmodFiniteField)
   r = nrows(x)
   c = ncols(x)
   S = Vector{GapObj}(undef,r)
   for i in 1:r
      S[i] = GAP.julia_to_gap([FieldElemHeckeToGap(x[i,j],F) for j in 1:c])
   end

   return GAP.julia_to_gap(S)
end

function MatJuliaToHecke(x::GAPGroupElem{MatrixGroup})
   F = GAP.Globals.FieldOfMatrixGroup(parent(x).X)
   return MatGapToHecke(x.X,F)
end

function Base.getindex(x::GAPGroupElem{MatrixGroup}, i::Int64, j::Int64)
   F = GAP.Globals.FieldOfMatrixGroup(parent(x).X)
   return FieldElemGapToHecke((x.X)[i,j], F)
end

function order(x::TempMatType)
   return GAP.Globals.Order(MatHeckeToGap(x,base_ring(parent(x))))
end

function order(x::fq_nmod)
   return GAP.Globals.Order(FieldElemHeckeToGap(x,parent(x)))
end

Base.show(io::IO, x::GAPGroupElem{MatrixGroup}) = print(io, MatJuliaToHecke(x))

function (G::MatrixGroup)(x::TempMatType)
   F = base_ring(x)
   y = MatHeckeToGap(x,F)
   @assert GAP.Globals.IN(y, G.X) "x does not belong to the group G"
   return group_element(G,MatHeckeToGap(x,F))
end

function base_ring(G::MatrixGroup)
   return FieldGapToHecke(GAP.Globals.FieldOfMatrixGroup(G.X))[1]
end

base_ring(x::GAPGroupElem{MatrixGroup}) = base_ring(parent(x))

########################################################################
#
# functions involving matrices
#
########################################################################

function Base.:in(x::TempMatType, G::MatrixGroup)
   F = base_ring(G)
   return GAP.Globals.IN(MatHeckeToGap(x,F),G.X) 
end

Base.:*(x::GAPGroupElem{MatrixGroup}, y::TempMatType) = MatJuliaToHecke(x)*y
Base.:*(x::TempMatType, y::GAPGroupElem{MatrixGroup}) = x*MatJuliaToHecke(y)

Base.:^(x::TempMatType, y::TempMatType) = inv(y)*x*y
Base.:^(x::GAPGroupElem{MatrixGroup}, y::TempMatType) = parent(x)(inv(y)*x*y)
Base.:^(x::TempMatType, y::GAPGroupElem{MatrixGroup}) = inv(y)*x*y

comm(x::TempMatType, y::TempMatType) = inv(x)*x^y
comm(x::GAPGroupElem{MatrixGroup}, y::TempMatType) = inv(x)*x^y
comm(x::TempMatType, y::GAPGroupElem{MatrixGroup}) = inv(x)*x^y

function conjugacy_class(G::MatrixGroup, g::TempMatType)
   F = base_ring(g)
   x = group_element(G,MatHeckeToGap(g,F))
   return _conjugacy_class(G, x, GAP.Globals.ConjugacyClass(G.X,x.X))
end

function isconjugate(G::MatrixGroup, x::Union{GAPGroupElem{MatrixGroup},TempMatType}, y::Union{GAPGroupElem{MatrixGroup},TempMatType})
   F = base_ring(G)
   if typeof(x)==TempMatType xgap = MatHeckeToGap(x,F) else xgap = x.X end
   if typeof(y)==TempMatType ygap = MatHeckeToGap(y,F) else ygap = y.X end
   if GAP.Globals.IsConjugate(G.X, xgap, ygap)
      K = GAP.Globals.FieldOfMatrixGroup(G.X)
      return true, group_element(G,GAP.Globals.RepresentativeAction(G.X, xgap, ygap))
   else
      return false, nothing
   end
end

function normalizer(G::MatrixGroup, x::TempMatType)
   F = base_ring(G)
   xgap = MatHeckeToGap(x,F)
   _as_subgroup(G, GAP.Globals.Normalizer(G.X,xgap))
end

function centralizer(G::MatrixGroup, x::TempMatType)
   F = base_ring(G)
   xgap = MatHeckeToGap(x,F)
   _as_subgroup(G, GAP.Globals.Centralizer(G.X,xgap))
end

function sub(G::MatrixGroup, L::Vector{TempMatType})
   F = base_ring(G)
   l = Vector{GapObj}(undef, length(L))
   for i in 1:length(L)
      l[i] = MatHeckeToGap(L[i],F)
   end
   elems_in_GAP = GAP.julia_to_gap(l)
   H = GAP.Globals.Group(elems_in_GAP)
   #H is the group. I need to return the inclusion map too
   return _as_subgroup(G, H)
end

function sub(L::TempMatType...)
   if length(L)==0 throw(ArgumentError("Empty list")) end
   F = base_ring(parent(L[1]))
   l=collect(L)
#   @assert all(x -> parent(x) == parent(l[1]), l)
   return sub(GL(nrows(L[1]),Int64(order(F))),l)
end

function quo(G::MatrixGroup, L::Vector{TempMatType})
   F = base_ring(G)
   l = Vector{GapObj}(undef, length(L))
   for i in 1:length(L)
      l[i] = MatHeckeToGap(L[i],F)
   end
   elems_in_gap = GAP.julia_to_gap(l)
   H = GAP.Globals.NormalClosure(G.X,GAP.Globals.Group(elems_in_gap))
   @assert GAP.Globals.IsNormal(G.X, H)
   H1 = MatrixGroup(H)
   return quo(G, H1)
end

function Base.:^(H::MatrixGroup, y::TempMatType)
   F = base_ring(parent(y))
   return MatrixGroup(H.X ^ MatHeckeToGap(y,F))
end

#########################################################################
#
# homomorphisms
#
#########################################################################

#TODO is it really necessary to write this function three times?

function hom(G::MatrixGroup, H::GAPGroup, gensG::Vector{TempMatType}, imgs::Vector)
  F = base_ring(G)
  vgens = GAP.julia_to_gap(GapObj[MatHeckeToGap(x,F) for x in gensG])
  vimgs = GAP.julia_to_gap(GapObj[x.X for x in imgs])
  mp = GAP.Globals.GroupHomomorphismByImages(G.X, H.X, vgens, vimgs)
  return GAPGroupHomomorphism{MatrixGroup, typeof(H)}(G, H, mp)
end

function hom(G::GAPGroup, H::MatrixGroup, gensG::Vector, imgs::Vector{TempMatType})
  F = base_ring(H)
  vgens = GAP.julia_to_gap(GapObj[x.X for x in gensG])
  vimgs = GAP.julia_to_gap(GapObj[MatHeckeToGap(x,F) for x in imgs])
  mp = GAP.Globals.GroupHomomorphismByImages(G.X, H.X, vgens, vimgs)
  return GAPGroupHomomorphism{typeof(G), MatrixGroup}(G, H, mp)
end

function hom(G::MatrixGroup, H::MatrixGroup, gensG::Vector{TempMatType}, imgs::Vector{TempMatType})
  F = base_ring(G)
  K = base_ring(H)
  vgens = GAP.julia_to_gap(GapObj[MatHeckeToGap(x,F) for x in gensG])
  vimgs = GAP.julia_to_gap(GapObj[MatHeckeToGap(x,K) for x in imgs])
  mp = GAP.Globals.GroupHomomorphismByImages(G.X, H.X, vgens, vimgs)
  return GAPGroupHomomorphism{MatrixGroup, typeof(H)}(G, H, mp)
end

#TODO is it really necessary to write the function two times?
function image(f::GAPGroupHomomorphism{MatrixGroup,T}, x::TempMatType) where T <: GAPGroup
   F = base_ring(domain(f))
   y = MatHeckeToGap(x,F)
   return group_element(codomain(f), GAP.Globals.Image(f.map,y))
end

function image(f::GAPGroupHomomorphism{MatrixGroup,MatrixGroup}, x::TempMatType)
   F = base_ring(domain(f))
   y = MatHeckeToGap(x,F)
   K = GAP.Globals.FieldOfMatrixGroup(codomain(f).X)
   return group_element(codomain(f),GAP.Globals.Image(f.map,y))
end

(f::GAPGroupHomomorphism{MatrixGroup,T})(x::TempMatType) where T <: GAPGroup = image(f, x)
(f::GAPGroupHomomorphism{MatrixGroup,MatrixGroup})(x::TempMatType) = image(f,x)
Base.:^(x::TempMatType,f::GAPGroupHomomorphism) = image(f,x)

function haspreimage(f::GAPGroupHomomorphism{T,MatrixGroup}, x::TempMatType) where T <: GAPGroup
  F = base_ring(codomain(f))
  r = GAP.Globals.PreImagesRepresentative(f.map, MatHeckeToGap(x,F))
  if r == GAP.Globals.fail
    return false, one(domain(f))
  else
    return true, group_element(domain(f), r)
  end
end

(f::GAPGroupElem{AutomorphismGroup{MatrixGroup}})(x::TempMatType) = apply_automorphism(f, x)
Base.:^(x::TempMatType,f::GAPGroupElem{AutomorphismGroup{MatrixGroup}}) = apply_automorphism(f, x)

function apply_automorphism(f::GAPGroupElem{AutomorphismGroup{MatrixGroup}}, x::TempMatType, check=true)
  A = parent(f)
  G = A.G
  F = base_ring(G)
  return group_element(G, GAP.Globals.Image(f.X,MatHeckeToGap(x,F)))
 # return MatGapToHecke(GAP.Globals.Image(f.X,MatHeckeToGap(x,F)),GAP.Globals.FieldOfMatrixGroup(G.X))
end
