module GModuleFromGap
using Oscar
using Hecke

import Oscar:gmodule

import AbstractAlgebra: Group, Module
import Base: parent

function GAP.gap_to_julia(::Type{QabElem}, a::GAP.GapObj) #which should be a Cyclotomic
  c = GAP.Globals.Conductor(a)
  E = abelian_closure(QQ)[2](c)
  z = parent(E)(0)
  co = GAP.Globals.CoeffsCyc(a, c)
  for i=1:c
    if !iszero(co[i])
      z += fmpq(co[i])*E^(i-1)
    end
  end
  return z
end

(::QabField)(a::GAP.GapObj) = GAP.gap_to_julia(QabElem, a)

function irreducible_modules(G::Oscar.GAPGroup)
  im = GAP.Globals.IrreducibleRepresentations(G.X)
  IM = GModule[] 
  K = abelian_closure(QQ)[1]
  for m in im
    z = map(x->matrix(map(y->map(K, y), m(x.X))), gens(G))
    F = free_module(K, nrows(z[1]))
    push!(IM, gmodule(G, map(x->hom(F, F, x), z)))
  end
  return IM
end

function minimize(a::AbstractArray{nf_elem})
  fl, c = Hecke.iscyclotomic_type(parent(a[1]))
  @assert fl
  for p = keys(factor(c).fac)
    while c % p == 0
      K, _ = cyclotomic_field(Int(div(c, p)), cached = false)
      b = similar(a)
      OK = true
      for x = eachindex(a)
        y = Hecke.force_coerce_cyclo(K, a[x], Val{false})
        if y == false
          OK = false
        else
          b[x] = y
        end
      end
      if OK
        a = b
        c = div(c, p)
      else
        break
      end
    end
  end
  return a
end

function minimize(a::MatElem{nf_elem})
  return matrix(minimize(a.entries))
end

function minimize(a::nf_elem)
  return minimize([a])[1]
end

function Oscar.conductor(a::nf_elem)
  return conductor(parent(minimize(a)))
end

function Oscar.conductor(a::QabElem)
  return conductor(data(a))
end

function irreducible_modules(::Type{AnticNumberField}, G::Oscar.GAPGroup)
  z = irreducible_modules(G)
  Z = GModule[]
  for m in z
  end
end

function gmodule(::typeof(CyclotomicField), C::GModule)
  @assert isa(base_ring(C), QabField)
  d = dim(C)
  l = 1
  for g = C.ac
    l = lcm(l, lcm(collect(map_entries(x->Hecke.iscyclotomic_type(parent(x.data))[2], mat(g)))))
  end
  K = cyclotomic_field(l, cached = false)[1]
  F = free_module(K, dim(C))
  return gmodule(group(C), [hom(F, F, map_entries(x->K(x.data), mat(x))) for x = C.ac])
end

import Base: ^
function ^(C::GModule{<:Any, Generic.FreeModule{nf_elem}}, phi::Map{AnticNumberField, AnticNumberField})
  F = free_module(codomain(phi), dim(C))
  return GModule(group(C), [hom(F, F, map_entries(phi, mat(x))) for x = C.ac])
end

function ^(C::GModule{<:Any, Generic.FreeModule{QabElem}}, phi::Map{QabField, QabField})
  F = free_module(codomain(phi), dim(C))
  return GModule(group(C), [hom(F, F, map_entries(phi, mat(x))) for x = C.ac])
end

function gmodule(::FlintRationalField, C::GModule{<:Any, Generic.FreeModule{nf_elem}})
  F = free_module(QQ, dim(C)*degree(base_ring(C)))
  return GModule(group(C), [hom(F, F, hvcat(dim(C), [representation_matrix(x) for x = transpose(mat(y))]...)) for y = C.ac])
end

function gmodule(k::Nemo.GaloisField, C::GModule{<:Any, Generic.FreeModule{fmpq}})
  F = free_module(k, dim(C))
  return GModule(group(C), [hom(F, F, map_entries(k, mat(x))) for x=C.ac])
end

function Gap(C::GModule{<:Any, Generic.FreeModule{gfp_elem}})
  z = AbstractAlgebra.get_special(C, :Gap)
  if z !== nothing
    return z
  end
  k = GAP.Globals.Z(Int(characteristic(base_ring(C))))
  one = k^0
  z = GAP.Globals.GModuleByMats(GAP.julia_to_gap([GAP.julia_to_gap(map(x->Int(x)*k, Matrix(lift(mat(x))))) for x = C.ac]), GAP.Globals.GF(Int(characteristic(base_ring(C)))))
  AbstractAlgebra.set_special(C, :Gap=>z)
  return z
end

function Oscar.isirreducible(C::GModule{<:Any, Generic.FreeModule{gfp_elem}})
  G = Gap(C)
  return GAP.Globals.MTX.IsIrreducible(G)
end

function isabsolutely_irreducible(C::GModule{<:Any, Generic.FreeModule{gfp_elem}})
  G = Gap(C)
  return GAP.Globals.MTX.IsAbsolutelyIrreducible(G)
end

function isdecomposable(C::GModule{<:Any, Generic.FreeModule{gfp_elem}})
  G = Gap(C)
  return !GAP.Globals.MTX.IsIndecomposable(G)
end

function Oscar.hom(C::T, D::T) where T <: GModule{<:Any, Generic.FreeModule{gfp_elem}}
  @assert base_ring(C) == base_ring(D)
  return GAP.Globals.MTX.BasisModuleHomomorphisms(Gap(C), Gap(D))
end

export irreducible_modules, isabsolutely_irreducible, isdecomposable

end #module GModuleFromGap

using .GModuleFromGap

export irreducible_modules, isabsolutely_irreducible, isdecomposable

