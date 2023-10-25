### Concrete implementation for tensor products of complexes

# As explained in `Types.jl`, we need to pass "factories" for the production of 
# the entries of a double complex and the morphisms between them to the 
# constructor. We first do the factory for the entries using the already 
# existing function for tensor products. 

# The concrete type contains all information necessary for the production.
mutable struct TensorProductFactory{ChainType} <: ChainFactory{ChainType}
  C1::ComplexOfMorphisms{<:ChainType}
  C2::ComplexOfMorphisms{<:ChainType}
end

# Then we override the call syntax as requested by the interface above.
function (fac::TensorProductFactory{<:ModuleFP})(dc::AbsDoubleComplexOfMorphisms, i::Int, j::Int)
  M = fac.C1[i]
  N = fac.C2[j]
  if iszero(M) || iszero(N)
    R = base_ring(M)
    is_graded(R) && return graded_free_module(R, 0)
    return FreeMod(R, 0)
  end
  return tensor_product(fac.C1[i], fac.C2[j])
end

# Now we move on to the factories for the maps. Here we have to provide 
# either one for the horizontal and the vertical morphisms.
mutable struct VerticalTensorMapFactory{MapType} <: ChainMorphismFactory{MapType}
  C1::ComplexOfMorphisms
  C2::ComplexOfMorphisms
end

# Again, we override the call syntax as required by the interface. Note that 
# here we need to make use of the first argument, the double complex itself, 
# in order to access the correct domains and codomains. 
function (fac::VerticalTensorMapFactory{<:ModuleFPHom})(dc::AbsDoubleComplexOfMorphisms, i::Int, j::Int)
  dom = dc[i, j]
  cod_ind = (i, j + (fac.C2.typ == :chain ? -1 : +1))
  cod = dc[cod_ind]
  (iszero(dom) || iszero(cod)) && return hom(dom, cod, elem_type(cod)[zero(cod) for i in 1:ngens(dom)])
  return tensor_product(dom, cod, [identity_map(fac.C1[i]), map(fac.C2, j)])
end

# Same for the horizontal maps.
mutable struct HorizontalTensorMapFactory{MapType} <: ChainMorphismFactory{MapType}
  C1::ComplexOfMorphisms
  C2::ComplexOfMorphisms
end

function (fac::HorizontalTensorMapFactory{<:ModuleFPHom})(dc::AbsDoubleComplexOfMorphisms, i::Int, j::Int)
  dom = dc[i, j]
  cod_ind = (i + (fac.C1.typ == :chain ? -1 : +1), j)
  cod = dc[cod_ind]
  (iszero(dom) || iszero(cod)) && return hom(dom, cod, elem_type(cod)[zero(cod) for i in 1:ngens(dom)])
  return tensor_product(dom, cod, [map(fac.C1, i), identity_map(fac.C2[j])])
end

# The user facing constructor for the tensor product of complexes 
# now takes the following rather simple form:
@doc raw"""
    tensor_product(C1::ComplexOfMorphisms{ChainType}, C2::ComplexOfMorphisms{ChainType}) where {ChainType}

Create the tensor product of two complexes `C1` and `C2` as a double complex.

In order for the generic implementation to work for your specific `ChainType` the following 
needs to be implemented.

  * `morphism_type(ChainType)` must produce the type of morphisms between objects of type `ChainType`;
  * the call signature for `function (fac::TensorProductFactory{ChainType})(dc::AbsDoubleComplexOfMorphisms, i::Int, j::Int)` needs to be overwritten for your specific instance of `ChainType` to produce the `(i, j)`-th entry of the double complex, i.e. the tensor product of `C1[i]` and `C2[j]`;
  * the call signature for `function (fac::HorizontalTensorMapFactory{ChainType})(dc::AbsDoubleComplexOfMorphisms, i::Int, j::Int)` needs to be overwritten to produce the map on tensor products `C1[i] ⊗ C2[j] → C1[i±1] ⊗ C2[j]` induced by the (co-)boundary map on `C1` (the sign depending on the `typ` of `C1`);
  * similarly for the `VerticalTensorMapFactory`.


See the file `experimental/DoubleComplex/src/tensor_products.jl` for examples.
"""
function tensor_product(C1::ComplexOfMorphisms{ChainType}, C2::ComplexOfMorphisms{ChainType}) where {ChainType}
  MorphismType = morphism_type(ChainType)
  result = DoubleComplexOfMorphisms(TensorProductFactory{ChainType}(C1, C2),
                                    HorizontalTensorMapFactory{MorphismType}(C1, C2),
                                    VerticalTensorMapFactory{MorphismType}(C1, C2),
                                    horizontal_typ=typ(C1),
                                    vertical_typ=typ(C2)
                                   )
  r1 = range(C1)
  result.right_bound = (typ(C1) == :chain ? first(r1) : last(r1))
  result.left_bound = (typ(C1) == :cochain ? first(r1) : last(r1))
  r2 = range(C2)
  result.upper_bound = (typ(C2) == :chain ? first(r2) : last(r2))
  result.lower_bound = (typ(C2) == :cochain ? first(r2) : last(r2))

  if is_complete(C1) || !isdefined(C1, :fill)
    result.extends_left = false
    result.extends_right = false
  elseif isdefined(C1, :fill) && typ(C1) == :chain # the filling function to extend to the right
    result.extends_right = true
    result.extends_left = false
  elseif isdefined(C1, :fill) && typ(C1) == :cochain # the filling function to extend to the right
    result.extends_right = false
    result.extends_left = true
  end

  if is_complete(C2) || !isdefined(C2, :fill)
    result.extends_down = false
    result.extends_up = false
  elseif isdefined(C2, :fill) && typ(C2) == :chain # the filling function to extend to the up
    result.extends_up = true
    result.extends_down = false
  elseif isdefined(C2, :fill) && typ(C2) == :cochain # the filling function to extend to the up
    result.extends_up = false
    result.extends_down = true
  end
  return result
end

### missing functionality
# A fallback in case nothing more specific is known
morphism_type(::Type{T}) where {T<:ModuleFP} = ModuleFPHom

zero(phi::ModuleFPHom) = hom(domain(phi), codomain(phi), [zero(codomain(phi)) for i in 1:ngens(domain(phi))])
zero_morphism(dom::ModuleFP, cod::ModuleFP) = hom(dom, cod, [zero(cod) for i in 1:ngens(dom)])
