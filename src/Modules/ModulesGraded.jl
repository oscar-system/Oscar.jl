export decoration, free_module_dec


###############################################################################
# FreeMod_dec constructors
###############################################################################

@doc Markdown.doc"""
    FreeMod_dec(R::CRing_dec, n::Int, name::String = "e"; cached::Bool = false) 

Construct a decorated (graded or filtered) free module over the ring `R` with rank `n`
with the standard degrees, that is the standard unit vectors have degree 0.
Additionally one can provide names for the generators. If one does 
not provide names for the generators, the standard names e_i are used for 
the standard unit vectors.
"""
function FreeMod_dec(R::CRing_dec, n::Int, name::String = "e"; cached::Bool = false) 
  return FreeMod_dec{elem_type(R)}(R, [Symbol("$name[$i]") for i=1:n], [decoration(R)[0] for i=1:n])
end

@doc Markdown.doc"""
    free_module_dec(R::CRing_dec, n::Int, name::String = "e"; cached::Bool = false)

Create the decorated free module $R^n$ equipped with its basis of standard unit vectors
and standard degrees, that is the standard unit vectors have degree 0.

The string `name` specifies how the basis vectors are printed. 

# Examples
```jldoctest
julia> R, (x,y) = grade(PolynomialRing(QQ, ["x", "y"])[1])
(Multivariate Polynomial Ring in x, y over Rational Field graded by 
  x -> [1]
  y -> [1], MPolyElem_dec{fmpq, fmpq_mpoly}[x, y])

julia> free_module_dec(R,3)
Decorated free module of rank 3 over Multivariate Polynomial Ring in x, y over Rational Field graded by 
  x -> [1]
  y -> [1]Multivariate Polynomial Ring in x, y over Rational Field graded by 
  x -> [1]
  y -> [1]^3([0])

```
"""
free_module_dec(R::CRing_dec, n::Int, name::String = "e"; cached::Bool = false) = FreeMod_dec(R, n, name, cached = cached)


@doc Markdown.doc"""
    FreeMod_dec(R::CRing_dec, d::Vector{GrpAbFinGenElem}, name::String = "e"; cached::Bool = false) 

Construct a decorated (graded or filtered) free module over the ring `R` 
with rank `n` where `n` is the length of `d`. `d` is the vector of degrees for the 
components, i.e. `d[i]` is the degree of `e[i]` where `e[i]` is the `i`th standard unit
vector of the free module.
Additionally one can provide names for the generators. If one does 
not provide names for the generators, the standard names e_i are used for 
the standard unit vectors.
"""
function FreeMod_dec(R::CRing_dec, d::Vector{GrpAbFinGenElem}, name::String = "e"; cached::Bool = false) 
  return FreeMod_dec{elem_type(R)}(R, [Symbol("$name[$i]") for i=1:length(d)],d)
end

@doc Markdown.doc"""
    free_module_dec(R::CRing_dec, d::Vector{GrpAbFinGenElem}, name::String = "e"; cached::Bool = false)

Create the decorated free module $R^n$ (`n` is the length of `d`)
equipped with its basis of standard unit vectors where the 
i-th standard unit vector has degree `d[i]`.

The string `name` specifies how the basis vectors are printed. 
"""
free_module_dec(R::CRing_dec, d::Vector{GrpAbFinGenElem}, name::String = "e"; cached::Bool = false) = FreeMod_dec(R, d, name, cached = cached)


function FreeMod_dec(F::FreeMod, d::Vector{GrpAbFinGenElem})
  return FreeMod_dec{elem_type(base_ring(F))}(F, d)
end


function AbstractAlgebra.extra_name(F::FreeMod_dec)
  t = get_attribute(F, :twist)
  if t !== nothing
    n = get_attribute(t[1], :name)
    if n !== nothing
      return "$n($(t[2]))"
    end
  end
  if length(Set(F.d)) == 1
    n = get_attribute(forget_decoration(F).R, :name)
    if n !== nothing
      return "$n^$(ngens(F))($(-F.d[1]))"
    end
  end
  return nothing
end

function show(io::IO, F::FreeMod_dec)
  @show_name(io, F)
  @show_special(io, F)

  print(io, "Decorated free module of rank $(rank(F)) over ")
  print(IOContext(io, :compact =>true), base_ring(F))

  i = 1
  while i < dim(F)
    d = F.d[i]
    j = 1
    while i+j <= dim(F) && d == F.d[i+j]
      j += 1
    end
    print(IOContext(io, :compact => true), base_ring(F), "^$j")
    print(IOContext(io, :compact => true), "(", -d, ")")
    if i+j < dim(F)
      print(io, " + ")
    end
    i += j
  end
end

function forget_decoration(F::FreeMod_dec)
  return F.F
end

@doc Markdown.doc"""
    base_ring(F::FreeMod_dec)

Return the underlying ring of `F`.
"""
base_ring(F::FreeMod_dec) = forget_decoration(F).R

@doc Markdown.doc"""
    rank(F::FreeMod_dec)

Return the rank of `F`.
"""
rank(F::FreeMod_dec) = rank(forget_decoration(F))

@doc Markdown.doc"""
    decoration(F::FreeMod_dec)

Return the vector of degrees of the standard unit vectors.
"""
decoration(F::FreeMod_dec) = F.d
decoration(R::MPolyRing_dec) = R.D

@doc Markdown.doc"""
    is_graded(F::FreeMod_dec)

Check if `F` is graded.
"""
is_graded(F::FreeMod_dec) = is_graded(base_ring(F))

@doc Markdown.doc"""
    is_filtered(F::FreeMod_dec)

Check if `F` is filtered.
"""
is_filtered(F::FreeMod_dec) = is_filtered(base_ring(F))

is_decorated(F::FreeMod_dec) = true

@doc Markdown.doc"""
    ==(F::FreeMod_dec, G::FreeMod_dec)

Return  `true` if `F` and `G` are equal, `false` otherwise.

Here, `F` and `G` are equal iff their base rings, ranks, decorations 
and names for printing the basis elements are equal.
"""
function Base.:(==)(F::FreeMod_dec, G::FreeMod_dec)
  return forget_decoration(F) == forget_decoration(G) && F.d == G.d
end

###############################################################################
# FreeModElem_dec constructors
###############################################################################

@doc Markdown.doc"""
    FreeModElem_dec(c::SRow{T}, parent::FreeMod_dec{T}) where T

Return the element of `F` whose coefficients with respect to the basis of
standard unit vectors of `F` are given by the entries of `c`.
"""
FreeModElem_dec(c::SRow{T}, parent::FreeMod_dec{T}) where T = FreeModElem_dec{T}(c, parent)

@doc Markdown.doc"""
    FreeModElem_dec(c::Vector{T}, parent::FreeMod_dec{T}) where T

Return the element of `F` whose coefficients with respect to the basis of
standard unit vectors of `F` are given by the entries of `c`.
"""
function FreeModElem_dec(c::Vector{T}, parent::FreeMod_dec{T}) where T
  @assert length(c) == rank(parent)
  sparse_coords = sparse_row(base_ring(parent), collect(1:rank(parent)), c)
  return FreeModElem_dec{T}(sparse_coords,parent)
end

#@doc Markdown.doc"""
#    (F::FreeMod_dec{T})(c::SRow{T}) where T
#
#Return the element of `F` whose coefficients with respect to the basis of
#standard unit vectors of `F` are given by the entries of `c`.
#"""
function (F::FreeMod_dec{T})(c::SRow{T}) where T
  return FreeModElem_dec(c, F)
end

#@doc Markdown.doc"""
#    (F::FreeMod_dec{T})(c::Vector{T}) where T
#
#Return the element of `F` whose coefficients with respect to the basis of
#standard unit vectors of `F` are given by the entries of `c`.
#"""
function (F::FreeMod_dec{T})(c::Vector{T}) where T 
  return FreeModElem_dec(c, F)
end

@doc Markdown.doc"""
    (F::FreeMod_dec)()

Return the zero element of `F`.
"""
function (F::FreeMod_dec)()
  return FreeModElem_dec(sparse_row(base_ring(F)), F)
end

@doc Markdown.doc"""
    FreeModElem(coords::SRow{T}, parent::FreeMod_dec{T}) where T <: CRingElem_dec

Return the element of `F` whose coefficients with respect to the basis of 
standard unit vectors of `F` are given by the entries of `c`.
"""
function FreeModElem(coords::SRow{T}, parent::FreeMod_dec{T}) where T <: CRingElem_dec
  return FreeModElem_dec{T}(coords, parent)
end

@doc Markdown.doc"""
    FreeModElem_dec(v::FreeModElem{T}, parent::FreeMod_dec{T}) where T <: CRingElem_dec

Lift `v` to the decorated module `parent`.
"""
function FreeModElem_dec(v::FreeModElem{T}, p::FreeMod_dec{T}) where T <: CRingElem_dec
  @assert forget_decoration(p) === parent(v)
  return FreeModElem_dec(coordinates(v), p)
end


elem_type(::Type{FreeMod_dec{T}}) where {T} = FreeModElem_dec{T}
parent_type(::Type{FreeModElem_dec{T}}) where {T} = FreeMod_dec{T}
elem_type(::FreeMod_dec{T}) where {T} = FreeModElem_dec{T}
parent_type(::FreeModElem_dec{T}) where {T} = FreeMod_dec{T}

@doc Markdown.doc"""
"""
function forget_decoration(v::FreeModElem_dec)
  return FreeModElem(coordinates(v),forget_decoration(parent(v)))
end


@doc Markdown.doc"""
    generator_symbols(F::FreeMod_dec)

Return the list of symbols of the standard unit vectors.
"""
function generator_symbols(F::FreeMod_dec)
  return generator_symbols(forget_decoration(F))
end
@enable_all_show_via_expressify FreeModElem_dec


@doc Markdown.doc"""
    degree_homogeneous_helper(u::FreeModElem_dec)

Compute the degree and homogeneity of `u` (simultaneously).
A tuple is returned: The first entry is the degree (or nothing if
the element has no degree); the second entry is true if `u` is 
homogeneous and false otherwise.
"""
function degree_homogeneous_helper(u::FreeModElem_dec)
  if iszero(u)
    return nothing, true
  end
  first = true
  homogeneous = true #only needed in filtered case
  F = parent(u)
  W = base_ring(F)
  ww = W.D[0]
  local w
  for (p,v) in coordinates(u)
    if !is_homogeneous(v)
      if is_graded(W)
        return nothing,false
      else
        homogeneous = false
      end
    end
    w = degree(v)+F.d[p]
    if first
      ww = w
      first = false
    elseif is_graded(W)
      if ww != w
        return nothing, false
      end
    else
      if ww != w
        homogeneous = false
      end
      if W.lt(ww, w) 
        ww = w
      end
    end
  end
  return ww, homogeneous
end

@doc Markdown.doc"""
    degree(a::FreeModElem_dec)

Return the degree of `a`. If `a` has no degree an error is thrown.
"""
function degree(a::FreeModElem_dec)
  d,_ = degree_homogeneous_helper(a)
  d === nothing ? error("elem has no degree") : return d
end

@doc Markdown.doc"""
    homogeneous_components(a::FreeModElem_dec)

Return the homogeneous components of `a` in a dictionary.
The keys are those group elements for which `a` has a component
having this element as its degree.
"""
function homogeneous_components(a::FreeModElem_dec)
  res = Dict{GrpAbFinGenElem, FreeModElem_dec}()
  F = parent(a)
  for (p,v) in coordinates(a)
    c = homogeneous_components(v)
    for (pp, vv) in c
      w = pp + F.d[p]
      if haskey(res, w)
        res[w] += vv*gen(F, p)
      else
        res[w] = vv*gen(F, p)
      end
    end
  end
  return res
end

@doc Markdown.doc"""
    homogeneous_component(a::FreeModElem_dec, g::GrpAbFinGenElem)

Return the homogeneous component of `a` which has degree `g`.
"""
function homogeneous_component(a::FreeModElem_dec, g::GrpAbFinGenElem)
  F = parent(a)
  x = zero(F)
  for (p,v) in coordinates(a)
    x += homogeneous_component(v, g-F.d[p])*gen(F, p)
  end
  return x
end

@doc Markdown.doc"""
    is_homogeneous(a::FreeModElem_dec)

Check if `a` is homogeneous.
"""
function is_homogeneous(a::FreeModElem_dec)
  return degree_homogeneous_helper(a)[2]
end

# Weight vector or function?
# Should we already grade ModuleGens?
# Should it be possible to construct ungraded SubQuo with graded elements? (I.e. should the constructors
# accept AbstractFreeMod and AbstractFreeModElem instead of FreeMod and FreeModElem?)
# proceed with FreeModHom_dec?



@doc Markdown.doc"""
    tensor_product(G::FreeMod_dec...; task::Symbol = :none)

Given decorated free modules $G_i$ compute the decorated tensor product 
$G_1\otimes \cdots \otimes G_n$.
If `task` is set to ":map", a map $\phi$ is returned that
maps tuples in $G_1 \times \cdots \times G_n$ to pure tensors
$g_1 \otimes \cdots \otimes g_n$. The map admits a preimage as well.
"""
function tensor_product(G::FreeMod_dec...; task::Symbol = :none)
  undecorated_tensor_product, tuple_to_pure = tensor_product(map(forget_decoration, G)...; task=:map)
  pure_to_tuple = inv(tuple_to_pure)
  d = [sum(map(degree, [FreeModElem_dec(elem,parent) for (elem,parent) in zip(pure_to_tuple(v),G)])) 
                                                    for v in gens(undecorated_tensor_product)]
  F = FreeMod_dec(undecorated_tensor_product, d)

  function pure(T::Tuple)
    return FreeModElem_dec(tuple_to_pure(map(forget_decoration, T)), F)
  end

  function inv_pure(e::FreeModElem_dec)
    a = pure_to_tuple(forget_decoration(e))
    return Tuple(FreeModElem_dec(elem,parent) for (elem,parent) in zip(a,G))
  end

  set_attribute!(F, :tensor_pure_function => pure, :tensor_generator_decompose_function => inv_pure)

  if task == :none
    return F
  end
  return F, MapFromFunc(pure, inv_pure, Hecke.TupleParent(Tuple([g[0] for g = G])), F)
end


###############################################################################
# FreeModuleHom_dec constructors
###############################################################################

FreeModuleHom_dec(F::FreeMod_dec{T}, G::ModuleFP_dec, a::Vector) where {T} = FreeModuleHom_dec{T}(F, G, a)

FreeModuleHom_dec(F::FreeMod_dec{T}, G::ModuleFP_dec, mat::MatElem{T}) where {T} = FreeModuleHom{T}(F, G, mat)

function forget_decoration_on_morphism(f::FreeModuleHom_dec)
  return f.f
end

function forget_decoration(f::FreeModuleHom_dec)
  F = forget_decoration(domain(f))
  G = forget_decoration(codomain(f))
  return hom(F, G, [forget_decoration(f(v)) for v in gens(domain(f))])
end

function matrix(a::FreeModuleHom_dec)
  return matrix(forget_decoration_on_morphism(a))
end

(h::FreeModuleHom_dec)(a::FreeModElem_dec) = image(h, a)

hom(F::FreeMod_dec{T}, G::ModuleFP_dec{T}, V::Vector{<:FreeModElem_dec}) where T = FreeModuleHom_dec(F, G, V) 
hom(F::FreeMod_dec{T}, G::ModuleFP_dec{T}, A::MatElem{T}) where T = FreeModuleHom_dec(F, G, A)


function hom(F::FreeMod_dec, G::FreeMod_dec)
  undecorated_hom, elem_to_hom = hom(forget_decoration(F), forget_decoration(G))
  d = [y-x for x in decoration(F) for y in decoration(G)]
  GH = FreeMod_dec(undecorated_hom, d)
  X = Hecke.MapParent(F, G, "homomorphisms")

  function im(v::FreeModElem_dec)
    return hom(F, G, [FreeModElem_dec(elem_to_hom(forget_decoration(v))(forget_decoration(u)),G) for u in gens(F)])
  end

  function pre(f::FreeModuleHom_dec)
    undecorated_v = inv(elem_to_hom)(forget_decoration(f))
    return FreeModElem_dec(undecorated_v, GH)
  end

  to_hom_map = Hecke.MapFromFunc(im, pre, GH, X)
  set_attribute!(GH, :show => Hecke.show_hom, :hom => (F, G), :module_to_hom_map => to_hom_map)
  return GH, to_hom_map
end
