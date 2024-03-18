 function map(FR::FreeResolution, i::Int)
  return map(FR.C, i)
end

function free_show(io::IO, C::ComplexOfMorphisms)
  name_mod = String[]
  rank_mod = Int[]

  rng = range(C)
  rng = first(rng):-1:0
  arr = ("<--", "--")

  R = Nemo.base_ring(C[first(rng)])
  R_name = AbstractAlgebra.get_name(R)
  if isnothing(R_name)
    R_name = "$R"
  end
 
  for i=reverse(rng)
    M = C[i]
    M_name = AbstractAlgebra.get_name(M)
    if isnothing(M_name)
      M_name = "$R_name^$(rank(M))"
    end
    push!(name_mod, M_name)
    push!(rank_mod, rank(M))
  end

  io = IOContext(io, :compact => true)
  N = get_attribute(C, :free_res)
  if N !== nothing
    print(io, "Free resolution")
    print(io, " of ", N)
  end
  print(io, "\n")

  pos = 0
  pos_mod = Int[]
  
  for i=1:length(name_mod)
    print(io, name_mod[i])
    push!(pos_mod, pos)
    pos += length(name_mod[i])
    if i < length(name_mod)
      print(io, " ", arr[1], arr[2], " ")
      pos += length(arr[1]) + length(arr[2]) + 2
    end
  end

  print(io, "\n")
  len = 0
  for i=1:length(name_mod)
    if i>1
      print(io, " "^(pos_mod[i] - pos_mod[i-1]-len))
    end
    print(io, reverse(rng)[i])
    len = length("$(reverse(rng)[i])")
  end
#  print(io, "\n")
end


@doc raw"""
    free_resolution(F::FreeMod)

Return a free resolution of `F`. The `length` and `algorithm`
keywords are here only for compatibility reasons with the other `free_resolution`
methods and have no effect on the computation.

# Examples
"""
function free_resolution(F::FreeMod; length::Int=0, algorithm::Symbol=:fres)
  res = presentation(F)
  set_attribute!(res, :show => free_show, :free_res => F)
  return FreeResolution(res)
end

@doc raw"""
    is_complete(FR::FreeResolution)

Return `true` if the free resolution `fr` is complete, otherwise return `false`.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"])
(Multivariate polynomial ring in 3 variables over QQ, QQMPolyRingElem[x, y, z])

julia> A = R[x; y]
[x]
[y]

julia> B = R[x^2; x*y; y^2; z^4]
[x^2]
[x*y]
[y^2]
[z^4]

julia> M = SubquoModule(A, B)
Subquotient of Submodule with 2 generators
1 -> x*e[1]
2 -> y*e[1]
by Submodule with 4 generators
1 -> x^2*e[1]
2 -> x*y*e[1]
3 -> y^2*e[1]
4 -> z^4*e[1]

julia> fr = free_resolution(M, length=1)
Free resolution of M
R^2 <---- R^6
0         1

julia> is_complete(fr)
false

julia> fr = free_resolution(M)
Free resolution of M
R^2 <---- R^6 <---- R^6 <---- R^2 <---- 0
0         1         2         3         4

julia> is_complete(fr)
true

```
"""
is_complete(FR::FreeResolution) = FR.C.complete

function chain_range(FR::FreeResolution)
  return Hecke.range(FR.C)
end

function map_range(FR::FreeResolution)
  return Hecke.map_range(FR.C)
end

function chain_range(C::ComplexOfMorphisms)
  return Hecke.range(C)
end

function map_range(C::ComplexOfMorphisms)
  return Hecke.map_range(C)
end


#= Fill functions (and helpers) for Hecke ComplexOfMorphismses in terms of free resolutions =#
function _get_last_map_key(cc::Hecke.ComplexOfMorphisms)
  return last(Hecke.map_range(cc))
end

function _extend_free_resolution(cc::Hecke.ComplexOfMorphisms, idx::Int)
# assuming a free res is a chain_complex, then it will be
# M_1 -> M_0 -> S -> 0
#the range is 1:-1:-2 or so
#thus
# - extending right is trivial - and doing zero only
# - extending lift is repeated pushfirst
# - the idx is only used to see how many maps are missing

  algorithm = get_attribute(cc, :algorithm)
  if algorithm === nothing
    algorithm = :fres
    set_attribute!(cc, :algorithm, :fres)
  end
  r = Hecke.map_range(cc)
  if idx < last(r)
    error("extending past the final zero not supported")
  end
  len_missing = idx - first(r)
  @assert len_missing > 0
  if cc.complete == true
    return map(cc, first(r))
  end

  kernel_entry          = image(cc.maps[1])[1]
  br                    = base_ring(kernel_entry)
  singular_free_module  = singular_module(ambient_free_module(kernel_entry))
  singular_kernel_entry = Singular.Module(base_ring(singular_free_module),
                              [singular_free_module(repres(g)) for g in gens(kernel_entry)]...)
  singular_kernel_entry.isGB = true

  len = len_missing + 1
  if algorithm == :fres
    res = Singular.fres(singular_kernel_entry, len, "complete")
  elseif algorithm == :lres
    error("LaScala's method is not yet available in Oscar.")
  elseif algorithm == :mres
    res = Singular.mres(singular_kernel_entry, len)
  elseif algorithm == :nres
    res = Singular.nres(singular_kernel_entry, len)
  else
    error("Unsupported algorithm $algorithm")
  end

  dom = domain(cc.maps[1])
  j   = 2

  #= get correct length of Singular sresolution =#
  slen = iszero(res[Singular.length(res)+1]) ? Singular.length(res) : Singular.length(res)+1
  #= adjust length for extension length in Oscar =#
  slen = slen > len ? len : slen

  br_name = AbstractAlgebra.get_name(base_ring(kernel_entry))
  if br_name === nothing
    br_name = "R"
  end


  while j <= slen
    rk = Singular.ngens(res[j])
    if is_graded(dom)
      codom = dom
      SM    = SubModuleOfFreeModule(codom, res[j])
      #generator_matrix(SM)
      #map = graded_map(codom, SM.matrix) # going via matrices does a lot of unnecessary allocation and copying!
      map = graded_map(codom, gens(SM); check=false)
      dom = domain(map)
      AbstractAlgebra.set_name!(dom, "$br_name^$rk")
    else
      codom = dom
      dom   = free_module(br, Singular.ngens(res[j]))
      SM    = SubModuleOfFreeModule(codom, res[j])
      AbstractAlgebra.set_name!(dom, "$br_name^$rk")
      #generator_matrix(SM)
      map = hom(dom, codom, gens(SM); check=false)
    end
    pushfirst!(cc, map) 
    j += 1
  end
  # Finalize maps.
  if slen < len
    Z = FreeMod(br, 0)
    AbstractAlgebra.set_name!(Z, "0")
    pushfirst!(cc, hom(Z, domain(cc.maps[1]), Vector{elem_type(domain(cc.maps[1]))}(); check=false))
    cc.complete = true
  end
  set_attribute!(cc, :show => free_show)
  maxidx = min(idx, first(Hecke.map_range(cc)))
  return map(cc, maxidx)
end

@doc raw"""
    free_resolution(M::SubquoModule{<:MPolyRingElem}; 
        ordering::ModuleOrdering = default_ordering(M),
        length::Int=0, algorithm::Symbol=:fres
      )

Return a free resolution of `M`.

If `length != 0`, the free resolution is only computed up to the `length`-th free module.
At the moment, options for `algorithm` are `:fres`, `:mres` and `:nres`. With `:mres` or `:nres`,
minimal free resolutions are returned.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"])
(Multivariate polynomial ring in 3 variables over QQ, QQMPolyRingElem[x, y, z])

julia> A = R[x; y]
[x]
[y]

julia> B = R[x^2; x*y; y^2; z^4]
[x^2]
[x*y]
[y^2]
[z^4]

julia> M = SubquoModule(A, B)
Subquotient of Submodule with 2 generators
1 -> x*e[1]
2 -> y*e[1]
by Submodule with 4 generators
1 -> x^2*e[1]
2 -> x*y*e[1]
3 -> y^2*e[1]
4 -> z^4*e[1]

julia> fr = free_resolution(M, length=1)
Free resolution of M
R^2 <---- R^6
0         1

julia> is_complete(fr)
false

julia> fr[4]
Free module of rank 0 over Multivariate polynomial ring in 3 variables over QQ

julia> fr
Free resolution of M
R^2 <---- R^6 <---- R^6 <---- R^2 <---- 0
0         1         2         3         4

julia> is_complete(fr)
true

julia> fr = free_resolution(M, algorithm=:fres)
Free resolution of M
R^2 <---- R^6 <---- R^6 <---- R^2 <---- 0
0         1         2         3         4
```

**Note:** Over rings other than polynomial rings, the method will default to a lazy, 
iterative kernel computation.
"""
function free_resolution(M::SubquoModule{<:MPolyRingElem}; 
                         ordering::ModuleOrdering = default_ordering(M),
                         length::Int=0, algorithm::Symbol=:fres)

  coefficient_ring(base_ring(M)) isa AbstractAlgebra.Field ||
      error("Must be defined over a field.")

  cc_complete = false

  #= Start with presentation =#
  pm = presentation(M)
  maps = [pm.maps[j] for j in 2:3]

  br = base_ring(M)
  kernel_entry          = image(pm.maps[1])[1]

  if ngens(kernel_entry) == 0
    cc = Hecke.ComplexOfMorphisms(Oscar.ModuleFP, pushfirst!(maps, pm.maps[1]), check = false, seed = -2)
    cc.fill     = _extend_free_resolution
    cc.complete = true
    return FreeResolution(cc)
  end

  singular_free_module  = singular_module(ambient_free_module(kernel_entry))
  singular_kernel_entry = Singular.Module(base_ring(singular_free_module),
                              [singular_free_module(repres(g)) for g in gens(kernel_entry)]...)

  #= This is the single computational hard part of this function =#
  if algorithm == :fres
    gbpres = Singular.std(singular_kernel_entry)
    res = Singular.fres(gbpres, length, "complete")
  elseif algorithm == :lres
    error("LaScala's method is not yet available in Oscar.")
    gbpres = singular_kernel_entry # or as appropriate, taking into account base changes
  elseif algorithm == :mres
    gbpres = singular_kernel_entry
    res = Singular.mres(gbpres, length)
  elseif algorithm == :nres
    gbpres = singular_kernel_entry
    res = Singular.nres(gbpres, length)
  else
    error("Unsupported algorithm $algorithm")
  end

  slen = iszero(res[Singular.length(res)+1]) ? Singular.length(res) : Singular.length(res)+1
  if length == 0 || slen < length
    cc_complete = true
  end

  if length != 0
    slen =  slen > length ? length : slen
  end

  br_name = AbstractAlgebra.get_name(base_ring(M))
  if br_name === nothing
    br_name = "R"
  end

  #= Add maps from free resolution computation, start with second entry
   = due to inclusion of presentation(M) at the beginning. =#
  j   = 1
  while j <= slen
    if is_graded(M)
      codom = domain(maps[1])
      rk    = Singular.ngens(res[j])
      SM    = SubModuleOfFreeModule(codom, res[j])
      #generator_matrix(SM)
      #ff = graded_map(codom, SM.matrix)
      ff = graded_map(codom, gens(SM); check=false)
      dom = domain(ff)
      AbstractAlgebra.set_name!(dom, "$br_name^$rk")
      insert!(maps, 1, ff)
      j += 1
    else
      codom = domain(maps[1])
      rk    = Singular.ngens(res[j])
      dom   = free_module(br, rk)
      SM    = SubModuleOfFreeModule(codom, res[j])
      #generator_matrix(SM)
      AbstractAlgebra.set_name!(dom, "$br_name^$rk")
      insert!(maps, 1, hom(dom, codom, gens(SM); check=false))
      j += 1
    end
  end
  if cc_complete == true
    # Finalize maps.
    if is_graded(domain(maps[1]))
      Z = graded_free_module(br, 0)
    else
      Z = FreeMod(br, 0)
    end
    AbstractAlgebra.set_name!(Z, "0")
    insert!(maps, 1, hom(Z, domain(maps[1]), Vector{elem_type(domain(maps[1]))}(); check=false))
  end

  cc = Hecke.ComplexOfMorphisms(Oscar.ModuleFP, maps, check = false, seed = -2)
  cc.fill     = _extend_free_resolution
  cc.complete = cc_complete
  set_attribute!(cc, :show => free_show, :free_res => M)
  set_attribute!(cc, :algorithm, algorithm)

  return FreeResolution(cc)
end

function free_resolution(M::SubquoModule{T}) where {T<:RingElem}
  # This generic code computes a free resolution in a lazy way.
  # We start out with a presentation of M and implement 
  # an iterative fill function to compute every higher term 
  # on request.
  R = base_ring(M)
  p = presentation(M)
  p.fill = function(C::Hecke.ComplexOfMorphisms, k::Int)
    # TODO: Use official getter and setter methods instead 
    # of messing manually with the internals of the complex.
    for i in first(chain_range(C)):k-1
      N = domain(map(C, i))

      if iszero(N) # Fill up with zero maps
        C.complete = true
        phi = hom(N, N, elem_type(N)[]; check=false)
        pushfirst!(C.maps, phi)
        continue
      end

      K, inc = kernel(map(C, i))
      nz = findall(x->!iszero(x), gens(K))
      F = FreeMod(R, length(nz))
      iszero(length(nz)) && AbstractAlgebra.set_name!(F, "0")
      phi = hom(F, C[i], iszero(length(nz)) ? elem_type(C[i])[] : inc.(gens(K)[nz]); check=false)
      pushfirst!(C.maps, phi)
    end
    return first(C.maps)
  end
  return p
end


@doc raw"""
    free_resolution_via_kernels(M::SubquoModule, limit::Int = -1)

Return a free resolution of `M`.

If `limit != -1`, the free resolution
is only computed up to the `limit`-th free module.

# Examples
"""
function free_resolution_via_kernels(M::SubquoModule, limit::Int = -1)
  p = presentation(M)
  mp = [map(p, j) for j in Hecke.map_range(p)]  
  while true
    k, mk = kernel(mp[1])
    nz = findall(x->!iszero(x), gens(k))
    if length(nz) == 0 
      if is_graded(domain(mp[1]))
        h = graded_map(domain(mp[1]), Vector{elem_type(domain(mp[1]))}(); check=false)
      else
        Z = FreeMod(base_ring(M), 0)
        AbstractAlgebra.set_name!(Z, "0")
        h = hom(Z, domain(mp[1]), Vector{elem_type(domain(mp[1]))}(); check=false)
      end
      insert!(mp, 1, h)
      break
    elseif limit != -1 && length(mp) > limit
      break
    end
    if is_graded(codomain(mk))
      g = graded_map(codomain(mk), collect(k.sub.gens)[nz]; check=false)
    else
      F = FreeMod(base_ring(M), length(nz))
      g = hom(F, codomain(mk), collect(k.sub.gens)[nz]; check=false)
    end
    insert!(mp, 1, g)
  end
  C = Hecke.ComplexOfMorphisms(ModuleFP, mp, check = false, seed = -2)
  #set_attribute!(C, :show => free_show, :free_res => M) # doesn't work
  return FreeResolution(C)
end

@doc raw"""
    free_resolution(I::MPolyIdeal; length::Int=0, algorithm::Symbol=:fres)

Compute a free resolution of `I`.

If `length != 0`, the free resolution is only computed up to the `length`-th free module.
At the moment, options for `algorithm` are `:fres`, `:mres` and `:nres`. With `:mres` or `:nres`,
minimal free resolutions are returned.

# Examples
"""
function free_resolution(I::MPolyIdeal;
                         length::Int=0, algorithm::Symbol=:fres)
  S = ideal_as_module(I)
  n = AbstractAlgebra.get_name(I)
  if n !== nothing
    AbstractAlgebra.set_name!(S, n)
  end
  return free_resolution(S, length = length, algorithm = algorithm)
end

@doc raw"""
    free_resolution(Q::MPolyQuoRing; length::Int=0, algorithm::Symbol=:fres)

Compute a free resolution of `Q`.

If `length != 0`, the free resolution is only computed up to the `length`-th free module.
At the moment, options for `algorithm` are `:fres`, `:mres` and `:nres`. With `:mres` or `:nres`,
minimal free resolutions are returned.

# Examples
"""
function free_resolution(Q::MPolyQuoRing;
                         length::Int=0, algorithm::Symbol=:fres)
  q = quotient_ring_as_module(Q)
  n = AbstractAlgebra.get_name(Q)
  if n !== nothing
    AbstractAlgebra.set_name!(q, n)
  end
  return free_resolution(q, length = length, algorithm = algorithm)
end
