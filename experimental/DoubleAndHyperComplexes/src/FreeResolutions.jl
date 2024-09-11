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
    R_name = "($R)"
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

  while j <= slen
    rk = Singular.ngens(res[j])
    if is_graded(dom)
      codom = dom
      SM    = SubModuleOfFreeModule(codom, res[j])
      #generator_matrix(SM)
      #map = graded_map(codom, SM.matrix) # going via matrices does a lot of unnecessary allocation and copying!
      map = graded_map(codom, gens(SM); check=false)
      dom = domain(map)
    else
      codom = dom
      dom   = free_module(br, Singular.ngens(res[j]))
      SM    = SubModuleOfFreeModule(codom, res[j])
      #generator_matrix(SM)
      map = hom(dom, codom, gens(SM); check=false)
    end
    pushfirst!(cc, map) 
    j += 1
  end
  # Finalize maps.
  if slen < len
    Z = FreeMod(br, 0)
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
        length::Int = 0, algorithm::Symbol = :fres
      )

Return a free resolution of `M`.

If `length != 0`, the free resolution is only computed up to the `length`-th free module.
Current options for `algorithm` are `:fres`, `:nres`, and `:mres`.

!!! note
    The function first computes a presentation of `M`. It then successively computes
    higher syzygy modules.

!!! note
    - If `algorithm == mres`, and `M` is positively graded, a minimal free resolution is returned.
    - If `algorithm == nres`, and `M` is positively graded, the function proceeds as above except 
      that it starts by computing a presentation which is not neccessarily minimal.
    In both cases, if `M` is not (positively) graded, the function still aims at returning an ''improved'' resolution.

!!! note
    If `algorithm == fres`, the function relies on an enhanced version of Schreyer's algorithm 
    [EMSS16](@cite). Typically, this is more efficient than the approaches above, but the 
    resulting resolution is far from being minimal.

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
Free module of rank 0 over R

julia> fr
Free resolution of M
R^2 <---- R^6 <---- R^6 <---- R^2 <---- 0
0         1         2         3         4

julia> is_complete(fr)
true
```

```
julia> R, (w, x, y, z) = graded_polynomial_ring(QQ, ["w", "x", "y", "z"]);

julia> Z = R(0)
0

julia> O = R(1)
1

julia> B = [Z Z Z O; w*y w*z-x*y x*z-y^2 Z];

julia> A = transpose(matrix(B));

julia> M = graded_cokernel(A)
Graded subquotient of submodule of R^2 generated by
1 -> e[1]
2 -> e[2]
by submodule of R^2 generated by
1 -> w*y*e[2]
2 -> (w*z - x*y)*e[2]
3 -> (x*z - y^2)*e[2]
4 -> e[1]

julia> FM1 = free_resolution(M)
Free resolution of M
R^2 <---- R^7 <---- R^8 <---- R^3 <---- 0
0         1         2         3         4

julia> betti_table(FM1)
       0  1  2  3
-----------------
-1   : -  1  -  -
0    : 2  -  -  -
1    : -  3  3  1
2    : -  3  5  2
-----------------
total: 2  7  8  3


julia> matrix(map(FM1, 1))
[1            0]
[0   -x*z + y^2]
[0   -w*z + x*y]
[0          w*y]
[0        x^2*z]
[0        w*x*z]
[0        w^2*z]

julia> FM2 = free_resolution(M, algorithm = :nres)
Free resolution of M
R^2 <---- R^4 <---- R^4 <---- R^2 <---- 0
0         1         2         3         4

julia> betti_table(FM2)
       0  1  2  3
-----------------
-1   : -  1  -  -
0    : 2  -  -  -
1    : -  3  -  -
2    : -  -  4  2
-----------------
total: 2  4  4  2


julia> matrix(map(FM2, 1))
[1            0]
[0   -x*z + y^2]
[0   -w*z + x*y]
[0          w*y]

julia> FM3 = free_resolution(M, algorithm = :mres)
Free resolution of M
R^1 <---- R^3 <---- R^4 <---- R^2 <---- 0
0         1         2         3         4

julia> betti_table(FM3)
       0  1  2  3
-----------------
0    : 1  -  -  -
1    : -  3  -  -
2    : -  -  4  2
-----------------
total: 1  3  4  2


julia> matrix(map(FM3, 1))
[-x*z + y^2]
[-w*z + x*y]
[       w*y]

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
  pm = algorithm == :mres ? _presentation_minimal(M, minimal_kernel=false) : presentation(M)
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
      insert!(maps, 1, ff)
      j += 1
    else
      codom = domain(maps[1])
      rk    = Singular.ngens(res[j])
      dom   = free_module(br, rk)
      SM    = SubModuleOfFreeModule(codom, res[j])
      #generator_matrix(SM)
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

###############################################################################
# Graded free resolutions
###############################################################################

function is_graded(resolution::FreeResolution{T}) where T
  C = resolution.C
 return all(is_graded(C[i]) for i in reverse(Hecke.range(C))) && all(is_graded(map(C, i)) for i in reverse(Hecke.map_range(C)))
end

###############################################################################
# Betti table
###############################################################################

@doc raw"""
    betti_table(F::FreeResolution)

Given a $\mathbb Z$-graded free resolution `F`, return the graded Betti numbers 
of `F` in form of a Betti table.

Alternatively, use `betti`.

# Examples
```julia
julia> R, (w, x, y, z) = graded_polynomial_ring(QQ, ["w", "x", "y", "z"]);

julia> I = ideal(R, [x*z, y*z, x*w^2, y*w^2])
ideal(x*z, y*z, w^2*x, w^2*y)

julia> A, _= quo(R, I)
(Quotient of multivariate polynomial ring by ideal with 4 generators, Map from
R to A defined by a julia-function with inverse)

julia> FA  = free_resolution(A)
Free resolution of A
R^1 <---- R^4 <---- R^4 <---- R^1 <---- 0
0         1         2         3         4

julia> betti_table(FA)
       0  1  2  3
------------------
0    : 1  -  -  -
1    : -  2  1  -
2    : -  2  3  1
------------------
total: 1  4  4  1

julia> R, (x, y) = graded_polynomial_ring(QQ, ["x", "y"]);

julia> I = ideal(R, [x, y, x+y]);

julia> M = quotient_ring_as_module(I);

julia> FM = free_resolution(M, algorithm = :nres);

julia> betti_table(FM)
       0  1  2
---------------
-1   : -  -  1
0    : 1  3  1
---------------
total: 1  3  2
```
"""
function betti_table(F::FreeResolution; project::Union{FinGenAbGroupElem, Nothing} = nothing, reverse_direction::Bool=false)
  generator_count = Dict{Tuple{Int, Any}, Int}()
  C = F.C
  rng = Hecke.map_range(C)
  n = first(rng)
  for i in 0:n
    module_degrees = F[i].d
    module_degrees === nothing && error("One of the modules in the graded free resolution is not graded.")
    for degree in module_degrees
      idx = (i, degree)
      generator_count[idx] = get(generator_count, idx, 0) + 1
    end
  end
  return BettiTable(generator_count, project = project, reverse_direction = reverse_direction)
end

function betti(b::FreeResolution; project::Union{FinGenAbGroupElem, Nothing} = nothing, reverse_direction::Bool = false)
  return betti_table(b; project, reverse_direction)
end

########################################################################
# Additional functionality
########################################################################
function hom_without_reversing_direction(F::FreeResolution, M::ModuleFP)
  return hom_without_reversing_direction(F.C, M)
end

### This is overwritten for compatibility with the ComplexOfMorphisms
function homology(C::FreeResolution)
  return homology(C.C)
end

########################################################################
# Flattenings of free resolutions
########################################################################
function free_resolution(
    M::SubquoModule{T}
  ) where {T <: FlattableRingElemType}
  flat = flatten(base_ring(M))
  M_flat, iso_M, iso_M_inv = flat(M)
  comp = free_resolution(M_flat) # assuming that this is potentially cached
  if !haskey(flat_counterparts(flat), comp)
    res_obj = ModuleFP[]
    isos = ModuleFPHom[]
    push!(res_obj, M)
    push!(isos, iso_M_inv)
    res_maps = Map[]
    for i in 0:first(chain_range(comp))
      F_flat = comp[i]
      @assert F_flat === domain(map(comp, i))
      @assert domain(last(isos)) === codomain(map(comp, i))
      F, iso_F_inv = _change_base_ring_and_preserve_gradings(inverse(flat), F_flat)
      iso_F = hom(F, F_flat, gens(F_flat), flat; check=false)
      push!(res_maps, hom(F, last(res_obj), last(isos).(map(comp, i).(iso_F.(gens(F)))); check=false))
      push!(res_obj, F)
      push!(isos, iso_F_inv)
    end
    comp_up = ComplexOfMorphisms(ModuleFP, reverse(res_maps), typ=:chain, seed=-1, check=false)
    comp_up.complete = true
    result = FreeResolution(comp_up)
    flat_counterparts(flat)[comp] = result
  end
  return flat_counterparts(flat)[comp]::FreeResolution
end

########################################################################
# Functionality for graded modules
########################################################################
@doc raw"""
    minimal_betti_table(F::FreeResolution{T}; check::Bool=true) where {T<:ModuleFP}

Given a graded free resolution `F` over a standard $\mathbb Z$-graded 
multivariate polynomial ring with coefficients in a field, return the
Betti table of the minimal free resolution arising from `F`.

!!! note
    The algorithm proceeds without actually minimizing the resolution.

# Examples
```julia
julia> R, (w, x, y, z) = graded_polynomial_ring(QQ, ["w", "x", "y", "z"]);

julia> I = ideal(R, [w^2-x*z, w*x-y*z, x^2-w*y, x*y-z^2, y^2-w*z]);

julia> A, _ = quo(R, I)
(Quotient of multivariate polynomial ring by ideal with 5 generators, Map from
R to A defined by a julia-function with inverse)

julia> FA = free_resolution(A)
Free resolution of A
R^1 <---- R^5 <---- R^6 <---- R^2 <---- 0
0         1         2         3         4

julia> betti_table(FA)
       0  1  2  3  
------------------
0    : 1  -  -  -  
1    : -  5  5  1  
2    : -  -  1  1  
------------------
total: 1  5  6  2  


julia> minimal_betti_table(FA)
       0  1  2  3  
------------------
0    : 1  -  -  -  
1    : -  5  5  -  
2    : -  -  -  1  
------------------
total: 1  5  5  1  
```
"""
function minimal_betti_table(res::FreeResolution{T}; check::Bool=true) where {T<:ModuleFP}
  @assert is_standard_graded(base_ring(res)) "resolution must be defined over a standard graded ring"
  @assert is_graded(res) "resolution must be graded"
  C = complex(res)
  @assert is_complete(res) "resolution must be complete"
  rng = range(C)
  # The following needs the resolution to be complete to be true
  res_length = first(rng)-1
  offsets = Dict{FinGenAbGroupElem, Int}()
  betti_hash_table = Dict{Tuple{Int, Any}, Int}()
  for i in 1:res_length+1
    phi = map(C, i)
    F = domain(phi)
    G = codomain(phi)
    dom_degs = unique!([degree(g; check) for g in gens(F)])
    cod_degs = unique!([degree(g; check) for g in gens(G)])
    for d in cod_degs
      d::FinGenAbGroupElem
      if d in dom_degs
        _, _, sub_mat = _constant_sub_matrix(phi, d; check)
        r = rank(sub_mat)
        c = ncols(sub_mat) - r - get(offsets, d, 0)
        !iszero(c) && (betti_hash_table[(i-1, d)] = c)
        offsets[d] = r
      else
        c = length(_indices_of_generators_of_degree(G, d; check)) - get(offsets, d, 0)
        !iszero(c) && (betti_hash_table[(i-1, d)] = c)
      end
    end
  end
  return BettiTable(betti_hash_table)
end

function hash_table(B::BettiTable) 
  return B.B
end

# TODO: Where is this called??? Adjust the use of `check` there!
function generators_of_degree(
    C::FreeResolution{T},
    i::Int,
    d::FinGenAbGroupElem;
    check::Bool=true
  ) where {T<:ModuleFP}
  F = C[i]
  return [g for g in gens(F) if degree(g) == d]
end

function _indices_of_generators_of_degree(
    F::FreeMod{T}, d::FinGenAbGroupElem; check::Bool=true
  ) where {T<:Union{MPolyDecRingElem{<:FieldElem}, 
                    MPolyQuoRingElem{<:MPolyDecRingElem{<:FieldElem}}}}
  return Int[i for (i, g) in enumerate(gens(F)) if degree(g; check) == d]
end

function _constant_sub_matrix(
    phi::FreeModuleHom{T, T},
    d::FinGenAbGroupElem;
    check::Bool=true
  ) where {RET<:Union{MPolyDecRingElem{<:FieldElem}, 
                      MPolyQuoRingElem{<:MPolyDecRingElem{<:FieldElem}}
                     }, T<:FreeMod{RET}}
  S = base_ring(domain(phi))::Union{MPolyDecRing, MPolyQuoRing{<:MPolyDecRingElem}}
  kk = coefficient_ring(S)::Field
  F = domain(phi)
  G = codomain(phi)
  ind_dom = _indices_of_generators_of_degree(F, d; check)
  ind_cod = _indices_of_generators_of_degree(G, d; check)
  m = length(ind_dom)
  n = length(ind_cod)
  result = zero_matrix(kk, m, n)
  img_gens = images_of_generators(phi)
  for (i, l) in enumerate(ind_dom)
    v = coordinates(img_gens[l])
    for (j, k) in enumerate(ind_cod)
      success, c = _has_index(v, k)
      !success && continue
      result[i, j] = _constant_coeff(c)
    end
  end
  return ind_dom, ind_cod, result
end

_constant_coeff(f::MPolyDecRingElem) = first(coefficients(f))
_constant_coeff(f::MPolyQuoRingElem) = first(coefficients(lift(f)))

# TODO: This will be provided soon from different sources.
function complex(F::FreeResolution) 
  return F.C
end

function base_ring(res::FreeResolution{T}) where {T<:ModuleFP}
  return base_ring(res[-1])
end

### Additional code to make the tests run again
#=
function map_range(comp::AbsHyperComplex)
  @assert isone(dim(comp))
  @assert has_upper_bound(comp, 1)
  @assert has_lower_bound(comp, 1)
  if direction(comp, 1) == :chain
    return upper_bound(comp, 1):-1:(lower_bound(comp, 1) + 1)
  else
    return lower_bound(comp, 1):(upper_bound(comp, 1) - 1)
  end
end
=#

function chain_range(comp::AbsHyperComplex)
  @assert isone(dim(comp))
  @assert has_upper_bound(comp, 1)
  @assert has_lower_bound(comp, 1)
  if direction(comp, 1) == :chain
    return upper_bound(comp, 1):-1:(lower_bound(comp, 1))
  else
    return lower_bound(comp, 1):(upper_bound(comp, 1))
  end
end

