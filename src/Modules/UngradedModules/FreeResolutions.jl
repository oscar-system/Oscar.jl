function map(FR::FreeResolution, i::Int)
  return map(FR.C, i)
end

function Base.show(io::IO, FR::FreeResolution)
  C = FR.C
  name_mod = String[]
  rank_mod = Int[]

  rng = range(C)
  ub = findfirst(i -> is_zero(C[i]), 2:first(map_range(C)))
  if ub == nothing
    rng = first(map_range(C)):-1:0
  else
    rng = ub+1:-1:0
  end

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

  io = terse(io)
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
  set_attribute!(res, :show => Hecke.pres_show, :free_res => F)
  return FreeResolution(res)
end

@doc raw"""
    is_complete(F::FreeResolution)

Return `true` if `F` is complete, `false` otherwise.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
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
Subquotient of submodule with 2 generators
  1: x*e[1]
  2: y*e[1]
by submodule with 4 generators
  1: x^2*e[1]
  2: x*y*e[1]
  3: y^2*e[1]
  4: z^4*e[1]

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

function _extend_free_resolution_to_the_left_by_zeros(cc::Hecke.ComplexOfMorphisms, idx::Int)
  r = Hecke.map_range(cc)
  if idx > first(r)
    for j = first(r):idx-1
      cod =  domain(first(cc.maps))
      phi = id_hom(cod)
      pushfirst!(cc.maps, phi)
    end
  end
end

function _extend_free_resolution(cc::Hecke.ComplexOfMorphisms, idx::Int)
# assuming a free res is a chain_complex, then it will be
# M_1 -> M_0 -> S -> 0
#the range is 1:-1:-2 or so
#thus
# - extending right is trivial - and doing zero only
# - extending left is repeated pushfirst
# - the idx is only used to see how many maps are missing

  algorithm = get_attribute(cc, :algorithm)
  if algorithm === nothing && base_ring(cc[-1]) isa MPolyRing
    algorithm = :fres
    set_attribute!(cc, :algorithm, :fres)
  end
  if algorithm === nothing && base_ring(cc[-1]) isa MPolyQuoRing
    algorithm = :sres
    set_attribute!(cc, :algorithm, :sres)
  end

  r = Hecke.map_range(cc)

  # extending to the right by zeros

  if idx < last(r)
    for j = last(r)-1:-1:idx
      cod =  codomain(last(cc.maps))
      phi = id_hom(cod)
      push!(cc.maps, phi)
      cc.seed = idx-1
    end
    return last(cc.maps)
  end

  len_missing = idx - first(r)
  @assert len_missing > 0
  if cc.complete == true
    _extend_free_resolution_to_the_left_by_zeros(cc, idx)
    return map(cc, first(r))
  end

  # The current version is a workaround; see #4998. Until we can provide a 
  # reasonable fix using a continuation of Schreyer's ordering, we 
  # provide a patch with new initiation of kernel computations.
  if algorithm == :fres || algorithm == :sres
    # In case of `:fres` (or `:sres` for quotient rings), 
    # we start the computation of the resolution 
    # from scratch. Continuation is not possible at the moment, because 
    # the information on the Schreyer orderings is lost in Singular internals. 
    # Continuation via kernel computations is, in most examples, more expensive 
    # than starting again from the beginning.
    return _extend_free_resolution_via_fres_or_sres(cc, idx, algorithm)
  end

  phi = first(cc.maps)
  K, inc_K = kernel(phi)
  R = base_ring(K)
  F = ambient_free_module(K)
  SF = singular_module(F)
  SK = singular_generators(algorithm == :fres || algorithm == :sres ? groebner_basis(K) : K.sub.gens) # We need to start with a gb for `fres`.

  if is_zero(K)
      cc.complete=true
      cod = domain(first(cc.maps))
      phi = is_graded(cod) ? graded_map(cod, ambient_representatives_generators(K); check=false) : hom(free_module(R, ngens(K)), cod, ambient_representatives_generators(K))
      pushfirst!(cc.maps, phi)
      _extend_free_resolution_to_the_left_by_zeros(cc, idx)
      return first(cc.maps)
  end

  if is_one(len_missing)
    cod = domain(first(cc.maps))
    phi = is_graded(cod) ? graded_map(cod, ambient_representatives_generators(K); check=false) : hom(free_module(R, ngens(K)), cod, ambient_representatives_generators(K))
    pushfirst!(cc.maps, phi)
    return first(cc.maps)
  end

  if algorithm == :fres
    error("this case should have been caught above")
    res = Singular.fres(SK, len_missing, "complete")
  elseif algorithm == :lres
    error("LaScala's method is not yet available in Oscar.")
  elseif algorithm == :mres
    res = Singular.mres(SK, len_missing)
  elseif algorithm == :nres
    res = Singular.nres(SK, len_missing)
  elseif algorithm == :sres
    error("this case should have been caught above")
    res = Singular.sres(SK, len_missing)
  else
    error("Unsupported algorithm $algorithm")
  end

  slen = Singular.length(res)
  slen =  slen > len_missing ? len_missing : slen #TODO check why the numbers may differ (Singular.jl)
  for j in 1:slen
    cod = domain(first(cc.maps))
    I = SubModuleOfFreeModule(cod, res[j]) # convert to an Oscar module
    phi = is_graded(cod) ? graded_map(cod, gens(I); check=false) : hom(free_module(R, ngens(I)), cod, gens(I))
    pushfirst!(cc.maps, phi)
    if is_zero(I)
      cc.complete = true
      _extend_free_resolution_to_the_left_by_zeros(cc, idx)
      return first(cc.maps)
    end
  end
  if is_zero(res[slen+1]) #TODO check why this is defined (Singular.jl)
    cod = domain(first(cc.maps))
    pushfirst!(cc.maps, is_graded(cod) ? graded_map(cod, elem_type(cod)[]; check=false) : hom(free_module(R, 0), cod, elem_type(cod)[]))
    cc.complete = true
    _extend_free_resolution_to_the_left_by_zeros(cc, idx)
    return cc.maps[len_missing > Singular.length(res) ? 1 : 2]
  end
  return first(cc.maps)
end

# For an explanation see the call to this method above.
# We need this dummy method, because the parent objects and 
# morphisms present in the partial resolution `cc` might have 
# been given to the user already. So even though we compute a 
# new resolution from scratch, it needs to be spliced with the 
# partial one which was already there. 
function _extend_free_resolution_via_fres_or_sres(cc::Hecke.ComplexOfMorphisms, idx::Int, algorithm::Symbol)
  M = cc[-1]
  res = free_resolution(M; algorithm, length=idx)
  i0 = first(range(cc))
  phi = map(res, i0+1)
  @assert ngens(codomain(phi)) == ngens(cc[i0]) "Schreyer resolutions are not compatible"
  Fi0 = cc[i0]
  img_gens = elem_type(cc[i0])[sum(c*Fi0[j] for (j, c) in coordinates(v); init=zero(Fi0)) for v in images_of_generators(phi)]
  pushfirst!(cc.maps, hom(res[i0+1], Fi0, img_gens))
  for j in first(range(cc))+1:first(range(res.C))
    pushfirst!(cc.maps, map(res, j))
  end
  if is_zero(domain(first(cc.maps)))
    cc.complete=true
  end
  _extend_free_resolution_to_the_left_by_zeros(cc, idx)
  return first(cc.maps)
end


@doc raw"""
    free_resolution(M::SubquoModule{T}; 
        length::Int = 0,
        algorithm::Symbol = T <: MPolyRingElem ? :fres : :sres) where {T <: Union{MPolyRingElem, MPolyQuoRingElem}}

Return a free resolution of `M`.

If `length != 0`, then the free resolution is only computed up to the `length`-th free module.
Current options for `algorithm` are `:mres` and `:nres`, as well as `:fres` (for multivariate
polynomial rings) and `:sres` (for quotients of multivariate polynomial rings), see below.

!!! note
    The function first computes a presentation of `M`. It then successively computes
    higher syzygy modules. In the first example of the examples section below, the free
    resolution is initially computed up to length 1:
    ```julia-repl
    julia> fr = free_resolution(M, length = 1)
    Free resolution of M
    R^2 <---- R^6
    0         1
    ```
    This resolution is not yet complete
    ```julia-repl
    julia> is_complete(fr)
    false
    ```
    Continuing the session as follows, the resolution is extended up to length 4:
    ```julia-repl
    julia> fr[4]
    Free module of rank 0 over R

    julia> fr
    Free resolution of M
    R^2 <---- R^6 <---- R^6 <---- R^2 <---- 0
    0         1         2         3         4
    ``` 
    As we already see from the output, the extended resolution is complete:
    ```julia-repl
    julia> is_complete(fr)
    true
    ```

!!! note
    If `M` is a module over a quotient of a polynomial ring, then the `length` keyword must
    be set to a nonzero value (recall that free resolutions over quotient rings may be infinite).

!!! note
    - If `algorithm` is set to `:mres`, and `M` is positively graded, then a minimal free resolution is returned.
    - If `algorithm` is set to `:nres`, and `M` is positively graded, then the function proceeds as above except
      that it starts by computing a presentation which is not necessarily minimal.
    In both cases, if `M` is not (positively) graded, then the function still aims at returning an ''improved'' resolution.

!!! note
    If the default options `algorithm = :fres` (for multivariate polynomial rings) or `algorithm = :sres` (for quotients of multivariate polynomial rings) are used, 
    then the function relies on enhanced versions of Schreyer's algorithm [EMSS16](@cite). This is often, but not always, more efficient
    than the methods specified by `algorithm = :mres` or `algorithm = :nres`. The resulting Schreyer resolution, however, may be far from minimal.
    Here is an extract from an OSCAR session illustrating this:
    ```julia-repl
    julia> FM = free_resolution(M)
    Free resolution of M
    Pn^44 <---- Pn^296 <---- Pn^808 <---- Pn^1019 <---- Pn^618 <---- Pn^169 <---- Pn^14 <---- 0
    0           1            2            3             4            5            6           7

    julia> free_resolution(M, algorithm = :mres)
    Free resolution of M
    Pn^6 <---- Pn^10 <---- Pn^4 <---- 0
    0          1           2          3

    ```

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
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
Subquotient of submodule with 2 generators
  1: x*e[1]
  2: y*e[1]
by submodule with 4 generators
  1: x^2*e[1]
  2: x*y*e[1]
  3: y^2*e[1]
  4: z^4*e[1]

julia> fr = free_resolution(M, length = 1)
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

```jldoctest
julia> R, (w, x, y, z) = graded_polynomial_ring(QQ, [:w, :x, :y, :z]);

julia> Z = R(0)
0

julia> O = R(1)
1

julia> B = [Z Z Z O; w*y w*z-x*y x*z-y^2 Z];

julia> A = transpose(matrix(B));

julia> M = graded_cokernel(A)
Graded subquotient of graded submodule of R^2 with 2 generators
  1: e[1]
  2: e[2]
by graded submodule of R^2 with 4 generators
  1: w*y*e[2]
  2: (w*z - x*y)*e[2]
  3: (x*z - y^2)*e[2]
  4: e[1]

julia> FM1 = free_resolution(M)
Free resolution of M
R^2 <---- R^7 <---- R^8 <---- R^3 <---- 0
0         1         2         3         4

julia> betti_table(FM1)
degree: 0  1  2  3
------------------
    -1: -  1  -  -
     0: 2  -  -  -
     1: -  3  3  1
     2: -  3  5  2
------------------
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
degree: 0  1  2  3
------------------
    -1: -  1  -  -
     0: 2  -  -  -
     1: -  3  -  -
     2: -  -  4  2
------------------
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
degree: 0  1  2  3
------------------
     0: 1  -  -  -
     1: -  3  -  -
     2: -  -  4  2
------------------
 total: 1  3  4  2

julia> matrix(map(FM3, 1))
[-x*z + y^2]
[-w*z + x*y]
[       w*y]

```
"""
function free_resolution(M::SubquoModule{T}; 
                         length::Int = 0,
                         algorithm::Symbol = T <: MPolyRingElem ? :fres : :sres) where {T <: Union{MPolyRingElem, MPolyQuoRingElem}}

  coefficient_ring(base_ring(M)) isa AbstractAlgebra.Field ||
      error("Must be defined over a field.")

  if T <: MPolyQuoRingElem
    !iszero(length) || error("Specify a length up to which a free resolution should be computed")
  end

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
    if algorithm == :mres && is_graded(cc[-1]) # TODO: include local case once :mres is available there
        set_attribute!(cc, :minimal=>true)
    end
    _extend_free_resolution_to_the_left_by_zeros(cc, length)
    return FreeResolution(cc)
  end

  singular_free_module  = singular_module(ambient_free_module(kernel_entry))
  singular_kernel_entry = Singular.Module(base_ring(singular_free_module),
                              [singular_free_module(repres(g)) for g in gens(kernel_entry)]...)

  #= This is the single computational hard part of this function =#
  if algorithm == :fres && T <: MPolyRingElem
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
  elseif algorithm == :sres && T <: MPolyQuoRingElem
    gbpres = Singular.std(singular_kernel_entry)
    res = Singular.sres(gbpres, length)
  else
    error("Unsupported algorithm $algorithm")
  end

  if algorithm == :sres
    ub = findfirst(i -> is_zero(res[i]), 2:Singular.length(res)+1)
    slen = ub == nothing ? Singular.length(res) : ub
  else
    slen = iszero(res[Singular.length(res)+1]) ? Singular.length(res) : Singular.length(res)+1
  end
  if length == 0 || slen < length
    cc_complete = true
  end

  if length != 0
    slen =  slen > length ? length : slen #TODO check why the numbers may differ (Singular.jl)
  end

  #= Add maps from free resolution computation, start with second entry
   = due to inclusion of presentation(M) at the beginning. =#
  j   = 1
  while j <= slen
    if is_graded(M)
      codom = domain(maps[1])
      rk    = Singular.ngens(res[j])
      SM    = Oscar.SubModuleOfFreeModule(codom, res[j])
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
  if cc_complete == true
    _extend_free_resolution_to_the_left_by_zeros(cc, length)
  end
  if algorithm == :mres && is_graded(cc[-1]) # TODO: include local case once :mres is available there
    set_attribute!(cc, :minimal=>true)
  end
  cc.fill     = _extend_free_resolution
  cc.complete = cc_complete
  set_attribute!(cc, :show => Hecke.pres_show, :free_res => M)
  set_attribute!(cc, :algorithm, algorithm)
  return FreeResolution(cc)
end

function free_resolution(M::SubquoModule{T}; length::Int=0) where {T<:RingElem}
  # This generic code computes a free resolution in a lazy way.
  # We start out with a presentation of M and implement
  # an iterative fill function to compute every higher term
  # on request.
  R = base_ring(M)
  p = presentation(M)
  p.fill = function(C::Hecke.ComplexOfMorphisms, k::Int)
    min_index = first(chain_range(C))
    target_index = length == 0 ? k-1 : max(min_index - (length - 1), k-1)

    for i in min_index:target_index
      N = domain(map(C, i))
      if iszero(N)
        C.complete = true
        phi = hom(N, N, elem_type(N)[]; check=false)
        pushfirst!(C.maps, phi)
        continue
      end
      K, inc = kernel(map(C, i))
      nz = findall(!is_zero, gens(K))
      F = FreeMod(R, length(nz))
      phi = hom(F, C[i], iszero(length(nz)) ? elem_type(C[i])[] : inc.(gens(K)[nz]); check=false)
      pushfirst!(C.maps, phi)
      if length != 0 && abs(i - min_index) + 1 >= length
        C.complete = false
        break
      end
    end
    return first(C.maps)
  end
  return FreeResolution(p)
end

@doc raw"""
    free_resolution(M::SubquoModule{T};
        length::Union{Int, PosInf} = 0,
        algorithm = :generic) where {T <: Union{MPolyLocRingElem, MPolyQuoLocRingElem}}

Return a free resolution of `M`.

If `length != 0`, then the free resolution is only computed up to the `length`-th free module.
The only current option for `algorithm` is `:generic`. The name `:generic` indicates that the
method makes use of the generic function `free_resolution_via_kernels`.

# Examples
```jldoctest
julia> R, x = polynomial_ring(QQ, :x => 1:12);

julia> M = matrix(R, 3, 4, x)
[x[1]    x[2]    x[3]    x[4]]
[x[5]    x[6]    x[7]    x[8]]
[x[9]   x[10]   x[11]   x[12]]

julia> U1 = complement_of_point_ideal(R, zeros(Int, 12));

julia> RL1, _ = localization(R, U1);

julia> I1 = ideal(RL1, minors(M, 3));

julia> M1 = quotient_ring_as_module(I1);

julia> FM1 = free_resolution(M1)
Free resolution of M1
RL1^1 <---- RL1^4 <---- RL1^3 <---- 0
0           1           2           3

julia> minimize(FM1)
Free resolution of M1
RL1^1 <---- RL1^4 <---- RL1^3 <---- 0
0           1           2           3

julia> U2 = complement_of_point_ideal(R, vcat([1], zeros(Int, 11)));

julia> RL2, _ = localization(R, U2);

julia> I2 = ideal(RL2, minors(M, 3));

julia> M2 = quotient_ring_as_module(I2);

julia> FM2 = free_resolution(M2)
Free resolution of M2
RL2^1 <---- RL2^4 <---- RL2^3 <---- 0
0           1           2           3

julia> minimize(FM2)
Free resolution of M2
RL2^1 <---- RL2^3 <---- RL2^2 <---- 0
0           1           2           3

```

"""
function free_resolution(M::SubquoModule{T};
                         length::Union{Int, PosInf} = 0,
                         algorithm = :generic) where {T <: Union{MPolyLocRingElem, MPolyQuoLocRingElem}}
  if T <: MPolyQuoLocRingElem
    !iszero(length) || error("Specify a length up to which a free resolution should be computed")
  end
  if length == 0
    length = -1
  end
  if algorithm == :generic
    F = free_resolution_via_kernels(M, length)
  else
    error("Unsupported algorithm $algorithm")
  end
  set_attribute!(F.C, :algorithm, algorithm)
  return F
end

@doc raw"""
    free_resolution_via_kernels(M::SubquoModule, limit::Int = -1)

Return a free resolution of `M`.

If `limit != -1`, the free resolution
is only computed up to the `limit`-th free module.

!!! note
    This is a generic function which can be applied whenever kernels of homomorphisms between free modules can be computed.
# Examples

```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z]);

julia> U = complement_of_point_ideal(R, [0, 0, 0]);

julia> RL, _ = localization(R, U);

julia> F = free_module(RL, 1)
Free module of rank 1 over localization of R at complement of maximal ideal of point (0, 0, 0)

julia> A = RL[x; y]
[x]
[y]

julia> B = RL[x^2; y^3; z^4]
[x^2]
[y^3]
[z^4]

julia> M = subquotient(F, A, B)
Subquotient of submodule with 2 generators
  1: x*e[1]
  2: y*e[1]
by submodule with 3 generators
  1: x^2*e[1]
  2: y^3*e[1]
  3: z^4*e[1]

julia> FM = free_resolution_via_kernels(M, 2)
Free resolution of M
RL^2 <---- RL^5 <---- RL^4
0          1          2

julia> is_complete(FM)
false

julia> FM[4]
Free module of rank 0 over localization of R at complement of maximal ideal of point (0, 0, 0)

julia> FM
Free resolution of M
RL^2 <---- RL^5 <---- RL^4 <---- RL^1 <---- 0
0          1          2          3          4

julia> is_complete(FM)
true

```

```jldoctest
julia> R, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z]);

julia> F = graded_free_module(R, 1)
Graded free module R^1([0]) of rank 1 over R

julia> A = R[x; y]
[x]
[y]

julia> B = R[x^2; y^3; z^4]
[x^2]
[y^3]
[z^4]

julia> M = subquotient(F, A, B)
Graded subquotient of graded submodule of F with 2 generators
  1: x*e[1]
  2: y*e[1]
by graded submodule of F with 3 generators
  1: x^2*e[1]
  2: y^3*e[1]
  3: z^4*e[1]

julia> FM = free_resolution_via_kernels(M, 2)
Free resolution of M
R^2 <---- R^5 <---- R^4
0         1         2

julia> FM[4]
Graded free module R^0 of rank 0 over R

julia> FM
Free resolution of M
R^2 <---- R^5 <---- R^4 <---- R^1 <---- 0
0         1         2         3         4

julia> betti(FM)
degree: 0  1  2  3
------------------
     1: 2  2  -  -
     2: -  1  -  -
     3: -  -  1  -
     4: -  2  2  -
     5: -  -  1  -
     6: -  -  -  1
------------------
 total: 2  5  4  1

```
"""
function free_resolution_via_kernels(M::SubquoModule, limit::Int = -1)
  p = presentation(M)
  mp = [map(p, j) for j in Hecke.map_range(p)]  
  while true
    if limit != -1 && length(mp) > limit+1
      break
    end
    k, mk = kernel(mp[1])
    nz = findall(!is_zero, gens(k))
    if length(nz) == 0
      if is_graded(domain(mp[1]))
        h = graded_map(domain(mp[1]), Vector{elem_type(domain(mp[1]))}(); check=false)
      else
        Z = FreeMod(base_ring(M), 0)
        h = hom(Z, domain(mp[1]), Vector{elem_type(domain(mp[1]))}(); check=false)
      end
      insert!(mp, 1, h)
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
  C.fill = _extend_free_resolution_via_kernels
  C.complete = false
  if is_zero(domain(C.maps[begin]))
    _extend_free_resolution_to_the_left_by_zeros(C, limit)
    C.complete = true
  end
  #set_attribute!(C, :show => free_show, :free_res => M) # doesn't work
  set_attribute!(C, :show => Hecke.pres_show, :free_res => M)
  return FreeResolution(C)
end

function _extend_free_resolution_via_kernels(cc::Hecke.ComplexOfMorphisms, idx::Int)
  r = Hecke.map_range(cc)

  if idx < last(r)
    for j = last(r)-1:-1:idx
      cod =  codomain(last(cc.maps))
      phi = id_hom(cod)
      push!(cc.maps, phi)
      cc.seed = idx-1
    end
    return last(cc.maps)
  end

  len_missing = idx - first(r)
  @assert len_missing > 0
  if cc.complete == true
     _extend_free_resolution_to_the_left_by_zeros(cc, idx)
    return map(cc, first(r))
  end
  mp = cc.maps
  R = base_ring(domain(mp[1]))
  while true
    if idx != -1 && length(mp) > idx+1
      break
    end
    k, mk = kernel(mp[1])
    nz = findall(!is_zero, gens(k))
    if length(nz) == 0
      if is_graded(domain(mp[1]))
        h = graded_map(domain(mp[1]), Vector{elem_type(domain(mp[1]))}(); check=false)
      else
        Z = FreeMod(R, 0)
        h = hom(Z, domain(mp[1]), Vector{elem_type(domain(mp[1]))}(); check=false)
      end
      insert!(mp, 1, h)
      break
    end
    if is_graded(codomain(mk))
      g = graded_map(codomain(mk), collect(k.sub.gens)[nz]; check=false)
    else
      F = FreeMod(R, length(nz))
      g = hom(F, codomain(mk), collect(k.sub.gens)[nz]; check=false)
    end
    insert!(mp, 1, g)
  end
  if is_zero(domain(first(cc.maps)))
    _extend_free_resolution_to_the_left_by_zeros(cc, idx)
    cc.complete=true
  end
  return first(cc.maps)
end


@doc raw"""
    free_resolution(I::Ideal{T};
        length::Int = 0,
        algorithm::Symbol = I isa MPolyIdeal ? :fres : :sres) where {T <: Union{MPolyRingElem, MPolyQuoRingElem}}

Consider `I` as a subquotient module over its base ring, and return a free resolution of that module.

If `length != 0`, the free resolution is only computed up to the `length`-th free module.
At the moment, options for `algorithm` are `:fres`, `nres`, and `:mres`. With `:mres`,
in the positively graded case, a minimal free resolution is returned.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, [:x, :y]);

julia> I = ideal(R, [x, y])^2;

julia> fI = free_resolution(I)
Free resolution of I
R^3 <---- R^2 <---- 0
0         1         2

julia> AC = augmented_complex(fI)
0 <---- I <---- R^3 <---- R^2 <---- 0

julia> AC[-1]
Submodule with 3 generators
  1: x^2*e[1]
  2: x*y*e[1]
  3: y^2*e[1]
represented as subquotient with no relations

```

```
julia> R, (w, x, y, z) = graded_polynomial_ring(QQ, [:w, :x, :y, :z]);

julia> S, _ = quo(R, ideal(R, [w+x+y+z]));

julia> I = ideal(S, [w, x, y, z]);

julia> FI1 = free_resolution(I, length = 3, algorithm = :nres)
Free resolution of I
S^4 <---- S^6 <---- S^3 <---- 0
0         1         2         3

julia> FI2 = free_resolution(I, length = 3, algorithm = :mres)
Free resolution of I
S^3 <---- S^3 <---- S^1 <---- 0
0         1         2         3

```
"""
function free_resolution(I::Ideal{T};
                         length::Int = 0,
			 algorithm::Symbol = I isa MPolyIdeal ? :fres : :sres) where {T <: Union{MPolyRingElem, MPolyQuoRingElem}}

  S = ideal_as_module(I)
  n = AbstractAlgebra.get_name(I)
  if n !== nothing
    AbstractAlgebra.set_name!(S, n)
  end
  return free_resolution(S, length = length, algorithm = algorithm)
end

@doc raw"""
    free_resolution(I::Ideal{T};
        length::Int = 0,
        algorithm = :generic) where {T <: Union{MPolyLocRingElem, MPolyQuoLocRingElem}}

Consider `I` as a subquotient module over its base ring, and return a free resolution of that module.

If `length != 0`, then the free resolution is only computed up to the `length`-th free module.
The only current option for `algorithm` is `:generic`. The name `:generic` indicates that the
method makes use of the generic function `free_resolution_via_kernels`.

# Examples
```jldoctest
julia> R, x = polynomial_ring(QQ, :x => 1:12);

julia> M = matrix(R, 3, 4, x)
[x[1]    x[2]    x[3]    x[4]]
[x[5]    x[6]    x[7]    x[8]]
[x[9]   x[10]   x[11]   x[12]]

julia> U = complement_of_point_ideal(R, zeros(Int, 12));

julia> A, _ = quo(R, ideal(R, [gens(R)[1]]))
(Quotient of multivariate polynomial ring by ideal (x[1]), Map: R -> A)

julia> AL, _ = localization(A, U);

julia> I = ideal(AL, minors(M, 3));

julia> FI = free_resolution(I, length = 3)
Free resolution of I
AL^4 <---- AL^7 <---- AL^4 <---- 0
0          1          2          3

```
"""
function free_resolution(I::Ideal{T};
             length::Int = 0, algorithm = :generic) where {T <: Union{MPolyLocRingElem, MPolyQuoLocRingElem}}
  S = ideal_as_module(I)
  n = AbstractAlgebra.get_name(I)
  if n !== nothing
    AbstractAlgebra.set_name!(S, n)
  end
  return free_resolution(S, length = length, algorithm = algorithm)
end

@doc raw"""
    free_resolution(Q::MPolyQuoRing{T};
        length::Int = 0,
        algorithm::Symbol = T <: MPolyRingElem ? :fres : :sres) where {T <: Union{MPolyRingElem, MPolyQuoRingElem}}

Consider `Q` as a subquotient module over its base ring, and return a free resolution of that module.

If `length != 0`, the free resolution is only computed up to the `length`-th free module.
At the moment, options for `algorithm` are `:fres`, `:nres`, and `:mres`. With `:nres` or `:mres`,
in the positively graded case, a minimal free resolution is returned.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, [:x, :y]);

julia> I = ideal(R, [x, y])^2;

julia> Q, _ = quo(R, I);

julia> fQ = free_resolution(Q)
Free resolution of Q
R^1 <---- R^3 <---- R^2 <---- 0
0         1         2         3

julia> AC = augmented_complex(fQ)
0 <---- Q <---- R^1 <---- R^3 <---- R^2 <---- 0

julia> AC[-1]
Subquotient of submodule with 1 generator
  1: e[1]
by submodule with 3 generators
  1: x^2*e[1]
  2: x*y*e[1]
  3: y^2*e[1]

```

```
julia> R, (w, x, y, z) = graded_polynomial_ring(QQ, [:w, :x, :y, :z]);

julia> S, _ = quo(R, ideal(R, [w+x+y+z]));

julia> I = ideal(S, [w, x, y, z]);

julia> Q, _ = quo(S, I)
(Quotient of multivariate polynomial ring by ideal (w + x + y + z, w, x, y, z), Hom: S -> Q)

julia> FI1 = free_resolution(Q, length = 5, algorithm = :nres)
Free resolution of Q
R^1 <---- R^4 <---- R^6 <---- R^4 <---- R^1 <---- 0
0         1         2         3         4         5

julia> FI2 = free_resolution(Q, length = 5, algorithm = :mres)
Free resolution of Q
R^1 <---- R^4 <---- R^6 <---- R^4 <---- R^1 <---- 0
0         1         2         3         4         5

```
"""
function free_resolution(Q::MPolyQuoRing{T};
             length::Int = 0,
             algorithm::Symbol = T <: MPolyRingElem ? :fres : :sres) where {T <: Union{MPolyRingElem, MPolyQuoRingElem}}

  q = quotient_ring_as_module(Q)
  n = AbstractAlgebra.get_name(Q)
  if n !== nothing
    AbstractAlgebra.set_name!(q, n)
  end
  return free_resolution(q, length = length, algorithm = algorithm)
end

@doc raw"""
    free_resolution(Q::MPolyQuoLocRing{T};
        length::Int = 0, algorithm = :generic) where T
Consider `Q` as a subquotient module over its base ring, and return a free resolution of that module.

If `length != 0`, then the free resolution is only computed up to the `length`-th free module.
The only current option for `algorithm` is `:generic`. The name `:generic` indicates that the
method makes use of the generic function `free_resolution_via_kernels`.

# Examples
```jldoctest
julia> R, x = polynomial_ring(QQ, :x => 1:12);

julia> M = matrix(R, 3, 4, x)
[x[1]    x[2]    x[3]    x[4]]
[x[5]    x[6]    x[7]    x[8]]
[x[9]   x[10]   x[11]   x[12]]

julia> U = complement_of_point_ideal(R, zeros(Int, 12));

julia> RL, _ = localization(R, U);

julia> I = ideal(RL, minors(M, 3));

julia> Q, _ = quo(RL, I);

julia> FQ = free_resolution(Q)
Free resolution of Q
RL^1 <---- RL^4 <---- RL^3 <---- 0
0          1          2          3

```
"""
function free_resolution(Q::MPolyQuoLocRing{T};
             length::Int = 0, algorithm = :generic) where T
  q = quotient_ring_as_module(Q)
  n = AbstractAlgebra.get_name(Q)
  if n !== nothing
    AbstractAlgebra.set_name!(q, n)
  end
  return free_resolution(q, length = length, algorithm = algorithm)
end

@doc raw"""
    augmented_complex(F::FreeResolution)

Return the augmented complex of `F`.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z]);

julia> A = R[x; y];

julia> B = R[x^2; x*y; y^2; z^4];

julia> M = SubquoModule(A, B);

julia> fr = free_resolution(M)
Free resolution of M
R^2 <---- R^6 <---- R^6 <---- R^2 <---- 0
0         1         2         3         4

julia> augmented_complex(fr)
0 <---- M <---- R^2 <---- R^6 <---- R^6 <---- R^2 <---- 0

```
"""
function augmented_complex(F::FreeResolution)
  ub = findfirst(i -> is_zero(F.C[i]), 2:first(map_range(F.C)))
  if isnothing(ub)
    AC = F.C[-2:first(map_range(F.C))]
  else
    AC = F.C[-2:ub+1]
  end
  set_attribute!(AC, :show => Hecke.pres_show)
  return AC
end

@doc raw"""
    length(F::FreeResolution)

Return the length of `F`.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z]);

julia> A = R[x; y];

julia> B = R[x^2; x*y; y^2; z^4];

julia> M = SubquoModule(A, B);

julia> fr = free_resolution(M, length = 2)
Free resolution of M
R^2 <---- R^6 <---- R^6
0         1         2

julia> length(fr)
2

julia> fr = free_resolution(M)
Free resolution of M
R^2 <---- R^6 <---- R^6 <---- R^2 <---- 0
0         1         2         3         4

julia> length(fr)
3

```
"""
function length(F::FreeResolution)
  ub = findfirst(i -> is_zero(F.C[i]), 2:first(map_range(F.C)))
  isnothing(ub) && return first(map_range(F.C))
  return ub::Int
end

