@doc raw"""
    root_symbols(n::Int) 
    
Return the list of all ADE root symbols up to permutation
whose combined rank is ``n``.

# Examples
```jldoctest
julia> root_symbols(3)
3-element Vector{Any}:
 [(:A, 3)]
 [(:A, 2), (:A, 1)]
 [(:A, 1), (:A, 1), (:A, 1)]

```
"""
function root_symbols(n::Int)
  result = []
  for p in partitions(n)
    tmp = Vector{Tuple{Symbol,Int}}[]
    for i in p
      symb = Tuple{Symbol,Int}[]
      push!(symb, (:A, i))
      if 4<=i
        push!(symb, (:D, i))
      end 
      if 6<=i<=8 
        push!(symb, (:E, i))
      end 
      push!(tmp, symb)
    end
    for i in Base.product(tmp...)
      g = collect(i)
      push!(result, g)
    end
  end 
  return result
end 

"""
    root_lattices(n::Int) -> Vector{ZZLat}
    
Return the list of all root lattices up to isomorphism whose rank is ``n``.
"""
function root_lattices(n::Int)
  result = []
  for p in partitions(n)
    tmp = Vector{QQMatrix}[]
    for i in p
      lat = QQMatrix[]
      push!(lat, Hecke._root_lattice_A(i))
      if 4<=i
        push!(lat, Hecke._root_lattice_D(i))
      end 
      if 6<=i<=8 
        push!(lat, Hecke._root_lattice_E(i))
      end 
      push!(tmp, lat)
    end
    for i in Base.product(tmp...)
      g = diagonal_matrix(collect(i))
      push!(result, integer_lattice(gram=g))
    end
  end 
  return result
end

# slooow
function _overlattice_orbits(L::ZZLat; even=true)
  d = ZZ(det(L))
  D = discriminant_group(L)
  idD = hom(D,D,gens(D))
  G = image_in_Oq(L)[1]
  orders = [i for i in divisors(d) if divides(d,i^2)[1]]
  result = ZZLat[]
  for ord in orders 
    @show ord, D
    b, l, p = is_prime_power_with_data(ord)
    if b && is_elementary(D, p)
      sg = Oscar._subgroups_orbit_representatives_and_stabilizers_elementary(idD, G, ord, p)    
    else 
      sg = Oscar._subgroups_orbit_representatives_and_stabilizers(idD, G, ord)    
    end
    for (S,_) in sg 
      M = cover(domain(S))
      if !is_integral(M) || (even && !is_even(M))
        continue
      end 
      push!(result,M)        
    end
  end
  return result
end

function root_overlattices(n::Int)
  result = ZZLat[]
  for R in root_lattices(n)
    for S in _overlattice_orbits(R)
      # only add the ones not adding new roots 
      RS = root_sublattice(S)
      if RS == R # this is terribly inefficient
        push!(result,S)
      end 
      if S != RS
        @show "new"
      end
    end 
  end 
  return result
end 

function orbit_representatives_and_stabilizers_magma(G::MatrixGroup{E}, k::Int) where E <: FinFieldElem
  g = gens(G) 
  d = degree(G)
  F = base_ring(G)
  V = vector_space(F, d)
  p, e = characteristic(F), degree(F)
  str = "K := GF($p, $e); G := MatrixGroup<$d, K | "
  for m in g
    mm = "[" * split(string([m[i,j] for i in 1:nrows(m) for j in 1:ncols(m)]), '[')[2]
    str *= mm
    str *= ", "
  end
  str = str[1:end-2]*" >; O := OrbitsOfSpaces(G, $k); S1 := Sprint([[[[M[j,k] : k in [1..NumberOfColumns(M)]] : j in [1..NumberOfRows(M)]] : M in Generators(Stabilizer(G, L[2]       ))] : L in O]); S2 := Sprint([[[M[k] : k in [1..NumberOfColumns(M)]] : M in Basis(L[2])] : L in O]); [S1, S2]"
  o = MagmaCall.interact() do stdout
    MagmaCall.putcmd(stdout, str)
    MagmaCall.readtotoken(String, stdout, missing)
    #|> Meta.parse |> eval
  end
  return o
  L1, L2 = o
  stabs = Vector{elem_type(G)}[elem_type(G)[G(matrix(F, d, d, reduce(vcat, v))) for v in bas] for bas in L1]
  stabs = [sub(G, bas)[1] for bas in stabs]
  orb = Vector{elem_type(V)}[elem_type(V)[V(F.(v)) for v in bas] for bas in L2]
  orb = [sub(V, bas)[1] for bas in orb]
  return [(orb[i], stabs[i]) for i in 1:length(orb)]
end   
