function exterior_power(M::SubquoModule, p::Int; cached::Bool=true)
  n = rank(ambient_free_module(M))
  R = base_ring(M)::base_ring_type(M)
  @req 0 <= p <= n "Exponent out of bounds"

  if cached
    powers = _exterior_powers(M)
    haskey(powers, p) && return powers[p]
  end

  if iszero(p)
    F = FreeMod(R, 1)
    result, _ = sub(F, [F[1]])
  else
    C = presentation(M)
    phi = map(C, 1)
    codomain(phi).S = function _get_symbol()
      return [Symbol"$e") for e in gens(M)]
    end
    result, mm = _exterior_power(phi, p)
  end

  function my_mult(u::Tuple{Vararg{SubquoModuleElem}})
    isempty(u) && return result[1] # only the case p=0
    @req all(x -> parent(x) === M, u) "elements must live in the same module"
    @req length(u) == p "need a $p-tuple of elements"
    return wedge(collect(u), parent=result)
  end
  function my_mult(u::SubquoModuleElem...)
    return my_mult(u)
  end

  function my_decomp(u::SubquoModuleElem)
    @req parent(u) === result "element does not belong to the correct module"
    k = findfirst(==(u), gens(result))
    @req !isnothing(k) "element must be a generator of the module"
    ind = combination(n, p, k)
    return Tuple(gen(M, i) for i in ind)
  end

  mult_map = MapFromFunc(Hecke.TupleParent(Tuple([zero(M) for f in 1:p])), result, my_mult, my_decomp)
  inv_mult_map = MapFromFunc(result, domain(mult_map), my_decomp, my_mult)

  set_attribute!(result, :multiplication_map, mult_map)
  set_attribute!(result, :wedge_pure_function, mult_map)
  set_attribute!(result, :wedge_generator_decompose_function, inv_mult_map)

  set_attribute!(result, :is_exterior_power, (M, p))
  cached && (_exterior_powers(M)[p] = (result, mult_map))

  # Set the variable names for printing
  orig_symb = ["$(e)" for e in ambient_representatives_generators(M)]
  new_symb = Symbol[]
  if iszero(p)
    new_symb = [Symbol("1")]
  else
    for ind in combinations(n, p)
      symb_str = orig_symb[ind[1]]
      for i in 2:p
        symb_str = symb_str * (is_unicode_allowed() ? "∧" : "^") * orig_symb[ind[i]]
      end
      push!(new_symb, Symbol(symb_str))
    end
  end

  symbols(result) = new_symb
  return result, mult_map
end

function _exterior_power(phi::FreeModuleHom, p::Int)
  F = codomain(phi)
  R = base_ring(F)
  rel = domain(phi)
  Fp, mult_map = exterior_power(F, p)
  Fq, _ = exterior_power(F, p-1)
  img_gens = [wedge(e, f) for e in gens(Fq) for f in phi.(gens(rel))]
  G = FreeMod(R, length(img_gens))
  psi = hom(G, Fp, img_gens)
  # TODO: the `cokernel` command does not have consistent output over all types of rings.
  return (base_ring(codomain(phi)) isa Union{MPolyLocRing, MPolyQuoLocRing} ? cokernel(psi)[1] : cokernel(psi)), mult_map
end
