#=
 Reference Clemens here


=#
function f4ncgb_version()
    ret_c = @ccall libf4ncgb.f4ncgb_version()::Cstring
    ret = unsafe_string(ret_c)
    return  ret
end

function f4ncgb_set_msg_printing(on::Bool)
  println("Setting message printing to ", on)
  @ccall libf4ncgb.f4ncgb_set_msg_printing(on::Cuchar)::Cvoid
end

f4ncgb_init() = @ccall libf4ncgb.f4ncgb_init()::Ptr{Cvoid}
f4ncgb_free(handle::Ptr{Cvoid}) = @ccall libf4ncgb.f4ncgb_free(handle::Ptr{Cvoid})::Cvoid
f4ncgb_get_state(handle::Ptr{Cvoid}) = @ccall libf4ncgb.f4ncgb_get_state(handle::Ptr{Cvoid})::Ptr{Cvoid}

function f4ncgb_add(handle::Ptr{Cvoid},
                    numerator::Int,
                    denominator::Int,
                    vars::Vector{UInt32})
    return @ccall libf4ncgb.f4ncgb_add(
      handle::Ptr{Cvoid},
      numerator::Clong,
      denominator::Clong,
      length(vars)::Csize_t,
      vars::Ptr{UInt32})::Cstring
end

function f4ncgb_add(handle::Ptr{Cvoid}, polynomial::FreeAssociativeAlgebraElem)
  for m in terms(polynomial)
    c = coeff(m, 1)
    d = denominator(c)
    n = numerator(c)
    f4ncgb_add(handle, Int(n), Int(d), UInt32.(m.exps[1]))
  end
  f4ncgb_end_poly(handle)
end

f4ncgb_end_poly(handle::Ptr{Cvoid}) = @ccall libf4ncgb.f4ncgb_end_poly(handle::Ptr{Cvoid})::Cstring

function f4ncgb_set_blocks(handle::Ptr{Cvoid}, block_lengths::Vector{UInt32})
    return @ccall libf4ncgb.f4ncgb_set_blocks(
      handle::Ptr{Cvoid},
      length(block_lengths)::Csize_t,
      block_lengths::Ptr{UInt32})::Cstring
end

function f4ncgb_set_nvars(handle::Ptr{Cvoid}, nvars::UInt32) 
  @ccall libf4ncgb.f4ncgb_set_nvars(handle::Ptr{Cvoid}, nvars::Csize_t)::Cstring
end

function f4ncgb_set_characteristic(handle::Ptr{Cvoid}, characteristic::UInt32)
  return @ccall libf4ncgb.f4ncgb_set_characteristic(
    handle::Ptr{Cvoid},
    characteristic::Csize_t)::Cstring
end

function f4ncgb_set_maxiter(handle::Ptr{Cvoid}, maxiter::UInt32)
  return @ccall libf4ncgb.f4ncgb_set_maxiter(
    handle::Ptr{Cvoid},
    maxiter::Csize_t)::Cstring
end

function f4ncgb_set_maxdeg(handle::Ptr{Cvoid}, maxdeg::UInt32)
  return @ccall libf4ncgb.f4ncgb_set_maxdeg(
    handle::Ptr{Cvoid},
    maxdeg::Csize_t)::Cstring
end

function f4ncgb_set_threads(handle::Ptr{Cvoid}, threads::UInt32)
  return @ccall libf4ncgb.f4ncgb_set_threads(
    handle::Ptr{Cvoid},
    threads::Csize_t)::Cstring
end

function f4ncgb_set_output_file(handle::Ptr{Cvoid}, filename::String)
  return @ccall libf4ncgb.f4ncgb_set_output_file(
    handle::Ptr{Cvoid},
    filename::Cstring)::Cstring
end

function f4ncgb_set_proof_file(handle::Ptr{Cvoid}, filename::String)
  return @ccall libf4ncgb.f4ncgb_set_proof_file(
    handle::Ptr{Cvoid},
    filename::Cstring)::Cstring
end

function f4ncgb_set_expanded_proof(handle::Ptr{Cvoid}, expanded::Bool)
  return @ccall libf4ncgb.f4ncgb_set_expanded_proof(
    handle::Ptr{Cvoid},
    expanded::Cuchar)::Cstring
end

function f4ncgb_set_tracer(handle::Ptr{Cvoid}, trace::Bool)
  return @ccall libf4ncgb.f4ncgb_set_tracer(
    handle::Ptr{Cvoid},
    trace::Cuchar)::Cstring
end

mutable struct f4ncgb_polys_helper
  gens::Vector{FreeAssociativeAlgebraElem{QQFieldElem}}
  current_poly::FreeAssociativeAlgebraElem{QQFieldElem}
  parent::FreeAssociativeAlgebra{QQFieldElem}
  function f4ncgb_polys_helper(ring::FreeAssociativeAlgebra{QQFieldElem})
    r = new()
    r.gens = FreeAssociativeAlgebraElem{QQFieldElem}[]
    r.parent = ring
    r.current_poly = zero(r.parent)
    return r
  end
end 

function add_cb(a::f4ncgb_polys_helper, monomial::FreeAssociativeAlgebraElem)
  current_poly += monomial
end

function parent(a::f4ncgb_polys_helper)
  return a.parent
end

function zero(a::f4ncgb_polys_helper)
  return zero(parent(a))
end

function add_cb(pa::Ptr{Nothing},
                pnumerator::Ptr{BigInt},
                pdenominator::Ptr{BigInt},
                varcount::Csize_t,
                pvars:: Ptr{UInt32})

  a = unsafe_pointer_to_objref(convert(Ptr{f4ncgb_polys_helper}, pa))
  P = parent(a)
  vars = unsafe_wrap(Array, pvars, varcount)
  numerator = unsafe_load(pnumerator)
  denominator = unsafe_load(pdenominator)
  monomial = P(QQ(numerator,denominator))
  for i in vars
    monomial *= P[Int(i)]
  end
  a.current_poly += monomial
  return nothing

end

function end_poly_cb(pa::Ptr{Nothing})
  a = unsafe_pointer_to_objref(convert(Ptr{f4ncgb_polys_helper}, pa))
  push!(a.gens, a.current_poly)
  a.current_poly = zero(parent(a))
  return nothing
end

function Base.show(io::IO, a::f4ncgb_polys_helper)
  print(io, "f4ncgb polynomial helper with $(length(a.gens)) elements")
end

function f4ncgb_solve(handle::Ptr{Cvoid}, ideal::f4ncgb_polys_helper)
  add_cb_c = @cfunction(add_cb, Cvoid, (Ptr{Cvoid}, Ptr{BigInt}, Ptr{BigInt}, Csize_t, Ptr{UInt32}))
  end_poly_cb_c = @cfunction(end_poly_cb, Cvoid, (Ptr{Cvoid},))

  return @ccall libf4ncgb.f4ncgb_solve(
    handle::Ptr{Cvoid},
    pointer_from_objref(ideal)::Ptr{Cvoid},
    add_cb_c::Ptr{Cvoid},
    end_poly_cb_c::Ptr{Cvoid})::Cint
end

function f4ncgb_groebner(x::Vector{<:FreeAssociativeAlgebraElem{QQFieldElem}})
  @req length(x) > 0 "Fix your shit"
  r = parent(first(x))
  handle = f4ncgb_init()
  f4ncgb_set_nvars(handle, UInt32(ngens(r)))
  f4ncgb_set_blocks(handle, UInt32[ngens(r)])
  f4ncgb_set_threads(handle, UInt32(1))
  f4ncgb_add.(Ref(handle), x)
  userdata = f4ncgb_polys_helper(r)
  f4ncgb_solve(handle, userdata)
  f4ncgb_free(handle)
  return userdata.gens
end

#=
S1 = quantum_symmetric_group(4);
Oscar.f4ncgb_groebner(gens(S1))
=#
