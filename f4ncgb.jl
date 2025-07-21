
using f4ncgb_jll

export f4ncgb_version,
       f4ncgb_set_msg_printing,
       f4ncgb_init,
       f4ncgb_add,
       f4ncgb_free,
       f4ncgb_get_state,
       f4ncgb_end_poly,
       f4ncgb_set_blocks,
       f4ncgb_set_nvars,
       f4ncgb_set_characteristic,
       f4ncgb_set_maxiter,
       f4ncgb_set_maxdeg,
       f4ncgb_set_threads,
       f4ncgb_set_output_file


#libf4ncgb = f4ncgb_jll.libf4ncgb_path

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
#=
handle = f4ncgb_init()

f4ncgb_add(handle, 1, 2, UInt32[1, 2, 3])
f4ncgb_set_output_file(handle, "output.txt")


=#
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

f4ncgb_end_poly(handle::Ptr{Cvoid}) = @ccall libf4ncgb.f4ncgb_end_poly(handle::Ptr{Cvoid})::Cstring

function f4ncgb_set_blocks(handle::Ptr{Cvoid}, block_lengths::Vector{UInt32})
    return @ccall libf4ncgb.f4ncgb_set_blocks(
      handle::Ptr{Cvoid},
      blocks::Csize_t,
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

#=
S1 = quantum_symmetric_group(4);
x1 = gens(S1)[1]
typeof(x1)
=#




