mutable struct WeymanCtx
  pfctx::Union{PushForwardCtx, ToricCtx, NewToricCtx, ToricCtxWithParams}
  cplx::AbsHyperComplex
  macro_modules::Dict
  weyman_inclusions::Dict{Tuple{Int, Int, Int}, Dict}
  weyman_projections::Dict{Tuple{Int, Int, Int}, Dict}

  function WeymanCtx(
      pfctx::Union{PushForwardCtx, ToricCtx, NewToricCtx, ToricCtxWithParams},
      cplx::AbsHyperComplex
    )
    @assert dim(cplx) == 1 "complex must be one-dimensional"
    @assert direction(cplx) == :cochain "complex must be cochain"
    return new(pfctx, cplx)
  end
end

pushforward_ctx(W::WeymanCtx) = W.pfctx
graded_complex(W::WeymanCtx) = W.cplx

mutable struct MacroMod
  weyman_ctx::WeymanCtx
  p::Int
  q::Int
  typ::Symbol
  function MacroMod(wctx::WeymanCtx, p::Int, q::Int, typ::Symbol)
    @assert (typ == :cochain || typ == :cohomology) "wrong typ"
    return new(wctx, p, q, typ)
  end
end

function get_macro_block!(wctx::WeymanCtx, p::Int, q::Int, typ::Symbol)
  @assert q >= 0
  if !isdefined(wctx, :macro_modules)
    wctx.macro_modules = Dict{Tuple{Int, Int, Symbol}, MacroMod}()
  end
  return get!(wctx.macro_modules::Dict{Tuple{Int, Int, Symbol}, MacroMod}, (p, q, typ)) do
    MacroMod(wctx, p, q, typ)
  end
end

mutable struct MacroVec
  macro_mod::MacroMod
  micro_vecs::Vector

  function MacroVec(mac_mod::MacroMod; check::Bool=false)
    return new(mac_mod)
  end
end

mutable struct MicroVec
  macro_mod::MacroMod
  j::Int
  e::Vector{Int}
  c::FreeModElem

  function MicroVec(mac_mod::MacroMod, j::Int, e::Vector{Int}, 
      c::FreeModElem;
      check::Bool=false
    )
    p, q = index(mac_mod)
    gc = graded_complex(mac_mod)
    F = gc[p]
    alpha = -degrees_of_generators(F)[j]
    ctx = pushforward_ctx(mac_mod)
    if typ(mac_mod) == :cochain
      @check parent(c) === ctx[e, alpha][-q] "wrong parent"
    elseif typ(mac_mod) == :cohomology
      @check parent(c) === simplified_strand(ctx, e, alpha)[-q] "wrong parent"
    end
    return new(mac_mod, j, e, c)
  end
  function MicroVec(mac_mod::MacroMod, j::Int, e::Vector{Int}; check::Bool=false)
    p, q = index(mac_mod)
    gc = graded_complex(mac_mod)
    F = gc[p]
    alpha = degrees_of_generators(F)[j]
    ctx = pushforward_ctx(mac_mod)
    c = typ(mac_mod) == :cochain ? zero(ctx[e, -alpha][-q]) : zero(simplified_strand(ctx, e, -alpha)[-q])
    return new(mac_mod, j, e, c)
  end
  function MicroVec(mac_mod::MacroMod, j::Int, e::Vector{Int}, k::Int; check::Bool=false)
    p, q = index(mac_mod)
    gc = graded_complex(mac_mod)
    F = gc[p]
    alpha = degrees_of_generators(F)[j]
    ctx = pushforward_ctx(mac_mod)
    c = typ(mac_mod) == :cochain ? ctx[e, -alpha][-q][k] : simplified_strand(ctx, e, -alpha)[-q][k]
    return new(mac_mod, j, e, c)
  end
end

weyman_ctx(v::MacroMod) = v.weyman_ctx
index(v::MacroMod) = (v.p, v.q)
function micro_vectors(v::MacroVec)
  if !isdefined(v, :micro_vecs)
    v.micro_vecs = Tuple{Int, MicroVec}[]
  end
  return v.micro_vecs::Vector{Tuple{Int, MicroVec}}
end
pushforward_ctx(v::MacroMod) = pushforward_ctx(weyman_ctx(v))
graded_complex(v::MacroMod) = graded_complex(weyman_ctx(v))
typ(v::MacroMod) = v.typ

macro_module(v::MacroVec) = v.macro_mod
weyman_ctx(v::MacroVec) = weyman_ctx(macro_module(v))
pushforward_ctx(v::MacroVec) = pushforward_ctx(macro_module(v))
typ(v::MacroVec) = typ(macro_module(v))
graded_complex(v::MacroVec) = graded_complex(macro_module(v))
index(v::MacroVec) = index(macro_module(v))

macro_module(v::MicroVec) = v.macro_mod
index(v::MicroVec) = v.j
macro_index(v::MicroVec) = index(macro_module(v))
typ(v::MicroVec) = typ(macro_module(v))
value(v::MicroVec) = v.c
direct_limit_index(v::MicroVec) = v.e
pushforward_ctx(v::MicroVec) = pushforward_ctx(macro_module(v))
graded_complex(v::MicroVec) = graded_complex(macro_module(v))

function degree(v::MicroVec)
  gc = graded_complex(v)
  j = index(v)
  p, q = macro_index(v)
  return -degrees_of_generators(gc[p])[j]
end


