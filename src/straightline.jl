## constants

const opmask    = 0xff00000000000000
const argmask   = 0x000000000fffffff
const inputmark = 0x0000000008000000
const tmpmark   = 0x0000000004000000
const argshift  = 28 # argmask == 2^argshift - 1


const showop = Dict{UInt64,Char}()

for (i, (op, unary, showchar)) in enumerate([(:uniplus      , true,  '+'),
                                             (:uniminus     , true,  '-'),
                                             (:plus         , false, '+'),
                                             (:minus        , false, '-'),
                                             (:times        , false, '*'),
                                             (:divide       , false, '/'),
                                             (:exponentiate , true,  '^')])
   isop = Symbol(:is, op)
   c = UInt64(i) << (2*argshift)
   if unary
      c |= 0x8000000000000000
   end
   @eval begin
      const $op = $c
      $isop(op) = op === $op
   end
   showop[c] = showchar
end


## raw manips

pack(op::UInt64, i, j) = op | ((i % UInt64) << argshift) | (j % UInt64)


## building SLPoly

function pushconst!(p::SLPoly{T}, c::T) where T
   push!(p.cs, c)
   l = lastindex(p.cs)
   @assert l < tmpmark
   l % UInt64
end

function pushop!(p::SLPoly, op, i, j=UInt64(0))
   @assert i <= argmask && j <= argmask
   push!(p.lines, pack(op, i, j))
   l = lastindex(p.lines)
   @assert l < tmpmark
   l % UInt64 | tmpmark
end
