## constants

const opmask    = 0xff00000000000000
const argmask   = 0x000000000fffffff
const inputmark = 0x0000000008000000
const tmpmark   = 0x0000000004000000
const argshift  = 28 # argmask == 2^argshift - 1


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
