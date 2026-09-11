struct Origami
  o::GapObj
  h::PermGroupElem
  v::PermGroupElem
  d::Int
end

struct Normal_stored_origami
  o::GapObj
  h::PermGroupElem
  v::PermGroupElem
  g::GAPGroup
end

# the bijection between the cycles of the bot/top permutation
# is given implicitly by using lists and considering the order
# of the cycles
struct CylinderDiagram
  bot::Vector{Vector{Int}}
  top::Vector{Vector{Int}}
  cycles_count::Int
  separatrix_count::Int
  function CylinderDiagram(bot::Vector{Vector{Int}}, top::Vector{Vector{Int}})
    max = 1
    for cycle in bot
      for i in cycle
        if i > max
          max = i
        end
      end
    end
    new(bot, top, length(bot), max + 1)
  end
end

# `cycles_count` and `separatrix_count` are derived from `bot`, so comparing
# `bot` and `top` is enough. Without this, containers fall back to `objectid`.
Base.:(==)(a::CylinderDiagram, b::CylinderDiagram) = a.bot == b.bot && a.top == b.top

function Base.hash(cd::CylinderDiagram, ha::UInt)
  return hash(cd.top, hash(cd.bot, ha))
end

struct CanonicalOrigamiKey
  h::Vector{Int}
  v::Vector{Int}
end

Base.:(==)(a::CanonicalOrigamiKey, b::CanonicalOrigamiKey) =
  a.h == b.h && a.v == b.v

function Base.hash(key::CanonicalOrigamiKey, ha::UInt)
  return hash(key.v, hash(key.h, ha))
end
