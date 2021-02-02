module BlockSys

using Oscar

struct BlockSystems
  n::Int
  l::Int
  cur::Array{Array{Int, 1}, 1}
  function BlockSystems(n::Int, l::Int)
    @assert n % l == 0
    return new(n, l, [collect((i-1)*l+1:i*l) for i=1:divexact(n, l)])
  end
end

function Base.iterate(B::BlockSystems)
  return B.cur, deepcopy(B.cur)
end

function Base.iterate(B::BlockSystems, st::Array{Array{Int, 1}})
  if B.l==1||B.l==B.n
    return nothing
  end
  i = length(B.cur)-1
  while true
    j = B.l
    while true
      if st[i][j] < B.n - B.l + j
        st[i][j] += 1
        free = Set(1:B.n)
        for l=1:i-1
          setdiff!(free, st[l])
        end
        if !(st[i][j] in free) 
          continue
        end
        if length(intersect(free, Set(st[i][j]+1:B.n)))<B.l-j
          continue
        end
        setdiff!(free, st[i][1:j])
        while j < B.l
          j += 1
          I = intersect(free, Set(st[i][j-1]:B.n))
          if isempty(I)
            break
          end
          st[i][j] = minimum(I)
          pop!(free, st[i][j])
        end
        i += 1
        while i <= length(st)
          for j=1:B.l
            st[i][j] = minimum(free)
            pop!(free, st[i][j])
          end
          i += 1
        end
        return deepcopy(st), st
      end
      j -= 1
      if j == 1
        i -= 1
        i == 0 && return nothing
        break
      end
    end
  end
end
Base.IteratorSize(::BlockSystems) = Base.SizeUnknown()


end
