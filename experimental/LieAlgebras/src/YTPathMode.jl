# Young Tableu Path Model

struct YTPathModel <: AbstractCrystal
  wt::WeightLatticeElem
end

struct YTPathModelElem <: AbstractCrystalElem
  parent::YTPathModel
  t::YoungTableau
end

function parent(t::YTPathModelElem)
  return t.parent
end

# Young Tableau path model implemenation

function (P::YTPathModel)(p::LSPathModelElem)
  @req P.wt == parent(p).wt "dominant weights must match"

  pp = dominant_path(P)
  i = Int.(word(max(p)))
  falpha!(pp, i, adapted_string(p, i))
  return pp
end

function yt_path_model(wt::WeightLatticeElem)
  return YTPathModel(wt)
end

function yt_path_model(R::RootSystem, wt::Vector{Int})
  return YTPathModel(WeightLatticeElem(R, wt))
end

function Base.show(io::IO, mime::MIME"text/plain", p::YTPathModelElem)
  show(io, mime, p.t)
end

function dominant_path(P::YTPathModel)
  rk = rank(root_system(P.wt))
  shape = zeros(Int, rk)
  shape[rk] = Int(P.wt[rk])
  for i in (rk - 1):-1:1
    shape[i] = shape[i + 1] + Int(P.wt[i])
  end

  return YTPathModelElem(
    P, young_tableau([fill(i, shape[i]) for i in 1:length(shape) if shape[i] > 0])
  )
end

function ealpha!(p::YTPathModelElem, i::Int)
  t = p.t.t

  i = (1, length(t[r]))
  m = 0
  h = 0
  for r in 1:length(t)
    for c in length(t[r]):-1:1
      if t[r][c] == i
        h += 1
      elseif t[r][c] == i + 1
        h -= 1
        if h <= m
          m = h
          i = (r, c)
        end
      end
    end
  end
end

function falpha!(p::YTPathModelElem, i::Int)
  t = p.t.t

  # find last time the global minimum is assumed
  k, l = 1, length(t[1]) + 1
  m = 0
  h = 0
  for r in 1:length(t)
    for c in length(t[r]):-1:1
      if t[r][c] == i
        h += 1
      elseif t[r][c] == i + 1
        h -= 1
      end

      if h <= m
        m = h
        k, l = r, c
      end
    end
  end

  if m < h
    if l == 1
      k += 1
      l = length(t[k]) + 1
    end
    t[k][l - 1] += 1
  end

  return p
end
