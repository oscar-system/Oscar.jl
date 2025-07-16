struct LeechPair
  leech::ZZLat
  G::MatrixGroup{QQFieldElem, QQMatrix}
  id::Int
  
  _info::Dict
  
  function LeechPair(id::Int, leech::ZZLat, gg::Vector{QQMatrix})
    path = @__DIR__
    l = readlines(path*"/query.csv")[id]
    vl = split(l, "|")
    info = Dict(:rank => parse(Int, vl[2]),
                :order => Hecke.parse(ZZRingElem, vl[3]),
                :alpha => parse(Int, vl[4]),
                :icoinvbar => parse(Int, vl[5]),
                :iinvbar => parse(Int, vl[6]),
                :ind => parse(Int, vl[7]),
                :hinv => parse(Int, vl[8]),
                :N => parse(Int, vl[9]),
                :type => vl[10]
           )
    return new(leech, matrix_group(gg), id, info)
  end
end
