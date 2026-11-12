################################################################################
# Wikidata items for serialized types
#
# Maps a serialized type name (see `encode_type`) to the Wikidata item of the
# mathematical object it represents. Entries are curated by hand and checked
# against their Wikidata label, since searching Wikidata by name is ambiguous
# (e.g. "ideal" also matches a piano composition). Types without an exact
# Wikidata item, such as `Polyhedron` (not necessarily bounded or 3-dimensional),
# are deliberately left out.

const wikidata_qids = Dict{String, String}(
  "Base.Int" => "Q12503",               # integer
  "ZZRingElem" => "Q12503",             # integer
  "QQFieldElem" => "Q1244890",          # rational number
  "FiniteField" => "Q603880",           # finite field
  "AbsSimpleNumField" => "Q616608",     # algebraic number field
  "PolyRing" => "Q1455652",             # polynomial ring
  "MPolyRing" => "Q1455652",            # polynomial ring
  "PolyRingElem" => "Q43260",           # polynomial
  "MPolyRingElem" => "Q43260",          # polynomial
  "MPolyIdeal" => "Q44649",             # ideal
  "MatElem" => "Q44337",                # matrix
  "Cone" => "Q2256541",                 # convex cone
  "PolyhedralComplex" => "Q60792057",   # polyhedral complex
  "LinearProgram" => "Q17232527",       # linear programming problem
  "SimplicialComplex" => "Q994399",     # simplicial complex
  "Matroid" => "Q898572",               # matroid
  "Graph{Undirected}" => "Q141488",     # graph
  "Graph{Directed}" => "Q1137726",      # directed graph
  "PhylogeneticTree" => "Q242125",      # phylogenetic tree
  "PermGroup" => "Q1412905",            # permutation group
  "FPGroup" => "Q124579556",            # finitely presented group
  "FinGenAbGroup" => "Q142140",         # finitely generated abelian group
  "NormalToricVariety" => "Q1528306",   # toric variety
)

"""
    wikidata_qid(T::Type)

Return the Wikidata item (e.g. `"Q1455652"`) of the mathematical object
represented by the serialized type `T`, or `nothing` if none is known.
"""
wikidata_qid(::Type{T}) where T = get(wikidata_qids, encode_type(T), nothing)
