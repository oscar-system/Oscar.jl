################################################################################
## basic functionality for determinantal germs, which needs to be overwritten
## for either correctness or performance
################################################################################


function codim(X::DeterminantalGerm) 
  m, n, t  = determinantal_type(X)
  return _expected_codim(m, n, t, _mat_type(X))
end


dim(X::DeterminantalGerm) = dim(ambient_germ(X)) - codim(X)


@doc raw"""
    ==(X::DeterminantalGerm, Y::DeterminantalGerm)

Return whether the determinantal germs `X` and `Y` are equal, i.e. they have the same determinantal structure and are given by the same defining_matrix.
# Examples:
```jldoctest
julia> R, (x,y,z) = QQ[:x, :y, :z];

julia> A = R[x 0 z;  0 y z]
[x   0   z]
[0   y   z]

julia> X_A = DeterminantalGerm(A, 2, [0,0,0]);


julia> B = R[0 y z;  x 0 z]
[0   y   z]
[x   0   z]

julia> X_B = DeterminantalGerm(B, 2, [0,0,0]);

julia> SpaceGerm(representative(X_A), point(X_A)) == SpaceGerm(representative(X_B), point(X_B))
true

julia> X_A == X_B
false
```
"""
function ==(X::DeterminantalGerm, Y::DeterminantalGerm)
  X === Y && return true
  _mat_type(X) == _mat_type(Y) || return false
  determinantal_type(X) == determinantal_type(Y) || return false
  ambient_coordinate_ring(X) === ambient_coordinate_ring(Y) || return false
  inverted_set(OO(X)) == inverted_set(OO(Y)) || return false
  return fraction.(defining_matrix(X)) == fraction.(defining_matrix(Y))
end
