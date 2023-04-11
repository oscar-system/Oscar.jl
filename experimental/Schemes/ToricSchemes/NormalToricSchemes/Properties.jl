#######################################
### 1: Properties of ToricCoveredScheme
#######################################

is_finalized(X::ToricCoveredScheme) = is_finalized(X.ntv)
is_normal(X::ToricCoveredScheme) = is_normal(X.ntv)
is_affine(X::ToricCoveredScheme) = is_affine(X.ntv)
is_projective(X::ToricCoveredScheme) = is_projective(X.ntv)
is_projective_space(X::ToricCoveredScheme) = is_projective_space(X.ntv)
is_smooth(X::ToricCoveredScheme) = is_smooth(X.ntv)
is_complete(X::ToricCoveredScheme) = is_complete(X.ntv)
has_torusfactor(X::ToricCoveredScheme) = has_torusfactor(X.ntv)
is_orbifold(X::ToricCoveredScheme) = is_orbifold(X.ntv)
is_simplicial(X::ToricCoveredScheme) = is_simplicial(X.ntv)
is_gorenstein(X::ToricCoveredScheme) = is_gorenstein(X.ntv)
is_q_gorenstein(X::ToricCoveredScheme) = is_q_gorenstein(X.ntv)
is_fano(X::ToricCoveredScheme) = is_fano(X.ntv)


#######################################
### 2: Properties of ToricSpec
#######################################

is_finalized(X::ToricSpec) = is_finalized(X.antv)
is_normal(X::ToricSpec) = is_normal(X.antv)
is_affine(X::ToricSpec) = is_affine(X.antv)
is_projective(X::ToricSpec) = is_projective(X.antv)
is_projective_space(X::ToricSpec) = is_projective_space(X.antv)
is_smooth(X::ToricSpec) = is_smooth(X.antv)
is_complete(X::ToricSpec) = is_complete(X.antv)
has_torusfactor(X::ToricSpec) = has_torusfactor(X.antv)
is_orbifold(X::ToricSpec) = is_orbifold(X.antv)
is_simplicial(X::ToricSpec) = is_simplicial(X.antv)
is_gorenstein(X::ToricSpec) = is_gorenstein(X.antv)
is_q_gorenstein(X::ToricSpec) = is_q_gorenstein(X.antv)
is_fano(X::ToricSpec) = is_fano(X.antv)
