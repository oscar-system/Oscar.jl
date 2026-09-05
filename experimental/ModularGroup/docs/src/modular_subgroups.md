# Modular subgroups

## Construction of modular subgroups

```@docs
modular_subgroup_via_right_action(s::PermGroupElem, t::PermGroupElem; check=true)
modular_subgroup_via_left_action(s::PermGroupElem, t::PermGroupElem; check=true)
s_right_perm(G::ModularGroup)
t_right_perm(G::ModularGroup)
r_right_perm(G::ModularGroup)
j_right_perm(G::ModularGroup)
```

## Generators and membership

```@docs
word_gens(G::ModularGroup)
s_t_decomposition(M::ZZMatrix)
is_word_element_of(w::FPGroupElem, G::ModularGroup)
coset_right_action_of(A::ZZMatrix, G::ModularGroup)
coset_left_action_of(A::ZZMatrix, G::ModularGroup)
```
