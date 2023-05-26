#####################################################
# 1. Properties for general global Tate models
#####################################################

@doc raw"""
    base_fully_specified(t::GlobalTateModel)

Return `true` is the Tate model has a concrete base space and `false` otherwise.

```jldoctest
julia> t = literature_tate_model(arxiv_id = "1109.3454", equ_nr = "3.5")
Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arxiv paper 1109.3454 (equ. 3.5)

julia> base_fully_specified(t)
false
```
"""
base_fully_specified(t::GlobalTateModel) = get_attribute(t, :base_fully_specified)



#####################################################
# 2. Properties for literature models
#####################################################

has_doi(t::GlobalTateModel) = has_attribute(t, :doi)
has_arxiv_id(t::GlobalTateModel) = has_attribute(t, :arxiv_id)
has_version(t::GlobalTateModel) = has_attribute(t, :version)
has_equ_nr(t::GlobalTateModel) = has_attribute(t, :equ_nr)
has_description(t::GlobalTateModel) = has_attribute(t, :description)
has_link(t::GlobalTateModel) = has_attribute(t, :link)
