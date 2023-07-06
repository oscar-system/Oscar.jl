#####################################################
# 1. Properties for literature models
#####################################################

has_arxiv_id(t::AbstractFTheoryModel) = has_attribute(t, :arxiv_id)
has_description(t::AbstractFTheoryModel) = has_attribute(t, :description)
has_doi(t::AbstractFTheoryModel) = has_attribute(t, :doi)
has_equation_number(t::AbstractFTheoryModel) = has_attribute(t, :equation_number)
has_link(t::AbstractFTheoryModel) = has_attribute(t, :link)
has_resolutions(t::AbstractFTheoryModel) = has_attribute(t, :resolutions)
has_version(t::AbstractFTheoryModel) = has_attribute(t, :version)