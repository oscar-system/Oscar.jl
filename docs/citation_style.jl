#
# This file is included by docs/make_work.jl to define the custom citation style `oscar_style`.
#
# It is heavily based upon
# https://juliadocs.org/DocumenterCitations.jl/v1.0.0/gallery/#Custom-style:-Citation-key-labels
#

import DocumenterCitations

const oscar_style = :oscar

# This is the only difference to the :alpha style: the label is the citation key instead of a number.
function DocumenterCitations.citation_label(style::Val{oscar_style}, entry, citations; _...)
  return entry.id
end

# The type of html tag to use for the bibliography
DocumenterCitations.bib_html_list_style(::Val{oscar_style}) = :dl

# The order of entries in the bibliography
DocumenterCitations.bib_sorting(::Val{oscar_style}) = :key  # sort by citation key

# The label in the bibliography
function DocumenterCitations.format_bibliography_label(::Val{oscar_style}, entry, citations)
  label = DocumenterCitations.citation_label(Val(oscar_style), entry, citations)
  return "[$label]"
end

# The long reference string in the bibliography
function DocumenterCitations.format_bibliography_reference(style::Val{oscar_style}, entry)
  return DocumenterCitations.format_labeled_bibliography_reference(style, entry; article_link_doi_in_title=true)
end

# The in-text citation format
function DocumenterCitations.format_citation(style::Val{oscar_style}, cit, entries, citations)
  return DocumenterCitations.format_labeled_citation(style, cit, entries, citations)
end
