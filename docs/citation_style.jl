#
# This file is included by docs/make_work.jl to define the custom citation style `oscar_style`.
#
# It is heavily based upon
# https://juliadocs.org/DocumenterCitations.jl/v1.0.0/gallery/#Custom-style:-Citation-key-labels
#

import DocumenterCitations

const oscar_style = :oscar

# The long reference string in the bibliography
function DocumenterCitations.format_bibliography_reference(style::Val{oscar_style}, entry)
  return DocumenterCitations.format_labeled_bibliography_reference(style, entry; article_link_doi_in_title=true)
end

# The label in the bibliography
function DocumenterCitations.format_bibliography_label(::Val{oscar_style}, entry, citations)
  return "[$(entry.id)]"
end

function DocumenterCitations.citation_label(style::Val{oscar_style}, entry, citations; _...)
  return entry.id
end

# The order of entries in the bibliography
DocumenterCitations.bib_sorting(::Val{oscar_style}) = :key  # sort by citation key

# The type of html tag to use for the bibliography
DocumenterCitations.bib_html_list_style(::Val{oscar_style}) = :dl

function DocumenterCitations.format_citation(
  style::Val{oscar_style}, cit, entries, citations
)
  # The only difference compared to `:alpha` is the citation label, which is
  # picked up automatically by redefining `citation_label` above.
  return DocumenterCitations.format_labeled_citation(style, cit, entries, citations)
end
