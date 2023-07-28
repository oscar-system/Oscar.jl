#
# This file is included by docs/make_work.jl to define the custom citation style `oscar_style`.
#
# It is heavily based upon
# https://juliadocs.org/DocumenterCitations.jl/v1.0.0/gallery/#Custom-style:-Citation-key-labels
#

import DocumenterCitations

# we use some (undocumented) internal helper functions for formatting...
using DocumenterCitations: format_names, tex2unicode, italicize_md_et_al

const oscar_style = :oscar

# The long reference string in the bibliography
function DocumenterCitations.format_bibliography_reference(::Val{oscar_style}, entry)
  return DocumenterCitations.format_bibliography_reference(:numeric, entry)
end

# The label in the bibliography
function DocumenterCitations.format_bibliography_label(::Val{oscar_style}, entry, citations)
  return "[$(entry.id)]"
end

# The order of entries in the bibliography
DocumenterCitations.bib_sorting(::Val{oscar_style}) = :nyt  # name, year, title

# The type of html tag to use for the bibliography
DocumenterCitations.bib_html_list_style(::Val{oscar_style}) = :dl

function DocumenterCitations.format_citation(
  ::Val{oscar_style}, entry, citations; note, cite_cmd, capitalize, starred
)
  link_text = isnothing(note) ? "[$(entry.id)]" : "[$(entry.id), $note]"
  if cite_cmd == :citet
    et_al = starred ? 0 : 1  # 0: no "et al."; 1: "et al." after 1st author
    names = tex2unicode(
      format_names(entry; names=:lastonly, and=true, et_al, et_al_text="*et al.*")
    )
    capitalize && (names = uppercasefirst(names))
    link_text = italicize_md_et_al("$names $link_text")
  end
  return link_text
end
