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
  if entry.id in ["OEIS", "Stacks"]
    # patch DocumenterCitations.format_bibliography_reference(:numeric, entry) for corporate authors (https://github.com/JuliaDocs/DocumenterCitations.jl/issues/44)
    authors = tex2unicode(format_names(entry; names=:full))
    title = DocumenterCitations.xtitle(entry)
    if !isempty(title)
      title = "<i>" * tex2unicode(title) * "</i>"
    end
    linked_title = DocumenterCitations.linkify(title, entry.access.url)
    published_in = DocumenterCitations.linkify(
      tex2unicode(DocumenterCitations.format_published_in(entry)),
      DocumenterCitations._doi_link(entry),
    )
    eprint = DocumenterCitations.format_eprint(entry)
    note = DocumenterCitations.format_note(entry)
    parts = String[]
    for part in (authors, linked_title, published_in, eprint, note)
      if !isempty(part)
        push!(parts, part)
      end
    end
    html = DocumenterCitations._join_bib_parts(parts)
    return html
  end
  return DocumenterCitations.format_bibliography_reference(:numeric, entry)
end

# The label in the bibliography
function DocumenterCitations.format_bibliography_label(::Val{oscar_style}, entry, citations)
  return "[$(entry.id)]"
end

# The order of entries in the bibliography
DocumenterCitations.bib_sorting(::Val{oscar_style}) = :key  # sort by citation key

# The type of html tag to use for the bibliography
DocumenterCitations.bib_html_list_style(::Val{oscar_style}) = :dl

function DocumenterCitations.format_citation(
  ::Val{oscar_style}, entry, citations; note, cite_cmd, capitalize, starred
)
  link_text = isnothing(note) ? "[$(entry.id)]" : "[$(entry.id), $note]"
  if cite_cmd == :citet
    et_al = starred ? 0 : 1  # 0: no "et al."; 1: "et al." after 1st author
    if entry.id in ["OEIS", "Stacks"]
      # patch for corporate authors (https://github.com/JuliaDocs/DocumenterCitations.jl/issues/44)
      names = tex2unicode(format_names(entry; names=:full))
    else
      names = tex2unicode(
        format_names(entry; names=:lastonly, and=true, et_al, et_al_text="*et al.*")
      )
    end
    capitalize && (names = uppercasefirst(names))
    link_text = italicize_md_et_al("$names $link_text")
  end
  return link_text
end
