# Synchronize the bibliography entries for GAP packages in
# docs/oscar_references.bib with the package versions shipped by the GAP.jl
# version in the current environment.
#
# Run from the root directory of the Oscar.jl repository with:
# > julia --project=. etc/update_gap_package_bib.jl
#
# Pass `--check` to only report outdated entries (exit code 1) without
# modifying any files.
#
# Entries are recognized by the field `note = {GAP package}`; the citation key
# is the package name. The author, title, year, month and url fields are
# regenerated from the `PackageInfo.g` record of the installed package.
#
# note: the @main function requires julia 1.11 or newer

using GAP

const oscardir = normpath(@__DIR__, "..")
const bibfile = joinpath(oscardir, "docs", "oscar_references.bib")

const entry_regex = r"^@Misc\{(?<key>[^,\s]+),\n(?<body>(?:  .*\n)*?)\}\n"m
const field_regex = r"^  (?<name>\w+)\s*=\s*\{?(?<value>.*?)\}?,?$"m

const months = ["jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec"]

function parse_fields(body::AbstractString)
  return Dict{String,String}(m[:name] => m[:value] for m in eachmatch(field_regex, body))
end

function package_info(name::AbstractString)
  infos = GAP.Globals.GAPInfo.PackagesInfo
  sym = Symbol(lowercase(name))
  hasproperty(infos, sym) || error("GAP package `$name` is not available in the current GAP installation")
  # the first entry is the one with the highest version
  return getproperty(infos, sym)[1]
end

# accepts dd/mm/yyyy and yyyy-mm-dd, the two date formats allowed in PackageInfo.g
function parse_date(date::AbstractString)
  m = match(r"^(\d{2})/(\d{2})/(\d{4})$", date)
  !isnothing(m) && return parse(Int, m[3]), parse(Int, m[2])
  m = match(r"^(\d{4})-(\d{2})-(\d{2})$", date)
  !isnothing(m) && return parse(Int, m[1]), parse(Int, m[2])
  error("cannot parse date `$date`")
end

function generate_entry(name::AbstractString)
  info = package_info(name)
  persons = [p for p in info.Persons if hasproperty(p, :IsAuthor) && p.IsAuthor]
  isempty(persons) && error("GAP package `$name` has no authors in its `PackageInfo.g`")
  authors = join(("$(String(p.LastName)), $(String(p.FirstNames))" for p in persons), " and ")
  year, month = parse_date(String(info.Date))
  title = "$(String(info.PackageName)), $(String(info.Subtitle)), Version $(String(info.Version))"
  url = String(info.PackageWWWHome)
  return """
    @Misc{$name,
      bibkey        = {$name},
      author        = {$authors},
      title         = {$title},
      note          = {GAP package},
      year          = {$year},
      month         = $(months[month]),
      url           = {$url}
    }
    """
end

function (@main)(args)
  check = "--check" in args
  bib = read(bibfile, String)
  updated = bib
  outdated = 0
  for m in eachmatch(entry_regex, bib)
    fields = parse_fields(m[:body])
    get(fields, "note", "") == "GAP package" || continue
    text = generate_entry(m[:key])
    text == m.match && continue
    outdated += 1
    println("$(m[:key]) is outdated")
    updated = replace(updated, m.match => text)
    new_fields = parse_fields(text)
    for field in ["author", "title", "year", "month", "url"]
      old_value, new_value = get(fields, field, ""), get(new_fields, field, "")
      old_value == new_value && continue
      if length(old_value) + length(new_value) < 80
        println("  $field: $old_value -> $new_value")
      else
        println("  $field: $old_value")
        println("  $(" "^length(field))  -> $new_value")
      end
    end
  end
  if outdated == 0
    println("All GAP package entries in docs/oscar_references.bib are up to date.")
    return 0
  end
  if check
    println("\n$outdated outdated entries. Run without `--check` to update them.")
    return 1
  end
  write(bibfile, updated)
  println("\nUpdated $outdated entries. Now run bibtool to standardize the bibliography:")
  println("  bibtool -r .bibtoolrsc docs/oscar_references.bib -o docs/oscar_references.bib")
  return 0
end
