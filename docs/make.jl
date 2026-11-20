# Build the manual for deployment; run as `julia --project=. docs/make.jl`.
# For a local build use `Oscar.@build_doc`.
using Oscar

Oscar.docs_env()

using OscarDocs

OscarDocs.build(Oscar; warnonly=false, local_build=false, doctest=false, open_browser=false)
