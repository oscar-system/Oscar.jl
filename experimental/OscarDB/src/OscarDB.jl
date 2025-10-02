module OscarDB

using ..Oscar
using ..Oscar.Serialization

import Oscar.Serialization: load_object, save_object, type_params

# Transitivesimplicialcomplex imports
import Oscar:
  simplicial_complex,
  dim,
  f_vector,
  automorphism_group,
  homology,
  betti_numbers,
  n_vertices

# for ca certificates
import NetworkOptions
import Mongoc: Mongoc, find, find_one

const OSCAR_DB = "oscar"
const OSCAR_DEV_DB = "oscar-dev"

"""
    Database

Type for referencing a specific database (usually the `polyDB`)
"""
struct Database
  mdb::Mongoc.Database
end

"""
    Collection

Type for referencing a specific collection.
`T<:Union{Polymake.BigObject, Mongoc.BSON}` defines the template and/or element types
returned by operations applied on objects of this type.
"""
struct Collection
  mcol::Mongoc.Collection
end

"""
    Cursor

Type for representing the results of a query.
Can be iterated, but the iterator can not be reset. For that, one has to query again.
"""
struct Cursor
  mcursor::Mongoc.Cursor{Mongoc.Collection}
end

@doc raw"""
    get_db()

Connect to the `OscarDB` and return `Database` instance.

The uri of the server can be set in advance by writing its `String` representation
into ENV["OSCARDB_TEST_URI"].
(used to connect to the github services container for testing)
# Examples
```julia-repl
julia> db = Oscar.OscarDB.get_db();

julia> typeof(db)
Oscar.OscarDB.Database
```
"""
function get_db(;dev=false)
  # we explicitly set the cacert file, otherwise we might get connection errors because the certificate cannot be validated
  username = "oscar-pub"
  password = "oShea%2FfooC%24ie6h"
  client = Mongoc.Client(get(ENV, "OSCARDB_TEST_URI", "mongodb://$username:$password@db.oscar-system.org/oscar?directConnection=true&authSource=users&tls=true&appName=mongosh+2.4.2&sslCertificateAuthorityFile=$(NetworkOptions.ca_roots_path())"))
  return Database(client[dev ? OSCAR_DEV_DB : OSCAR_DB])
end

#TODO add examples
"""
    getindex(db::Database, name::AbstractString)

Return a `Oscar.OscarDB.Collection` instance
from `db` with the given `name`.
Sections and collections in the name are connected with the '.' sign.
# Examples
"""
Base.getindex(db::Database, name::AbstractString) = Collection(db.mdb[name])

"""
    length(c::Collection{T}, d::Dict=Dict())

Count documents in a collection `c` matching the criteria given by `d`.

# Examples
Same as above, but faster.
```julia-repl
julia> db = Oscar.OscarDB.get_db();

julia> tscit = Oscar.OscarDB.find(db["TransitiveSimplicialComplexes"], Dict("data.betti_numbers" => ["0", "0", "0", "1"]));

julia> length(tscit)
63
```
"""
function Base.length(c::Collection, d::Dict=Dict())
  return Base.length(c.mcol, Mongoc.BSON(d))
end

@doc raw"""
    find(c::Collection{T}, d::Dict=Dict(); opts::Union{Nothing, Dict})

Returns an iterator over documents in the collection `c` which match the criteria given by `d`.
Apply search options `opts`.

# Examples

Here we show how to find all vertex-transitive 3-dimensional rational homology spheres in that database (i.e., with up to 15 vertices). 

The dimension $d$ of the manifold is here implicitly given as the length of the vector of Betti numbers less one.

```julia-repl
julia> db = Oscar.OscarDB.get_db();

julia> tscit = Oscar.OscarDB.find(db["TransitiveSimplicialComplexes"], Dict("data.betti_numbers" => ["0", "0", "0", "1"]));

julia> length([tsc for tsc in tscit])
63
```
"""
function Mongoc.find(c::Collection, d::Dict=Dict();
                     opts::Union{Nothing, Dict}=nothing)
  return Cursor(Mongoc.find(c.mcol, Mongoc.BSON(d); options=(isnothing(opts) ? nothing : Mongoc.BSON(opts))))
end

@doc raw"""
    find_one(c::Collection{T}, d::Dict=Dict(); opts::Union{Nothing, Dict})

Return one document from a collection `c` matching the criteria given by `d`.
`T` can be chosen from `Polymake.BigObject` and `Mongoc.BSON`.
Apply search options `opts`.

# Examples

Here we show how to find one vertex-transitive combinatorial surface with reduced rational Betti numbers $\tilde\beta_0 = \beta_1 = \beta_2 = 0$.

The dimension of the manifold is here implicitly given as the length of the vector of Betti numbers less one.

```jldoctest
julia> db = Oscar.OscarDB.get_db();

julia> tsc = Oscar.OscarDB.find_one(db["TransitiveSimplicialComplexes"], Dict("data.betti_numbers" => ["0", "0", "0"]));
```
In this case we find the unique six-vertex triangulation of the real projective plane.
```jldoctest
julia> facets(tsc.complex)
10-element Vector{Set{Int64}}:
 Set([2, 3, 1])
 Set([6, 2, 1])
 Set([5, 3, 1])
 Set([5, 4, 1])
 Set([4, 6, 1])
 Set([4, 2, 3])
 Set([5, 4, 2])
 Set([5, 6, 2])
 Set([4, 6, 3])
 Set([5, 6, 3])
```
"""
function find_one(c::Collection, d::Dict=Dict(); opts::Union{Nothing, Dict}=nothing)
  p = Mongoc.find_one(c.mcol, Mongoc.BSON(d); options=(isnothing(opts) ? nothing : Mongoc.BSON(opts)))
  return isnothing(p) ? nothing : parse_document(p)
end

# not sure we need this?
# I am guessing this is to handle when there are many different types of
# databases.
"""
    load(q::String, c::String, s::Symbol = :OscarDB)

Load an entry given by the query `q` in the collection `c` in the database `s`.
Currently, only the Oscar DB (indicated by `:OscarDB`) is supported.
"""
# function Oscar.load(q::String, c::String, s::Symbol = :OscarDB)
#   if s != :OscarDB
#     println("Not implemented for non Oscar databases!")
#     @assert s == :OscarDB
#   end
#   db = get_db()
#   collection = getindex(db, c)
#   result = OscarDB.find_one(collection, Dict("_id" => q))
#   return result
# end

#TODO clean the docs of this function up
"""
    parse_document(bson::Mongoc.BSON)

Create an Oscar object from the data given by `bson`.
"""
function parse_document(bson::Mongoc.BSON)
  # TODO should accept override p
  str = Mongoc.as_json(bson)
  return Oscar.load(IOBuffer(str))
end

# Iterator

Base.IteratorSize(::Type{<:Cursor}) = Base.SizeUnknown()
Base.IteratorSize(::Type{<:Collection}) = Base.SizeUnknown()

# functions for `BSON` iteration
Base.iterate(cursor::Cursor, state::Nothing=nothing) =
  iterate(cursor.mcursor, state)

Base.iterate(coll::Collection) =
  iterate(coll.mcol)

Base.iterate(coll::Collection, state::Mongoc.Cursor) =
  iterate(coll.mcol, state)

# shows information about a specific Collection
function Base.show(io::IO, coll::Collection)
  db = Database(coll.mcol.database)
  print(io, typeof(coll), ": ", coll.mcol.name)
end

include("exports.jl")

include("LeechPairs.jl")
include("TransitiveSimplicialComplex.jl")
include("serialization.jl")


end # Module

using .OscarDB

include("exports.jl")

