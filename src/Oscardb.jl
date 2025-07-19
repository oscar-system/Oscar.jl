module Oscardb

using ..Oscar
using ..Oscar: load

# for ca certificates
using NetworkOptions

using Mongoc

using URIs

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

Type containing the results of a query.
Can be iterated, but the iterator can not be reset. For this cause, one has to query again.
"""
struct Cursor
   mcursor::Mongoc.Cursor{Mongoc.Collection}
end

"""
      get_db()

Connect to the `OscarDB` and return `Database` instance.

The uri of the server can be set in advance by writing its `String` representation
into ENV["OSCARDB_TEST_URI"].
(used to connect to the github services container for testing)
# Examples
```julia-repl
julia> db = Oscar.Oscardb.get_db();

julia> typeof(db)
Oscar.Oscardb.Database
```
"""
function get_db(;dev=false)
  # we explicitly set the cacert file, otherwise we might get connection errors because the certificate cannot be validated
  username = "oscar-pub"
  password = escapeuri(raw"oShea/fooC$ie6h")
  client = Mongoc.Client(get(ENV, "OSCARDB_TEST_URI", "mongodb://$username:$password@db.oscar-system.org/oscar?directConnection=true&authSource=users&tls=true&appName=mongosh+2.4.2&sslCertificateAuthorityFile=$(NetworkOptions.ca_roots_path())"))
  return Database(client[dev ? OSCAR_DEV_DB : OSCAR_DB])
end

"""
      getindex(db::Database, name::AbstractString)

Return a `Polymake.Polydb.Collection{Polymake.BigObject}` instance
from `db` with the given `name`.
Sections and collections in the name are connected with the '.' sign.
# Examples
```julia-repl
julia> db = Polymake.Polydb.get_db();

julia> collection = getindex(db, "Polytopes.Lattice.SmoothReflexive")
Polymake.Polydb.Collection{Polymake.BigObject}: Polytopes.Lattice.SmoothReflexive

julia> collection = db["Matroids.Small"]
Polymake.Polydb.Collection{Polymake.BigObject}: Matroids.Small
```
"""
Base.getindex(db::Database, name::AbstractString) = Collection(db.mdb[name])

"""
      length(c::Collection{T}, d::Dict=Dict())

Count documents in a collection `c` matching the criteria given by `d`.

# Examples
```julia-repl
julia> db = Polymake.Polydb.get_db();

julia> collection = db["Polytopes.Lattice.SmoothReflexive"];

julia> query = Dict("DIM"=>3, "N_FACETS"=>5);

julia> length(collection, query)
4
```
"""
function Base.length(c::Collection, d::Dict=Dict())
   return Base.length(c.mcol, Mongoc.BSON(d))
end

"""
      find(c::Collection{T}, d::Dict=Dict(); opts::Union{Nothing, Dict})

Search a collection `c` for documents matching the criteria given by `d`.
Apply search options `opts`.
# Examples
```julia-repl
julia> db = Polymake.Polydb.get_db();

julia> collection = db["Polytopes.Lattice.SmoothReflexive"];

julia> query = Dict("DIM"=>3, "N_FACETS"=>5);

julia> results = Polymake.Polydb.find(collection, query);

julia> typeof(results)
Polymake.Polydb.Cursor{Polymake.BigObject}
```
"""
function Mongoc.find(c::Collection, d::Dict=Dict();
                     opts::Union{Nothing, Dict}=nothing)
   return Cursor(Mongoc.find(c.mcol, Mongoc.BSON(d); options=(isnothing(opts) ? nothing : Mongoc.BSON(opts))))
end

"""
      find_one(c::Collection{T}, d::Dict=Dict(); opts::Union{Nothing, Dict})

Return one document from a collection `c` matching the criteria given by `d`.
`T` can be chosen from `Polymake.BigObject` and `Mongoc.BSON`.
Apply search options `opts`.
# Examples
```julia-repl
julia> db = Polymake.Polydb.get_db();

julia> collection = db["Polytopes.Lattice.SmoothReflexive"];

julia> query = Dict("DIM"=>5, "N_FACETS"=>8);

julia> opts = Dict("skip"=>13);

julia> pm_object = Polymake.Polydb.find_one(collection, query, opts=opts);

julia> typeof(pm_object)
Polymake.LibPolymake.BigObjectAllocated

julia> collection_bson = Polymake.Polydb.Collection{Mongoc.BSON}(collection);

julia> pm_object = Polymake.Polydb.find_one(collection_bson, query, opts=opts);

julia> typeof(pm_object)
Mongoc.BSON
```
"""
function find_one(c::Collection, d::Dict=Dict(); opts::Union{Nothing, Dict}=nothing)
   p = Mongoc.find_one(c.mcol, Mongoc.BSON(d); options=(isnothing(opts) ? nothing : Mongoc.BSON(opts)))
   return isnothing(p) ? nothing : parse_document(p)
end

"""
      load(q::String, c::String, s::Symbol = :Oscardb)

Load an entry given by the query `q` in the collection `c` in the database `s`.
Currently, only the Oscar DB (indicated by `:Oscardb`) is supported.
"""
function Oscar.load(q::String, c::String, s::Symbol = :Oscardb)
    if s != :Oscardb
        println("Not implemented for non Oscar databases!")
        @assert s == :Oscardb
    end
    db = get_db()
    collection = getindex(db, c)
    result = Oscardb.find_one(collection, Dict("_id" => q))
    return result
end

"""
      parse_document(bson::Mongoc.BSON)
Create an Oscar object from the data given by `bson`.
!!!!!Note that examples may need to be fixed!!!!!
# Examples
```julia-repl
julia> db = Oscar.Oscardb.get_db();

julia> collection = Oscar.Oscardb.Collection{Mongoc.BSON}(db["Polytopes.Lattice.SmoothReflexive"])
Oscar.Oscardb.Collection{Mongoc.BSON}: Polytopes.Lattice.SmoothReflexive

julia> bson = collect(Oscar.Oscardb.find(collection, "DIM"=>3, "N_FACETS"=>5))[1];

julia> bo = Oscar.Oscardb.parse_document(bson);

julia> typeof(bo) 
Polymake.BigObjectAllocated
```
"""
function parse_document(bson::Mongoc.BSON)
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

   

end
