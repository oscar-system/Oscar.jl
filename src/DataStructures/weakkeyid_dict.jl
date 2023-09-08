# Weak key dict using object-id hashing/equality
# see also https://github.com/JuliaCollections/DataStructures.jl/pull/402
# Based on Julia's WeakKeyIdDict

# Type to wrap a WeakRef to furbish it with objectid comparison and hashing.
#
# Note that various getter and setter functions below all need to explicitly
# use `WeakRefForWeakDict(key)` instead of `key` because the automatism that
# works for `WeakRef` does not work here: for `WeakRef` the hash function is
# simply that of the wrapped object, and comparing a `WeakRef` to a value
# automatically unwraps. But this does not work for `WeakRefForWeakDict`
# because we use a custom hash function based on the `objectid` (this is
# important because it allows efficient use of objects as keys even if
# there is no effective hash function for those objects).
struct WeakRefForWeakDict
    w::WeakRef
    WeakRefForWeakDict(wr::WeakRef) = new(wr)
    WeakRefForWeakDict(@nospecialize(v)) = new(WeakRef(v))
end

Base.:(==)(wr1::WeakRefForWeakDict, wr2::WeakRefForWeakDict) = wr1.w.value===wr2.w.value
Base.hash(wr::WeakRefForWeakDict, h::UInt) = Base.hash_uint(3h - objectid(wr.w.value))

"""
    WeakKeyIdDict([itr])

`WeakKeyIdDict()` constructs a hash table where the keys are weak
references to objects which may be garbage collected even when
referenced in a hash table.

See [`Dict`](@ref) for further help.  Note, unlike [`Dict`](@ref),
`WeakKeyIdDict` does not convert keys on insertion, as this would imply the key
object was unreferenced anywhere before insertion.

See also [`WeakRef`](@ref), [`WeakKeyDict`](@ref).
"""
mutable struct WeakKeyIdDict{K,V} <: AbstractDict{K,V}
    ht::Dict{WeakRefForWeakDict,V}
    lock::ReentrantLock
    finalizer::Function
    dirty::Bool

    # Constructors mirror Dict's
    function WeakKeyIdDict{K,V}() where V where K
        t = new(Dict{WeakRefForWeakDict,V}(), ReentrantLock(), identity, 0)
        t.finalizer = k -> t.dirty = true
        return t
    end
end
function WeakKeyIdDict{K,V}(kv) where V where K
    h = WeakKeyIdDict{K,V}()
    for (k,v) in kv
        h[k] = v
    end
    return h
end
WeakKeyIdDict{K,V}(p::Pair) where V where K = setindex!(WeakKeyIdDict{K,V}(), p.second, p.first)
function WeakKeyIdDict{K,V}(ps::Pair...) where V where K
    h = WeakKeyIdDict{K,V}()
    sizehint!(h, length(ps))
    for p in ps
        h[p.first] = p.second
    end
    return h
end
WeakKeyIdDict() = WeakKeyIdDict{Any,Any}()

WeakKeyIdDict(kv::Tuple{}) = WeakKeyIdDict()
Base.copy(d::WeakKeyIdDict) = WeakKeyIdDict(d)

WeakKeyIdDict(ps::Pair{K,V}...)           where {K,V} = WeakKeyIdDict{K,V}(ps)
WeakKeyIdDict(ps::Pair{K}...)             where {K}   = WeakKeyIdDict{K,Any}(ps)
WeakKeyIdDict(ps::(Pair{K,V} where K)...) where {V}   = WeakKeyIdDict{Any,V}(ps)
WeakKeyIdDict(ps::Pair...)                            = WeakKeyIdDict{Any,Any}(ps)

function WeakKeyIdDict(kv)
    try
        Base.dict_with_eltype((K, V) -> WeakKeyIdDict{K, V}, kv, eltype(kv))
    catch
        if !Base.isiterable(typeof(kv)) || !all(x->isa(x,Union{Tuple,Pair}),kv)
            throw(ArgumentError("WeakKeyIdDict(kv): kv needs to be an iterator of tuples or pairs"))
        else
            rethrow()
        end
    end
end

function _cleanup_locked(h::WeakKeyIdDict)
    if h.dirty
        h.dirty = false
        idx = Base.skip_deleted_floor!(h.ht)
        while idx != 0
            if h.ht.keys[idx].w.value === nothing
                Base._delete!(h.ht, idx)
            end
            idx = Base.skip_deleted(h.ht, idx + 1)
        end
    end
    return h
end

Base.sizehint!(d::WeakKeyIdDict, newsz) = sizehint!(d.ht, newsz)
Base.empty(d::WeakKeyIdDict, ::Type{K}, ::Type{V}) where {K, V} = WeakKeyIdDict{K, V}()

Base.IteratorSize(::Type{<:WeakKeyIdDict}) = Base.SizeUnknown()

Base.islocked(wkh::WeakKeyIdDict) = islocked(wkh.lock)
Base.lock(wkh::WeakKeyIdDict) = lock(wkh.lock)
Base.unlock(wkh::WeakKeyIdDict) = unlock(wkh.lock)
Base.lock(f, wkh::WeakKeyIdDict) = lock(f, wkh.lock)
Base.trylock(f, wkh::WeakKeyIdDict) = trylock(f, wkh.lock)

function Base.setindex!(wkh::WeakKeyIdDict{K}, v, key) where K
    !isa(key, K) && throw(ArgumentError("$key is not a valid key for type $K"))
    # 'nothing' is not valid both because 'finalizer' will reject it,
    # and because we therefore use it as a sentinel value
    key === nothing && throw(ArgumentError("`nothing` is not a valid WeakKeyIdDict key"))
    lock(wkh) do
        _cleanup_locked(wkh)
        k = getkey(wkh.ht, WeakRefForWeakDict(key), nothing)
        if k === nothing
            finalizer(wkh.finalizer, key)
            k = WeakRefForWeakDict(key)
        else
            k.w.value = key
        end
        wkh.ht[k] = v
    end
    return wkh
end
function Base.get!(wkh::WeakKeyIdDict{K}, key, default) where {K}
    v = lock(wkh) do
        k = WeakRefForWeakDict(key)
        if key !== nothing && haskey(wkh.ht, k)
            wkh.ht[k]
        else
            wkh[k] = default
        end
    end
    return v
end
function Base.get!(default::Base.Callable, wkh::WeakKeyIdDict{K}, key) where {K}
    v = lock(wkh) do
        k = WeakRefForWeakDict(key)
        if key !== nothing && haskey(wkh.ht, k)
            wkh.ht[k]
        else
            wkh[k] = default()
        end
    end
    return v
end

function Base.getkey(wkh::WeakKeyIdDict{K}, kk, default) where K
    k = lock(wkh) do
        local k = getkey(wkh.ht, WeakRefForWeakDict(kk), nothing)
        k === nothing && return nothing
        return k.w.value
    end
    return k === nothing ? default : k::K
end

Base.map!(f, iter::Base.ValueIterator{<:WeakKeyIdDict})= map!(f, values(iter.dict.ht))

function Base.get(wkh::WeakKeyIdDict{K}, key, default) where {K}
    key === nothing && throw(KeyError(nothing))
    lock(wkh) do
        return get(wkh.ht, WeakRefForWeakDict(key), default)
    end
end
function Base.get(default::Base.Callable, wkh::WeakKeyIdDict{K}, key) where {K}
    key === nothing && throw(KeyError(nothing))
    lock(wkh) do
        return get(default, wkh.ht, WeakRefForWeakDict(key))
    end
end
function Base.pop!(wkh::WeakKeyIdDict{K}, key) where {K}
    key === nothing && throw(KeyError(nothing))
    lock(wkh) do
        return pop!(wkh.ht, WeakRefForWeakDict(key))
    end
end
function Base.pop!(wkh::WeakKeyIdDict{K}, key, default) where {K}
    key === nothing && return default
    lock(wkh) do
        return pop!(wkh.ht, WeakRefForWeakDict(key), default)
    end
end
function Base.delete!(wkh::WeakKeyIdDict, key)
    key === nothing && return wkh
    lock(wkh) do
        delete!(wkh.ht, WeakRefForWeakDict(key))
    end
    return wkh
end
function Base.empty!(wkh::WeakKeyIdDict)
    lock(wkh) do
        empty!(wkh.ht)
    end
    return wkh
end
function Base.haskey(wkh::WeakKeyIdDict{K}, key) where {K}
    key === nothing && return false
    lock(wkh) do
        return haskey(wkh.ht, WeakRefForWeakDict(key))
    end
end
function Base.getindex(wkh::WeakKeyIdDict{K}, key) where {K}
    key === nothing && throw(KeyError(nothing))
    lock(wkh) do
        return getindex(wkh.ht, WeakRefForWeakDict(key))
    end
end
Base.isempty(wkh::WeakKeyIdDict) = length(wkh) == 0
function Base.length(t::WeakKeyIdDict)
    lock(t) do
        _cleanup_locked(t)
        return length(t.ht)
    end
end

function Base.iterate(t::WeakKeyIdDict{K,V}, state...) where {K, V}
    return lock(t) do
        while true
            y = iterate(t.ht, state...)
            y === nothing && return nothing
            wkv, state = y
            k = wkv[1].w.value
            GC.safepoint() # ensure `k` is now gc-rooted
            k === nothing && continue # indicates `k` is scheduled for deletion
            kv = Pair{K,V}(k::K, wkv[2])
            return (kv, state)
        end
    end
end

Base.filter!(f, d::WeakKeyIdDict) = Base.filter_in_one_pass!(f, d)
