if !isdefined(AbstractAlgebra, Symbol("@enable_all_show_via_expressify"))

# Only when AbstractAlgebra.expressify(a::T; context = nothing) has been
# defined may these be used.
# AA defines Base.show for "text/latex" and "text/html" for a general set
# of x, but for backward compatibility it is not defined for general x and
# "text/plain" or the mime-less version.
# Rationale: when neither Base.show nor AA.expressify is defined for T, then,
# since expressify calls Base.show for backward compatibility, a definition of
# Base.show in terms of expressify would give a stack overflow.

macro enable_all_show_via_expressify(T)
  return quote
    function Base.show(io::IO, x::$(esc(T)))
       AbstractAlgebra.show_via_expressify(io, x)
    end

    function Base.show(io::IO, mi::MIME"text/plain", x::$(esc(T)))
       AbstractAlgebra.show_via_expressify(io, mi, x)
    end

    function Base.show(io::IO, mi::MIME"text/latex", x::$(esc(T)))
       AbstractAlgebra.show_via_expressify(io, mi, x)
    end

    function Base.show(io::IO, mi::MIME"text/html", x::$(esc(T)))
       AbstractAlgebra.show_via_expressify(io, mi, x)
    end
  end
end
end

# Several interfaces (expressify, iteration, ...) require a single object. Use
# OscarPair as an easy way to pass multiple objects without creating an official
# type for the combinations: polynomial + ordering, old iter + new iter, ...
struct OscarPair{S, T}
  first::S
  second::T
end

# A non-compiletime preference
function allow_unicode(flag::Bool)
  old_flag = is_unicode_allowed()
  @set_preferences!("unicode" => flag)
  return old_flag
end
export allow_unicode

function is_unicode_allowed()
  return @load_preference("unicode", default = false)
end

function with_unicode(f::Function)
  old_allow_unicode = allow_unicode(true);
  f()
  allow_unicode(old_allow_unicode);
end
