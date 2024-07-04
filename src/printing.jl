# Several interfaces (expressify, iteration, ...) require a single object. Use
# OscarPair as an easy way to pass multiple objects without creating an official
# type for the combinations: polynomial + ordering, old iter + new iter, ...
struct OscarPair{S, T}
  first::S
  second::T
end
