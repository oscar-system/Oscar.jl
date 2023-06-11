function _eval_poly(E::Expr, vars)
  @assert E.head == :call
  if E.args[1] == :+
    return reduce(+, (_eval_poly(E.args[i], vars) for i in 2:length(E.args)))
  elseif E.args[1] == :*
    return reduce(*, (_eval_poly(E.args[i], vars) for i in 2:length(E.args)))
  elseif E.args[1] == :-
    if length(E.args) == 2
      return -_eval_poly(E.args[2], vars)
    else
      @assert length(E.args) == 3
      return _eval_poly(E.args[2], vars) - _eval_poly(E.args[3], vars)
    end
  elseif E.args[1] == :^
    return _eval_poly(E.args[2], vars)^_eval_poly(E.args[3], vars)
  elseif E.args[1] == ://
    @assert E.args[2] isa Number && E.args[3] isa Number
    return E.args[2]//E.args[3]
  end
end

function _eval_poly(E::Symbol, vars)
  return vars[E]
end

function _eval_poly(E::Number, vars)
  return E
end

function eval_poly(s::String, R)
  symR = symbols(R) # Symbol[]
  genR = gens(R)
  return R(_eval_poly(Meta.parse(s), Dict(symR[i] => genR[i] for i in 1:length(symR))))
end

# Example
# julia> Qx, (x1, x2) = QQ["x1", "x2"];
#
# julia> eval_poly("-x1 - 3//5*x2^3 + 5 - 3", Qx)
# -x1 - 3//5*x2^3 + 2
