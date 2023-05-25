########################################################################
# Internal sanity checks for input
#
# Many functions in Oscar make use of an optional argument 
# `check::Bool=true`. If enabled, the function will run expensive 
# sanity checks on the input using the `@check` macro.
#
# For internal calls to such functions with input which has already 
# been verified, one should then set this argument to `false`. However, 
# this is often forgotten when programming. In order to systematically 
# snoop for such leaks, the following global variable can be set to `true` 
# and then tests can be run which are supposed to *not* trigger any such 
# expensive input checks. If an internal test is nevertheless run with 
# @check, then an error is thrown and the faulty call can be disabled 
# following the subsequent stacktrace. 
########################################################################

global const THROW_ERROR_FOR_INTERNAL_CHECKS = false

macro check()
  quote
    if $(esc(:(check)))
      THROW_ERROR_FOR_INTERNAL_CHECKS && error("internal checks are enabled")
    end
  end
end

macro check(cond)
  quote
    if $(esc(:(check)))
      THROW_ERROR_FOR_INTERNAL_CHECKS && error("internal checks are enabled")
      $(esc(cond)) != false || error("internal tests failed")
    end
  end
end

macro check(cond, msg)
  quote
    if $(esc(:(check)))
      THROW_ERROR_FOR_INTERNAL_CHECKS && error("internal checks are enabled")
      $(esc(cond)) != false || error($(esc(msg)))
    end
  end
end

