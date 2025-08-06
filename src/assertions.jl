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

@doc raw"""
    @check(cond)
    @check(cond, msg)

This macro can be used to run expensive internal checks for sanity of 
the input to a function. It needs a variable `check` in the local scope 
and in case `check==true` the condition `cond` gets evaluated. 
Whenever `cond` evaluates to anything else but `false` the test is 
accepted; otherwise an error with the message `msg` is thrown. 

Usually such checks are disabled for internal calls to functions with 
verified input by setting an optional argument `check::Bool=true` to 
`false`. However, this is often forgotten when programming new stuff. 
To snoop for such leaks, one can set the global variable 

    THROW_ERROR_FOR_INTERNAL_CHECKS

in `src/assertions.jl` to `true` and run tests which are not supposed 
to trigger any checks by the `@check` macro. If any such unexpected 
test is nevertheless run, an error is thrown and one can disable the 
faulty internal test following the subsequent stacktrace. 

# Examples
```
  function my_fun(a::Int; check::Bool=true)
    @check a==5 "input is invalid"
    return a
  end;
  
  function my_fun2(s::String; check::Bool=true)
    b = my_fun(length(s)) # Forgotten check=check in the call of `my_fun` here
    return b
  end;
  
  my_fun2("hi!"); # Will throw an error message
  
  my_fun(4, check=false); # Will not throw an error message
  
  my_fun(5); # Will not throw an error message

  my_fun(4); # Will throw an error message
  
  my_fun2("hello", check=false); # Will process fine in general; but throws an error when THROW_ERROR_FOR_INTERNAL_CHECKS == true
```
"""
macro check(cond, msg)
  quote
    if $(esc(:(check)))
      THROW_ERROR_FOR_INTERNAL_CHECKS && error("internal checks are enabled")
      $(esc(cond)) != false || error($(esc(msg)))
    end
  end
end

