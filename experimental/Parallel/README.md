# Parallel, resp. distributed computations in Oscar

## Provide simple tools to deploy embarrassingly parallel tasks to workers

The standard julia methods for parallelization do not naturally go well together with 
Oscar datatypes. The main reason is that the `parent` structure needs to be maintained 
when communicating between different processes. In this experimental package we aim 
to providing support for these things in a hopefully convenient way for the user and 
the programmer.

## Status

This is currently work in progress with no guarantee for long term support. 
However, if you would like to use it or should you find out that the methods 
provided here do not fully serve your needs, feel free to reach out to 
`@antonydellavecchia` and `@HechtiDerLachs`. 
