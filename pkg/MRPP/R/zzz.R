.First.lib <- function(lib, pkg) library.dynam("mrpp", pkg, lib)

.Last.lib <- function(libpath) library.dynam.unload("mrpp", libpath)
