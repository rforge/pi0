.First.lib <- function(lib, pkg) library.dynam("MRPP", pkg, lib)

.Last.lib <- function(libpath) library.dynam.unload("MRPP", libpath)
