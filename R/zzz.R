# ==========================================================================
# Biobase package initialization
# ==========================================================================
.onLoad <- function(libname, pkgname) {
   require(methods)
   require(Rsamtools)||stop("Package Rsamtools not available. Please install the package")
   library.dynam( "BICseq", package = "BICseq",lib.loc=libname) 
}

.onUnload <- function( libpath ) {
   library.dynam.unload( "BICseq", libpath )
}
