# load_libraries.r

load_libraries <- function()
{
  suppressPackageStartupMessages( {
    library( digest )
    library( expm )
    library( gridExtra )   
    library( gtools )
    library( modeest )
    library( rstan )
    library( stringr )     
    library( tidyverse )
  } )
  # can be used depending on the tidyverse / readr library version loaded
  try( options( readr.show_col_types = FALSE ) )  
  
  # the package posterior needed but better to avoid loading its functions 
  # and masking the default ones ( e.g., mad() )
  # similar for ggpubr
  stopifnot( nchar( system.file( package = "posterior" ) ) > 0 ) 
  stopifnot( nchar( system.file( package = "ggpubr" ) ) > 0 )
}
