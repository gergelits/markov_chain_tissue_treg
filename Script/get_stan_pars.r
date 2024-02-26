# get_stan_pars.r

get_stan_pars <- function( stan_pars_v, ps )
{
  stan_pars_script <- sprintf( "set_model_parameters_%s.r", stan_pars_v )
  source( sprintf( "%s/%s", ps$CODE_PATH, stan_pars_script ), local = TRUE )
  
  stan_pars <- list()
  stan_pars$F_MODEL_VER_EXT <- sprintf( "_%s", stan_pars_v )
  stan_pars$PB_WEEK0I <- WEEK0
  
  stan_pars$HOST_DONOR_RATE_ <- HDR
  stan_pars$USE_HDR2_ <- U_HDR2
  stan_pars$USE_WGHS_ <- U_WGHS
  
  stan_pars$Q12_IS_FIXED <- 0
  stan_pars$Q12_VALUE <- Q12_VALUE
  stan_pars$BRF <- BRF
  stan_pars$MAX_FLOW_TO_N3_ <- MAX_FLOW_TO_N3_FACTOR
  
  stan_pars$EPS <- 0.00595
  stan_pars$NORMALIZE <- TRUE
  stan_pars$STATES_5 <- FALSE
  stan_pars$id <- stan_pars_v
  
  return( stan_pars )
}
