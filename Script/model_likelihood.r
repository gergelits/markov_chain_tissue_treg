# model_likelihood.r

calculate_model_loglik <- function( m, chain = NA )
{
  chain_dir <- get_path( "chain_dir", m = m, 
                         sel_type = "each", sel_chains = chain )
  Q <- read_Q_matrix( chain_dir = chain_dir )
  af <- read_af_matrix( m = m, chain = chain )
  
  model_vars_dir <- get_path( "model_vars_dir", m = m )
  
  load( file.path( model_vars_dir, "pb.RData" ) )
  load( file.path( model_vars_dir, "model_vars_.RData" ) )

    
  alpha <- list()
  for ( j in 1 : wn_ )
  {
    a = a0i_ %*% expm( 7 * ( w_[ j ] - w0i_ ) * Q )    
    alpha[[ j ]] <- list()
    for ( k in 1 : tn_ ) 
    {
      alpha[[ j ]][[ k ]] = c( 
        a[ ( ( k-1 ) * pn_ + 1 ) : ( k * pn_ ) ],
        a[ ( ( tn_ + k - 1 ) * pn_ + 1 ) : ( ( tn_ + k ) * pn_ ) ] )
      alpha[[ j ]][[ k ]] = alpha[[ j ]][[ k ]] * 
        ( af[ j, k ] / sum( alpha[[ j ]][[ k ]] ) )
    }
  }

  tidyr::expand_grid( wk = wi_ %>% unique,
                      tis = ti_ %>% unique,
                      loglik = 0 ) %>% 
    arrange( wk, tis ) ->
    loglik_parts
  
  for ( i in 1 : n_ )
  {
    # Blood, wk = 0 is not counted
    if ( ! ( wi_[ i ] == 1 && ti_[ i ] == 1 ) )
    {
      lik_i <- ( ddirichlet( x_[ i, ], alpha[[ wi_[ i ] ]][[ ti_[ i ] ]] ) )
      stopifnot( lik_i > 0 )
      loglik_parts$loglik[ loglik_parts$wk == wi_[ i ] & 
                             loglik_parts$tis == ti_[ i ] ] <- 
        loglik_parts$loglik[ loglik_parts$wk == wi_[ i ] & 
                               loglik_parts$tis == ti_[ i ] ] + log( lik_i )
    }
  }
  loglik <- sum( loglik_parts$loglik )
  loglik_parts %>% 
    dplyr::rename( model_part = tis ) %>% 
    group_by( model_part ) %>% 
    summarise( ., loglik = sum( loglik ) ) ->
    loglik_model_parts
  
  res <- list( loglik, loglik_model_parts )
  names( res ) <- c( "loglik", "loglik_model_parts" )
  return( res )
}
