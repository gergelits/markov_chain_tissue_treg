# estimate_one_markov_model.r

estimate_one_markov_model <- function( 
    m = NULL,
    use_analysis_done_as_input = FALSE,
    sel_models_file_name )
{
  ps <- m$ps
  model_name <- m$model_name
  
  # fit markov model
  
  # NOTE: possible extension with ps$ANALYSIS_DONE_PATH
  model_dir <- get_path( type = "model_dir", m = m ) 
  if ( !dir.exists( model_dir ) ) { dir.create( model_dir, recursive = TRUE ) }
  
  # LOAD_PREMODEL_VARS_AND_DATA_FROM_markov_premodel_prints
    model_vars_dir <- get_path( type = "model_vars_dir", m = m )  
    load( file = sprintf( "%s/model_vars_.RData", model_vars_dir ) )
    load( file = sprintf( "%s/pb.RData", model_vars_dir ) )
    load( file = sprintf( "%s/other_vars.RData", model_vars_dir ) )
    
    # decide whether to run new MCMC simulations  
    MCMC_is_done <- file.exists( sprintf( 
      "%s/parabio_fit%i.rda", model_dir, m$mcmc_pars$mcmc.iter.n ) )                         
    if ( MCMC_is_done ) {
      load( sprintf( "%s/parabio_fit%i.rda", model_dir, m$mcmc_pars$mcmc.iter.n ) )
    } else {
      # run new MCMC simulations
      
      # copy model first to avoid access from multiple sites ( e.g., on cluster )
      file.copy( from = sprintf( "%s/%s.stan", ps$CODE_PATH, model_name ),
                 to = sprintf( "%s/%s.stan", model_dir, model_name ) )
      cat( "\nTranslating from Stan to C++ and compiling in C++ ...\n" )
      
      # translate from Stan to C++ and compile in C++
      parabio.model <- stan_model( 
        sprintf( "%s/%s.stan", model_dir, model_name ), verbose = FALSE )        
      
      # estimate the model by new MCMC simulation:
      cat( "\n############\n" )
      cat( "New MCMC simulation has started.\n" )
      cat( "############\n\n" )
      data_list <- list( host_donor_rate_ = host_donor_rate_,
                         use_wghs_ = use_wghs_, 
                         use_hdr2_ = use_hdr2_, hdr2_lower_ = hdr2_lower_,
                         max_flow_to_N3_ = max_flow_to_N3_,
                         q12_ = q12_, q12_fixed_ = q12_fixed_, 
                         b_rate_rel_to_blood_naive_ = b_rate_rel_to_blood_naive_,                   
                         n_ = n_, wn_ = wn_, pn_ = pn_, tn_ = tn_,
                         w0i_ = w0i_, w_ = w_, wi_ = wi_, ti_ = ti_, x_ = x_,    
                         al_ = al_, a0i_ = a0i_ )
      model_sim <- purrr::quietly( sampling )(                                   
        diagnostic_file = sprintf( "%s/diagnostic_file.csv", model_dir ),
        sample_file = sprintf( "%s/sample_file.csv", model_dir ),
        verbose = TRUE,
        object = parabio.model, data = data_list, 
        iter = m$mcmc_pars$mcmc.iter.n, 
        warmup = m$mcmc_pars$mcmc.warmup.n, 
        chains = m$mcmc_pars$mcmc.chain.n,        
        control = list( adapt_delta = m$mcmc_pars$sampling.adapt_delta,
                        max_treedepth = 10 ) )                                    
      # MCMC simulation finished
      
      parabio.fit <- model_sim$result
      save( parabio.fit, file = sprintf(
        "%s/parabio_fit%i.rda", model_dir, m$mcmc_pars$mcmc.iter.n ) )
      # MCMC simulation saved
      
      # save warnings, output, and messages
      sink( sprintf( "%s/parabio_fit_sampling_and_warnings.txt", model_dir ) )
      cat( "Warning messages:\n", model_sim$warnings, "\n" )
      cat( "Output:\n",           model_sim$output,   "\n" )
      cat( "Messages:\n",         model_sim$messages, "\n" )
      sink()
    }
    if ( PRINT_PARABIO_FIT_ALL_CHAINS <- TRUE ) {
      sink( file = sprintf( "%s/parabio_fit_print.txt", model_dir ) )
      width.backup <- getOption( "width" ); options( width = 200 ) 
      print( parabio.fit )
      options( width = width.backup )
      sink()
    }
    
    
    # Which chains to get results for:
    sel_chains_set_max <- select_mcmc_chains( 
      parabio.fit = parabio.fit, 
      mcmc.chain.n = m$mcmc_pars$mcmc.chain.n, 
      type = "max" ) %>%
      paste( collapse = "_" )
    
    sel_chains_table <- tibble( 
      sel_chains_set = c( 1 : m$mcmc_pars$mcmc.chain.n, sel_chains_set_max ),
      sel_chains_type = c( rep( "each", m$mcmc_pars$mcmc.chain.n ), "max" ) )
    
    sel_chains_set_3mad <- select_mcmc_chains(
      parabio.fit = parabio.fit, mcmc.chain.n = m$mcmc_pars$mcmc.chain.n, 
      type = "3mad" ) %>%
      paste( collapse = "_" )
    if ( sel_chains_set_3mad != "" ) 
    {
      sel_chains_table <- bind_rows( sel_chains_table, 
                                     tibble( sel_chains_set = sel_chains_set_3mad, 
                                             sel_chains_type = "3mad" ) )
    }
    
    write_csv( sel_chains_table, sprintf( "%s/sel_chains_table.csv", model_dir ) )
    
    # either all chain, which takes ~5 minutes:
    for( i in 1 : nrow( sel_chains_table ) )  {
      
      # or only max lp__ chain to save time:
      # for( i in ( which( sel_chains_table$sel_chains_type == "max" ) ) )
      {
    plot_par_densities_and_calc_Q(
      m = m,
      parabio.fit = parabio.fit,
      sel_chains = sel_chains_table$sel_chains_set[ i ],
      sel_type = sel_chains_table$sel_chains_type[ i ] ) ->   
      Q  
    
    test_eqeq(
      m = m,
      sel_chains = sel_chains_table$sel_chains_set[ i ],
      sel_type = sel_chains_table$sel_chains_type[ i ] ) -> 
      final_fit_eqeq
    
    # if debugging needed write where the equilibrium is not met
    if ( !final_fit_eqeq & sel_chains_table$sel_chains_set[ i ] == "1" & 
         sel_chains_table$sel_chains_type[ i ] == "each" ) {
      for ( iter in sample( 1 : 100, 10, replace = FALSE ) )
      {
        test_eqeq(
          m = m,
          parabio.fit = parabio.fit, 
          sel_chains = sel_chains_table$sel_chains_set[ i ],
          sel_type = sel_chains_table$sel_chains_type[ i ],
          iter = iter )           
      }
    }
      
    get_Qcounts(
      m = m,
      sel_chains = sel_chains_table$sel_chains_set[ i ],
      sel_type = sel_chains_table$sel_chains_type[ i ] )
    
    get_and_plot_trajectories(
      m = m,
      sel_chains = sel_chains_table$sel_chains_set[ i ],
      sel_type = sel_chains_table$sel_chains_type[ i ],
      sel_models_file_name = sel_models_file_name )
  }
  }
}
