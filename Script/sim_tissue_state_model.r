# sim_tissue_state_model.r

sim_combined_dwell <- function( m, n_sim = 1000, incl_past = TRUE )
{
  chain_dir <- get_path( "chain_dir", m = m, sel_type = "max",
                         sel_chains = get_max_chain( m ) )
  fig_fname_start <- get_fig_fname_start( m = m )
  figs_dwell_dir <- get_path( "figs_m_dwell", m = m )
  sim_file <- sprintf( "%s/%s_%s_%isims_steps_%s.csv", 
                       figs_dwell_dir, fig_fname_start, 
                       tolower( m$tissue ), n_sim,
                       ifelse( incl_past, "past_fut", "fut" ) )
  sim_done <- file.exists( sim_file )
  Q <- read_Q_matrix( chain_dir )
  if ( is.null( Q ) ) { return( NULL ) }
  
  if ( ! sim_done ) {
    al_ <- get_al_( m = m )
    source( sprintf( "%s/%s", m$ps$CODE_PATH, "tissue_states_model.r" ), 
            local = TRUE )
    m456_ex <- get_tissue_states_model_exit( Q = Q, al_ = al_ )
    m456_en <- get_tissue_states_model_entry( Q = Q, al_ = al_ )
    sim_f <- sim_tissue_states_dwell_time_exit( n_sim = n_sim, m456 = m456_ex )
    
    if ( incl_past ) {
      sim <- sim_tissue_states_past_future( m456_en = m456_en, sim_future = sim_f )
    } else { sim <- sim_f }
    sim %>% write_csv( sim_file )
  } else {
    sim <- read_csv( sim_file )
  }
  
  return( sim )
}


plot_combined_dwell <- function( m, sim, ylim_max = NA, export = TRUE )
{
  # assign order to each simulation (order by total time)
  sim %>%
    group_by( i_sim ) %>% 
    summarise( total_t = sum( time ) ) %>% 
    mutate( i_sim_order = rank( -total_t ) ) -> 
    sim_steps_ordered
  
  # output:
  # paths
  fig_fname_start <- get_fig_fname_start( m = m )
  figs_dwell_dir <- get_path( "figs_m_dwell", m = m )
  
  # table / text
  prop_exit <- sim %>% 
    filter( occur == "exit" ) %>% 
    group_by( state ) %>% 
    summarise( n = n() ) %>% 
    mutate( freq = n / sum( n ) )
  
  # Check to have all "state" in table even if it equals 0.
  prop_exit %>% 
    bind_rows( tibble( 
      state = c( "blood_s5", "blood_s6", "die_s5", "die_s6", "exit_s4" ),
      n = rep( 0, 5 ),
      freq = rep( 0, 5 ) ) ) %>% 
    group_by( state ) %>% 
    filter( n == max( n ) ) %>% unique() %>% 
    ungroup() %>% arrange( state ) -> 
    prop_exit
  
  n_sim <- max( sim$i_sim )
  if ( export ) {
    tibble( tissue = tolower( m$tissue ),
            median = median( sim_steps_ordered$total_t ) ) %>% 
      write_csv( sprintf( "%s/%s_%s_%isims_median_dwell.csv", figs_dwell_dir, 
                          fig_fname_start, tolower( m$tissue ), n_sim ) )
    
    prop_exit %>% 
      write_csv( sprintf( "%s/%s_%s_%isims_prop_exit.csv", figs_dwell_dir, 
                          fig_fname_start, tolower( m$tissue ), n_sim ) )
  }
    
  # figures
  cm_to_in <- 1 / 2.54   
  f <- cm_to_in * 1.7
  pal <- hue_pal()( 3 )
  pal_vector <- pal[ c( rep( c( 2, 3 ), 100 ), 1 ) ]
  names( pal_vector ) <- sprintf( 
    "o%s", stringr::str_pad( c( 1 : ( length( pal_vector ) - 1 ), 0 ), 
                             3, pad = "0" ) )

  prop_die <- prop_exit %>% filter( state %in% c( "die_s5", "die_s6" ) ) %>% 
    dplyr::select( freq ) %>% sum()
  text_die <- sprintf( "die %s%% cells", prop_die * 100 )
  sim_steps_ordered <- 
    sim_steps_ordered %>% filter( !is.na( total_t ) )
  q95 <- round( quantile( sim_steps_ordered$total_t, 0.95 ) )
  
  # Plot only subset of 1000 as readable. While estimated from more.
  sim_sub <- sim %>% filter( i_sim <= 1000 )
  sim_steps_ordered_sub <- sim_sub %>%
    group_by( i_sim ) %>% 
    summarise( total_t = sum( time ) ) %>% 
    mutate( i_sim_order = rank( -total_t ) )
  n_sim <- max( sim_sub$i_sim )
  
  sim_sub %>%
    left_join( ., sim_steps_ordered_sub %>% dplyr::select( i_sim, i_sim_order ),
               by = c( "i_sim" = "i_sim" ) ) %>%
    mutate( occur = factor( 
      occur, levels = sort( unique( occur ), decreasing = TRUE ) ) ) %>% 

    ggplot( aes( x = i_sim_order, y = time, fill = occur ) ) +
    geom_bar( position = "stack", stat = "identity" ) + 
    labs( title = m$tissue, 
          x = "Sampled cell index", 
          y = "Dwell time (days)",
          fill = "Cell state" ) + 
    ylim( 0, ylim_max ) + 
    coord_flip() +
    annotate( geom = "text", label = text_die, x = 0.75 * n_sim, y = q95 - 10 ) +
    geom_hline( yintercept = q95 ) + 
    annotate( geom = "text", x = 0.9 * n_sim, y = q95 - 10, 
              label = sprintf( "q95=%s", q95 %>% as.character() ) ) + 
    scale_fill_manual( values = pal_vector,
                       guide = guide_legend( reverse = TRUE ),
                       breaks = c( "o000", "o001", "o002" ),
                       labels = c( "naive", "activated", "CD69+" ) ) -> 
    ggfig_steps; ggfig_steps
  
  if ( export ) {
    pdf_name <- sprintf( "%s/%s_%s_%isims_dwell_steps.pdf", figs_dwell_dir, 
                         fig_fname_start, tolower( m$tissue ), n_sim )
    ggexport_my( ggfig_steps, filename = pdf_name,
                 width = 21 * cm_to_in, height = 6 * cm_to_in )
  }

  return( ggfig_steps )
}
