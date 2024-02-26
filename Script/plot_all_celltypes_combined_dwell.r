# plot_all_celltypes_combined_dwell.r

if ( FALSE ) {
  source( "universal_script_setup.r" )
}


d_res <- tibble(
  celltype = c( "Tconv", "CD8", "Treg", "NK", "B" ),
  mcmc = c( "i1001", "i1001", "i1001", "i1001", "i1001" ),
  # mcmc = c( "i10020", "i10020", "i5004", "i10020", "i10020" ),
  # model_id = c( "m0025", "m0025", "m0001", "m0025", "m0025" )
  model_id = c( "m0057", "m0057", "m0057", "m0057", "m0057" )
)


read_medians <- function( d_res )
{
  d_medians <- tibble( tissue = NULL, median = NULL )
  for ( i in 1 : nrow( d_res ) ) {
    celltype <- d_res$celltype[ i ]
    mcmc <- d_res$mcmc[ i ]
    model_id <- d_res$model_id[ i ]
    
    for ( tissue_j in get_TISSUES( celltype = celltype ) )
    {
      stan_pars_v <- "v01"; ps <- PS
      m <- model( ps = ps,
                  celltype = celltype,
                  tissue = tissue_j, 
                  mcmc_pars = get_mcmc_pars( mcmc_pars_v = mcmc, ps = ps ),
                  stan_pars_v = stan_pars_v,
                  model_name = sprintf( "model_%s", model_id ),
                  model_ver = sprintf( "%s_t%s", stan_pars_v, N_TISSUES = 1 ),
                  max_lp_chain = 1 )
      
      fig_fname_start <- get_fig_fname_start( m = m )
      figs_dwell_dir <- get_path( "figs_m_dwell", m = m )
      
      i_j_file <- sprintf( "%s/%s_%s_%ssims_median_dwell.csv", figs_dwell_dir, 
        fig_fname_start, tolower( m$tissue ), 10000 ) 
      if ( file.exists( i_j_file ) ) {
        read_csv( i_j_file ) %>% 
          mutate( celltype = celltype ) %>% 
          bind_rows( d_medians, . ) ->
          d_medians
      }
    }
  }
  return( d_medians )
}

d_medians <- read_medians( d_res = d_res )

plot_all_comb_dwell <- function( d_medians )
{
  d_medians %>% 
    left_join( ., get_pmc()$dTissueAllOrderedGroup.f %>% 
                 dplyr::select( tissue.all, f.tissue, f.tissue.group ) %>% 
                 mutate( lower.tissue.all = tolower( tissue.all ) ), 
               by = c( "tissue" = "lower.tissue.all" ) ) ->
    d_medians_groups
  
  xintercept1 <- d_medians_groups$f.tissue[ 
    d_medians_groups$f.tissue.group == "BoneMarrow" ] %>% unique() %>% length()
  xintercept2 <- d_medians_groups$f.tissue[ 
    d_medians_groups$f.tissue.group == "Lymphoid" ] %>% unique() %>% length()
  xintercept3 <- d_medians_groups$f.tissue[ 
    d_medians_groups$f.tissue.group == "Non-lymphoid" ] %>% unique() %>% length()
  
  
  cm_to_in <- 1 / 2.54   
  f <- cm_to_in * 1.7
  
  d_medians_groups %>%  
    mutate( median2 = ifelse( median < 1.1, 1.1, median ) ) %>%
    ggplot( aes( x = f.tissue, y = median, fill = celltype ) ) +
    geom_bar( stat = "identity", position = "dodge" ) +
    geom_vline( xintercept = c( xintercept1 + 0.5, 
                                xintercept1 + xintercept2 + 0.5,
                                xintercept1 + xintercept2 + xintercept3 + 0.5 ),
                linetype = "dashed" ) +
    labs( y = "Median Dwell Time (days)", fill = "Cell type" ) +    
    theme_classic() +
    theme( plot.title = element_text( size = 12 * f, hjust = 0 ), 
           axis.title.x = element_blank(),
           axis.title.y = element_text( size = 12 * f, hjust = 0.5 ),  
           axis.text.x = element_text( angle = 45, hjust = 1, size = 12 * f ),
           axis.text.y = element_text( size = 12 * f ),  
           panel.border = element_blank(),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           legend.text = element_text( size = 12 * f ),
           legend.title = element_text( size = 12 * f ),
           strip.text.x = element_text( size = 12 * f ) ) -> gg_all_comb_dwell
  return( gg_all_comb_dwell )
}
gg_all <- plot_all_comb_dwell( d_medians = d_medians )

gg_all

# d_medians %>% View

