# visualize_results.r

visualize_results <- function( 
    m_tisNA = NULL,
    sel_models_file_name )
{
  plot_figures_3_5(
    m_tisNA = m_tisNA,
    sel_models_file_name = sel_models_file_name )
  
  # rare issue with annotation_logticks() at first run
  try(
    plot_figures_2_4(
      m_tisNA = m_tisNA,
      sel_models_file_name = sel_models_file_name ) ) 
  
  try(
    plot_figures_2_4(
      m_tisNA = m_tisNA,
      sel_models_file_name = sel_models_file_name,
      qij_I = 2L ) )
  
  write_dwell_time(
    m_tisNA = m_tisNA,
    sel_models_file_name = sel_models_file_name )
  
  plot_dt_distributions( m_tisNA = m_tisNA, cellstate = "cd69p" )
  if ( ! CALC_TTP ) {
    plot_dt_distributions( m_tisNA = m_tisNA, cellstate = "naive" )
    plot_dt_distributions( m_tisNA = m_tisNA, cellstate = "activ" )
  }
  
  plot_aggr_tissuegroup_cellstate_barplot( m_tisNA = m_tisNA )
  
  # For 9-states flow figure
  write_aggregated_group_flows( 
    m_tisNA = m_tisNA,
    sel_models_file_name = sel_models_file_name,
    write_csv = TRUE )
  
  get_aggr_tissuegroup_cellstate_areas( 
    m_tisNA = m_tisNA,
    write_csv = TRUE )
  
  get_aggr_total_counts_areas( 
    m_tisNA = m_tisNA,
    write_csv = TRUE )
  
  merge_fig5_subfigures( 
    m_tisNA = m_tisNA,
    sel_models_file_name = sel_models_file_name )
  
  plot_all_tissue_trajectories(
    m_tisNA = m_tisNA,
    sel_models_file_name = sel_models_file_name )
}
