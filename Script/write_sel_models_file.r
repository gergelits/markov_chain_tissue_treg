# write_sel_models_file.r

write_sel_models_file <- function( 
    celltype, ps, model.name, model.ver, sel_models_file_name )
{
  model.dir <- sprintf( "%s_%i%s", model.name, 
                        mcmc.chain.n * ( mcmc.iter.n - mcmc.warmup.n ),
                        model.ver )
  pmc <- get_parabiosis_const()
  
  dSelected_models <- tibble( 
    celltype = celltype,
    tissue = pmc$dTissueAllOrderedGroup$tissue.all.ordered, 
    model.id = NA, 
    model.dir.name = NA )
  
  for ( tissue in pmc$dTissueAllOrderedGroup$tissue.all.ordered ) {
    output.dir <- sprintf( 
      "%s/%s/%s/%s", ps$ANALYSIS_PATH, celltype, tissue, model.dir )
    if (
      # this celltype-tissue-model is ready to be plotted
      file.exists(
        sprintf( "%s/parabio_fit%i.rda", output.dir, mcmc.iter.n ) ) &
      file.exists( 
        sprintf( "%s/parabio_fit_HDI.csv", output.dir ) )
    ) {
      dSelected_models$model.id[ dSelected_models$tissue == tissue ] <- 
        str_match( model.dir, "markov_model__(.*?)_")[ , 2 ]
      dSelected_models$model.dir.name[ dSelected_models$tissue == tissue ] <- 
        model.dir
  }
  dSelected_models %>% write.table( 
    sprintf( "%s/%s", ps$ANALYSIS_PATH, sel_models_file_name ), 
    sep = ",", row.names = FALSE )
  }
}
 