# Gene-level ORA: exclude multi-gene protein groups by default; expanded is sensitivity.
prepare_go_gene_sets <- function(result, background_rows, foreground_rows,
                                 mapping_policy=c('unique','expanded')) {
  mapping_policy <- match.arg(mapping_policy)
  symbol_map <- map_symbol_to_entrez(split_symbols(result$Genes[background_rows]))
  records <- lapply(seq_len(nrow(result)),function(i) {
    symbols <- split_symbols(result$Genes[i])
    ids <- unique(symbol_map$ENTREZID[symbol_map$SYMBOL %in% symbols])
    included <- background_rows[i] && if(mapping_policy=='unique') length(ids)==1L else length(ids)>0L
    data.frame(ProteinAccessions=result$ProteinAccessions[i],Genes=result$Genes[i],
               Tested=background_rows[i],SelectedForeground=foreground_rows[i],
               MappedEntrezIDs=paste(ids,collapse=';'),MappedGeneCount=length(ids),
               IncludedByMapping=included,MappingPolicy=mapping_policy)
  })
  audit <- do.call(rbind,records)
  build <- function(selected) {
    symbols <- split_symbols(result$Genes[selected & audit$IncludedByMapping])
    unique(symbol_map[symbol_map$SYMBOL %in% symbols,c('SYMBOL','ENTREZID'),drop=FALSE])
  }
  background <- build(background_rows)
  foreground <- build(foreground_rows)
  foreground <- foreground[foreground$ENTREZID %in% background$ENTREZID,,drop=FALSE]
  list(background=background,foreground=foreground,audit=audit,
       background_groups=sum(background_rows & audit$IncludedByMapping),
       foreground_groups=sum(foreground_rows & audit$IncludedByMapping))
}
