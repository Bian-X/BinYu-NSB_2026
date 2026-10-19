# Validate identities independently of report column order and filename indices.
validate_region_samples <- function(columns, region, groups=paste0('G',1:6)) {
  clean <- sub('\\.raw\\.PG\\.Quantity$','',sub('^\\[[0-9]+\\] +','',columns))
  pattern <- paste0('^Astral_MSP2600441_[0-9]+_(G[1-6])_([1-4])_',region,'(_[0-9]{14})?$')
  if(any(!grepl(pattern,clean)) || anyDuplicated(clean)) stop('Invalid or duplicate sample identifiers in ',region)
  identities <- sub(pattern,'\\1_\\2',clean)
  expected <- as.vector(outer(groups,1:4,paste,sep='_'))
  if(anyDuplicated(identities) || length(identities)!=length(expected) || !setequal(identities,expected)) {
    stop('Expected unique groups x replicates 1-4 for ',region)
  }
  invisible(list(sample_id=clean,identity=identities))
}
